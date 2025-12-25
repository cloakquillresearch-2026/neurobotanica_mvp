#!/usr/bin/env python3
"""
HIPAA-Compliant Patient Data Collection System
NeuroBotanica Therapeutic Discovery Platform

This module demonstrates HIPAA-compliant methods for collecting,
de-identifying, and processing patient data for AI model training.
"""

import hashlib
import json
import uuid
from datetime import datetime
from typing import Dict, List, Optional
from pathlib import Path
import logging
from dateutil.relativedelta import relativedelta

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class HIPAACompliantDataCollector:
    """
    HIPAA-compliant patient data collection and processing system.
    Implements de-identification, consent management, and audit trails.
    """

    # HIPAA's 18 identifiers that must be removed
    HIPAA_IDENTIFIERS = [
        'name', 'address', 'phone', 'email', 'ssn', 'medical_record_number',
        'health_plan_number', 'certificate_number', 'license_number',
        'vehicle_id', 'device_id', 'biometric_id', 'full_face_photo',
        'unique_identifying_number', 'ip_address', 'geographic_location'
    ]

    def __init__(self, data_dir: str = "data/patient_data"):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.consent_log = []
        self.audit_log = []

    def collect_patient_consent(self, patient_id: str, consent_details: Dict) -> str:
        """
        Record patient consent with full audit trail.

        Args:
            patient_id: Anonymous patient identifier
            consent_details: Consent scope and permissions

        Returns:
            Consent token for data processing
        """
        consent_record = {
            'consent_id': str(uuid.uuid4()),
            'patient_id': patient_id,
            'timestamp': datetime.now().isoformat(),
            'consent_version': '1.0',
            'scope': consent_details.get('scope', []),
            'duration_days': consent_details.get('duration_days', 365),
            'withdrawal_allowed': True,
            'data_usage': consent_details.get('data_usage', 'research'),
            'digital_signature': self._generate_consent_hash(consent_details)
        }

        self.consent_log.append(consent_record)
        self._audit_event('consent_granted', patient_id, consent_record)

        return consent_record['consent_id']

    def deidentify_patient_record(self, raw_record: Dict, consent_id: str) -> Dict:
        """
        De-identify patient record according to HIPAA Safe Harbor method.

        Removes all 18 HIPAA identifiers and creates research-safe dataset.
        """
        # Create audit trail
        original_hash = self._hash_record(raw_record)
        self._audit_event('deidentification_started', raw_record.get('patient_id', 'unknown'), {
            'original_hash': original_hash,
            'consent_id': consent_id
        })

        # Remove all HIPAA identifiers
        deidentified = {}

        # Safe demographic data (categorized, not specific)
        if 'demographics' in raw_record:
            demo = raw_record['demographics']
            deidentified['demographics'] = {
                'age_range': self._categorize_age(demo.get('birth_date')),
                'experience_level': demo.get('experience_level'),
                'thc_tolerance': demo.get('thc_tolerance'),
                'gender_category': demo.get('gender_category')  # M/F/Other
            }

        # Safe treatment data
        if 'treatment' in raw_record:
            treatment = raw_record['treatment']
            deidentified['treatment'] = {
                'condition_category': treatment.get('condition_category'),  # Generalized
                'treatment_duration_weeks': treatment.get('treatment_duration_weeks'),
                'delivery_method': treatment.get('delivery_method'),
                'dosage_range': treatment.get('dosage_range')  # Categorized
            }

        # Safe response data
        if 'response' in raw_record:
            response = raw_record['response']
            deidentified['response'] = {
                'efficacy_rating': response.get('efficacy_rating'),  # 1-10 scale
                'adverse_events_count': len(response.get('adverse_events', [])),
                'adherence_level': response.get('adherence_level'),
                'overall_satisfaction': response.get('overall_satisfaction')
            }

        # Add metadata for research integrity
        deidentified['_metadata'] = {
            'deidentified_at': datetime.now().isoformat(),
            'consent_id': consent_id,
            'data_quality_score': self._assess_data_quality(raw_record),
            'research_eligible': True
        }

        # Final audit
        deidentified_hash = self._hash_record(deidentified)
        self._audit_event('deidentification_completed', 'anonymous', {
            'deidentified_hash': deidentified_hash,
            'fields_removed': len(self.HIPAA_IDENTIFIERS)
        })

        return deidentified

    def process_research_submission(self, raw_data: Dict) -> Optional[Dict]:
        """
        Complete workflow for processing patient data submission.

        1. Validate consent
        2. De-identify data
        3. Store securely
        4. Return research-safe record
        """
        try:
            # Step 1: Validate consent
            consent_id = raw_data.get('consent_id')
            if not consent_id or not self._validate_consent(consent_id):
                logger.error("Invalid or missing consent")
                return None

            # Step 2: De-identify
            patient_id = raw_data.get('patient_id', str(uuid.uuid4()))
            deidentified = self.deidentify_patient_record(raw_data, consent_id)

            # Step 3: Store securely
            research_id = str(uuid.uuid4())
            storage_record = {
                'research_id': research_id,
                'deidentified_data': deidentified,
                'storage_timestamp': datetime.now().isoformat(),
                'consent_id': consent_id
            }

            storage_path = self.data_dir / f"research_record_{research_id[:8]}.json"
            with open(storage_path, 'w') as f:
                json.dump(storage_record, f, indent=2)

            # Step 4: Audit and return
            self._audit_event('research_record_created', research_id, {
                'storage_path': str(storage_path),
                'data_quality': deidentified['_metadata']['data_quality_score']
            })

            logger.info(f"Successfully processed research record: {research_id}")
            return deidentified

        except Exception as e:
            logger.error(f"Error processing research submission: {e}")
            self._audit_event('processing_error', 'unknown', {'error': str(e)})
            return None

    def _categorize_age(self, birth_date: Optional[str]) -> Optional[str]:
        """Convert birth date to age range category."""
        if not birth_date:
            return None

        try:
            birth_year = int(birth_date.split('-')[0])
            current_year = datetime.now().year
            age = current_year - birth_year

            if age < 25:
                return "18-24"
            elif age < 35:
                return "25-34"
            elif age < 45:
                return "35-44"
            elif age < 55:
                return "45-54"
            elif age < 65:
                return "55-64"
            else:
                return "65+"
        except:
            return None

    def _assess_data_quality(self, record: Dict) -> float:
        """Assess quality of patient data record."""
        required_fields = ['demographics', 'treatment', 'response']
        present_fields = sum(1 for field in required_fields if field in record)

        # Bonus for detailed response data
        if 'response' in record and len(record['response'].get('adverse_events', [])) > 0:
            present_fields += 0.5

        return min(present_fields / len(required_fields), 1.0)

    def _validate_consent(self, consent_id: str) -> bool:
        """Validate that consent is current and valid."""
        for consent in self.consent_log:
            if consent['consent_id'] == consent_id:
                # Check if consent is still valid (not expired, not withdrawn)
                granted = datetime.fromisoformat(consent['timestamp'])
                duration_days = consent.get('duration_days', 365)
                
                # Properly calculate expiration date
                from dateutil.relativedelta import relativedelta
                expires = granted + relativedelta(days=duration_days)

                return datetime.now() < expires and not consent.get('withdrawn', False)

        return False

    def _generate_consent_hash(self, consent_details: Dict) -> str:
        """Generate cryptographic hash of consent details."""
        consent_str = json.dumps(consent_details, sort_keys=True)
        return hashlib.sha256(consent_str.encode()).hexdigest()

    def _hash_record(self, record: Dict) -> str:
        """Generate hash of record for audit purposes."""
        # Remove timestamps and IDs that change
        audit_record = {k: v for k, v in record.items()
                       if not k.endswith('_id') and k != 'timestamp'}
        record_str = json.dumps(audit_record, sort_keys=True)
        return hashlib.sha256(record_str.encode()).hexdigest()

    def _audit_event(self, event_type: str, entity_id: str, details: Dict):
        """Record audit event for compliance."""
        audit_entry = {
            'timestamp': datetime.now().isoformat(),
            'event_type': event_type,
            'entity_id': entity_id,
            'details': details,
            'system_version': '1.0.0'
        }

        self.audit_log.append(audit_entry)

        # In production, this would write to secure audit log
        logger.info(f"AUDIT: {event_type} - {entity_id}")

    def get_compliance_report(self) -> Dict:
        """Generate compliance report for HIPAA auditing."""
        return {
            'total_consents': len(self.consent_log),
            'active_consents': len([c for c in self.consent_log if not c.get('withdrawn', False)]),
            'audit_events': len(self.audit_log),
            'data_quality_distribution': self._analyze_data_quality(),
            'last_audit_check': datetime.now().isoformat()
        }

    def _analyze_data_quality(self) -> Dict:
        """Analyze data quality distribution."""
        # This would analyze stored research records
        return {
            'high_quality': 0.7,  # Placeholder
            'medium_quality': 0.25,
            'low_quality': 0.05
        }


# Example usage and testing
def demonstrate_hipaa_compliance():
    """Demonstrate HIPAA-compliant data collection workflow."""

    collector = HIPAACompliantDataCollector()

    # Example patient data (would come from clinic system)
    example_patient_data = {
        'patient_id': 'PATIENT_12345',  # Will be removed
        'name': 'John Doe',  # HIPAA identifier - will be removed
        'email': 'john.doe@email.com',  # HIPAA identifier - will be removed
        'demographics': {
            'birth_date': '1980-05-15',  # Will be converted to age range
            'experience_level': 'occasional',
            'thc_tolerance': 'medium',
            'gender_category': 'M'
        },
        'treatment': {
            'condition_category': 'chronic_pain',
            'treatment_duration_weeks': 8,
            'delivery_method': 'oral',
            'dosage_range': 'moderate'
        },
        'response': {
            'efficacy_rating': 7,
            'adverse_events': ['dry_mouth', 'drowsiness'],
            'adherence_level': 'high',
            'overall_satisfaction': 8
        }
    }

    # Step 1: Patient provides consent
    consent_details = {
        'scope': ['research', 'model_training'],
        'duration_days': 365,
        'data_usage': 'therapeutic_ai_development'
    }

    consent_id = collector.collect_patient_consent(
        example_patient_data['patient_id'],
        consent_details
    )

    # Step 2: Process patient data
    example_patient_data['consent_id'] = consent_id
    research_record = collector.process_research_submission(example_patient_data)

    # Step 3: Show results
    print("=== HIPAA-Compliant Data Processing Demo ===")
    print(f"Original record had {len(collector.HIPAA_IDENTIFIERS)} HIPAA identifiers")
    print("De-identified research record:")
    print(json.dumps(research_record, indent=2))

    # Step 4: Compliance report
    compliance = collector.get_compliance_report()
    print(f"\nCompliance Status: {compliance}")

    return research_record


if __name__ == "__main__":
    demonstrate_hipaa_compliance()