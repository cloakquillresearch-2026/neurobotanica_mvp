"""
Demographic Bias Correction for NeuroBotanica Platform
AI-Enhanced Fairness Adjustments with Cultural Preservation

This module applies demographic corrections (age, gender, genetics) to dosing and profiles.
Integrates with OmniPath for TK consent and HIPAA-compliant logging.

Performance Target: <1s per correction, <5% bias variance.
"""

import sqlite3
from typing import Dict, List, Optional
from datetime import datetime
import os

# Cloudflare D1 connection
DB_PATH = "neurobotanica.db"  # Local for dev

class DemographicBiasCorrection:
    def __init__(self, db=None):
        self.db = db  # Cloudflare D1 database instance
        self.cache = {}  # KV caching for demographic factors

    def _get_connection(self):
        """Get database connection - PostgreSQL for Railway, SQLite for local."""
        db_url = os.getenv('DATABASE_URL')
        if db_url and not db_url.startswith('sqlite'):
            # Railway PostgreSQL
            import psycopg2
            return psycopg2.connect(db_url)
        else:
            # Local SQLite
            conn = sqlite3.connect(DB_PATH)
            conn.row_factory = sqlite3.Row
            return conn

    def _execute_query(self, query: str, params: tuple = ()):
        """Execute query on D1, PostgreSQL, or SQLite."""
        if self.db:
            # Cloudflare D1
            stmt = self.db.prepare(query)
            if params:
                stmt = stmt.bind(*params)
            return stmt.all()
        else:
            # Railway PostgreSQL or local SQLite
            conn = self._get_connection()
            if 'psycopg2' in str(type(conn)):
                # PostgreSQL - replace ? with %s
                query = query.replace('?', '%s')
                cursor = conn.cursor()
                cursor.execute(query, params)
                if query.strip().upper().startswith('SELECT'):
                    columns = [desc[0] for desc in cursor.description] if cursor.description else []
                    results = cursor.fetchall()
                    # Convert to dict-like for compatibility
                    dict_results = []
                    for row in results:
                        dict_results.append(dict(zip(columns, row)))
                    cursor.close()
                    conn.close()
                    return dict_results
                else:
                    conn.commit()
                    cursor.close()
                    conn.close()
                    return None
            else:
                # SQLite
                cursor = conn.cursor()
                cursor.execute(query, params)
                if query.strip().upper().startswith('SELECT'):
                    results = cursor.fetchall()
                    conn.close()
                    return results
                else:
                    conn.commit()
                    conn.close()
                    return None

    def apply_corrections(self, base_dose_mg: float, compound_id: str, demographics: Dict, customer_tier: str = "computational_only") -> Dict:
        """
        Apply demographic corrections to base dose.
        
        Args:
            base_dose_mg: Initial dose (e.g., 10.0).
            compound_id: NeuroBotanica compound ID.
            demographics: Dict with 'age', 'gender', 'genetics' (e.g., {'age': 65, 'gender': 'female'}).
            customer_tier: For TK access.
        
        Returns:
            Dict with adjusted dose, factors applied, and bias metrics.
        """
        adjustments = self._query_adjustments(compound_id, demographics)
        adjusted_dose = base_dose_mg * adjustments.get("cyp450_adjustment", 1.0) * adjustments.get("dosing_adjustment", 1.0)
        
        # TK Check: If TK-enhanced and adjustments involve TK data, verify consent
        if customer_tier == "tk_enhanced" and adjustments.get("traditional_knowledge_source"):
            consent_valid = self._verify_consent(adjustments.get("consent_id"))
            if not consent_valid:
                adjusted_dose = base_dose_mg  # Fallback to base
        
        # Calculate bias variance (enhanced with population data)
        bias_variance = self._calculate_bias_variance(adjusted_dose, demographics)
        
        # Log to audit with community attribution
        self._log_audit("bias_correction_applied", {"compound_id": compound_id, "demographics": demographics, "consent_id": adjustments.get("consent_id")}, bias_variance < 0.05)
        
        return {
            "adjusted_dose_mg": round(adjusted_dose, 2),
            "factors_applied": adjustments,
            "bias_variance": bias_variance,
            "processing_time_ms": 0  # Measure in prod
        }

    def _query_adjustments(self, compound_id: str, demographics: Dict) -> Dict:
        """Query demographic factors with KV caching."""
        key = f"{compound_id}_{demographics.get('age', 'unknown')}_{demographics.get('gender', 'unknown')}"
        if key in self.cache:
            return self.cache[key]
        
        query = """
        SELECT cyp450_adjustment, dosing_adjustment
        FROM neurobotanica_demographic_factors
        WHERE demographic_category = ? AND demographic_value = ? AND applicable_compounds LIKE ?
        """
        category = "age" if "age" in demographics else "gender"
        value = demographics.get(category, "unknown")
        results = self._execute_query(query, (category, value, f"%{compound_id}%"))
        if results:
            if self.db:
                result = {"cyp450_adjustment": results[0].cyp450_adjustment, "dosing_adjustment": results[0].dosing_adjustment}
            else:
                result = {"cyp450_adjustment": results[0][0], "dosing_adjustment": results[0][1]}
        else:
            result = {"cyp450_adjustment": 1.0, "dosing_adjustment": 1.0}
        self.cache[key] = result
        return result

    def _verify_consent(self, consent_id: str) -> bool:
        """Verify TK consent. Triggers policy rules on failure."""
        query = "SELECT consent_status FROM omnipath_consent_artifacts WHERE consent_id = ?"
        results = self._execute_query(query, (consent_id,))
        if results:
            if self.db:
                status = results[0].consent_status
            else:
                status = results[0]["consent_status"]
        else:
            status = None
        if status != "active":
            # Trigger policy halt
            self._trigger_policy_halt("consent_revoked", {"consent_id": consent_id})
            return False
        return True

    def _trigger_policy_halt(self, violation_type: str, details: Dict):
        """Trigger 2-second global halt for policy violations."""
        # Log violation
        self._log_audit("policy_violation", details, False)
        # In prod, implement halt logic
        print(f"Policy violation: {violation_type} - Global halt triggered")  # Placeholder

    def _calculate_bias_variance(self, adjusted_dose: float, demographics: Dict) -> float:
        """Enhanced bias variance using population averages."""
        # Mock population data (replace with real in prod)
        population_averages = {
            "age_65_female": 7.2,  # mg (adjusted for 65+ female)
            "age_30_male": 12.0,
            "default": 10.0
        }
        key = f"age_{demographics.get('age', 'unknown')}_{demographics.get('gender', 'unknown')}"
        avg_dose = population_averages.get(key, population_averages["default"])
        return abs(adjusted_dose - avg_dose) / avg_dose

    def _log_audit(self, event_type: str, event_data: Dict, result: bool):
        """Log to omnipath_audit_log with community attribution."""
        # Add community attribution for TK events
        if "consent_id" in str(event_data):
            attribution = self._get_community_attribution(event_data.get("consent_id"))
            event_data["community_attribution"] = attribution
        
        query = """
        INSERT INTO omnipath_audit_log (event_type, application, event_data, result, hipaa_deidentified)
        VALUES (?, 'neurobotanica', ?, ?, 1)
        """
        if self.db:
            stmt = self.db.prepare(query).bind(event_type, str(event_data), "fair" if result else "biased", 1)
            stmt.run()
        else:
            conn = self._get_connection()
            if 'psycopg2' in str(type(conn)):
                query = query.replace('?', '%s')
                cursor = conn.cursor()
                cursor.execute(query, (event_type, str(event_data), "fair" if result else "biased"))
                conn.commit()
                cursor.close()
                conn.close()
            else:
                conn.execute(query, (event_type, str(event_data), "fair" if result else "biased"))
                conn.commit()
                conn.close()

    def _get_community_attribution(self, consent_id: str) -> str:
        """Fetch community attribution for TK data."""
        query = "SELECT tk_labels FROM omnipath_consent_artifacts WHERE consent_id = ?"
        results = self._execute_query(query, (consent_id,))
        if results:
            if self.db:
                return results[0].tk_labels
            else:
                return results[0]["tk_labels"]
        return ""

    def close(self):
        pass  # No persistent connection to close

# Example usage
if __name__ == "__main__":
    corrector = DemographicBiasCorrection()
    result = corrector.apply_corrections(10.0, "cbd", {"age": 65, "gender": "female"}, "tk_enhanced")
    print(result)
    corrector.close()