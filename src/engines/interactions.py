"""
Drug-Drug Interaction Checker for NeuroBotanica Platform
AI-Enhanced Drug Interaction Analysis with Cultural Preservation

This module checks for CYP450 enzyme interactions between cannabinoids/compounds and other drugs.
Integrates with OmniPath for TK consent and HIPAA-compliant logging.

Performance Target: <1s per check, 95% accuracy on known interactions.
"""

import sqlite3
from typing import Dict, List, Optional
from datetime import datetime
import os

# Cloudflare D1 connection
# For Workers: import { D1Database } from '@cloudflare/workers-types'; use env.DB
DB_PATH = "neurobotanica.db"  # Local for dev; use D1 endpoint in prod

class DrugInteractionChecker:
    def __init__(self, db=None):
        self.db = db  # Cloudflare D1 database instance (passed from Worker)
        self.cache = {}  # KV caching placeholder

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

    def check_interactions(self, compound_ids: List[str], drug_names: List[str], customer_tier: str = "computational_only") -> Dict:
        """
        Check for drug interactions.
        
        Args:
            compound_ids: List of NeuroBotanica compound IDs (e.g., ['cbd', 'thc']).
            drug_names: List of drug names (e.g., ['warfarin', 'clopidogrel']).
            customer_tier: 'computational_only' or 'tk_enhanced' for TK access.
        
        Returns:
            Dict with warnings, severity, and recommendations.
        """
        if not drug_names:
            return {"warnings": [], "total_warnings": 0, "processing_time_ms": 0}
        interactions = []
        for compound_id in compound_ids:
            for drug_name in drug_names:
                interaction = self._query_interaction(compound_id, drug_name)
                if interaction:
                    # TK Check: If TK-enhanced tier and interaction involves TK data, verify consent
                    if customer_tier == "tk_enhanced" and interaction.get("requires_consent_check"):
                        consent_valid = self._verify_consent(interaction.get("consent_id"))
                        if not consent_valid:
                            continue  # Skip if consent revoked
                    interactions.append({
                        "compound": compound_id,
                        "drug": drug_name,
                        "severity": interaction["severity_level"],
                        "effect": interaction["clinical_effect"],
                        "recommendation": interaction["adjustment_recommendation"],
                        "evidence_level": interaction["evidence_level"],
                        "citations": interaction["citations"]
                    })
        # Log to audit for HIPAA compliance
        self._log_audit("drug_interaction_check", {"compound_ids": compound_ids, "drug_names": drug_names}, len(interactions) > 0)
        return {
            "warnings": interactions,
            "total_warnings": len(interactions),
            "processing_time_ms": 0  # Placeholder; measure in prod
        }

    def _query_interaction(self, compound_id: str, drug_name: str) -> Optional[Dict]:
        """Query D1 for interaction data with KV caching."""
        key = f"{compound_id}_{drug_name}"
        if key in self.cache:
            return self.cache[key]
        
        query = """
        SELECT severity_level, clinical_effect, adjustment_recommendation, evidence_level, citations, requires_consent_check, consent_id
        FROM neurobotanica_drug_interactions
        WHERE compound_id = ? AND drug_name = ?
        """
        results = self._execute_query(query, (compound_id, drug_name))
        if results:
            if self.db:
                # D1 returns list of objects
                row = results[0]
                result = {
                    "severity_level": row.severity_level,
                    "clinical_effect": row.clinical_effect,
                    "adjustment_recommendation": row.adjustment_recommendation,
                    "evidence_level": row.evidence_level,
                    "citations": row.citations,
                    "requires_consent_check": row.requires_consent_check,
                    "consent_id": row.consent_id
                }
            else:
                # PostgreSQL or SQLite return dict-like
                row = results[0]
                result = {
                    "severity_level": row["severity_level"],
                    "clinical_effect": row["clinical_effect"],
                    "adjustment_recommendation": row["adjustment_recommendation"],
                    "evidence_level": row["evidence_level"],
                    "citations": row["citations"],
                    "requires_consent_check": row["requires_consent_check"],
                    "consent_id": row["consent_id"]
                }
        else:
            result = None
        self.cache[key] = result
        return result

    def _verify_consent(self, consent_id: str) -> bool:
        """Verify TK consent via OmniPath. Triggers policy rules on failure."""
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
        # In prod, implement halt logic (e.g., delay response)
        print(f"Policy violation: {violation_type} - Global halt triggered")  # Placeholder

    def _log_audit(self, event_type: str, event_data: Dict, result: bool):
        """Log to omnipath_audit_log for HIPAA with community attribution."""
        # Add community attribution for TK events
        if "consent_id" in str(event_data):
            attribution = self._get_community_attribution(event_data.get("consent_id"))
            event_data["community_attribution"] = attribution
        
        query = """
        INSERT INTO omnipath_audit_log (event_type, application, event_data, result, hipaa_deidentified)
        VALUES (?, 'neurobotanica', ?, ?, 1)
        """
        if self.db:
            # D1
            stmt = self.db.prepare(query).bind(event_type, str(event_data), "approved" if result else "denied", 1)
            stmt.run()
        else:
            conn = self._get_connection()
            if 'psycopg2' in str(type(conn)):
                query = query.replace('?', '%s')
                cursor = conn.cursor()
                cursor.execute(query, (event_type, str(event_data), "approved" if result else "denied"))
                conn.commit()
                cursor.close()
                conn.close()
            else:
                conn.execute(query, (event_type, str(event_data), "approved" if result else "denied"))
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

# Example usage (for testing)
if __name__ == "__main__":
    checker = DrugInteractionChecker()
    result = checker.check_interactions(["cbd"], ["warfarin"], "tk_enhanced")
    print(result)
    checker.close()