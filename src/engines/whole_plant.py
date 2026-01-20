"""
Whole-Plant Analysis Engine for NeuroBotanica Platform
Aggregates cannabinoid profiles, computes ratios, integrates TK sources.

This module analyzes whole-plant formulations for optimal cannabinoid ratios.
Integrates with OmniPath for TK consent and HIPAA-compliant logging.

Performance Target: <2s analysis, 95% match to data.
"""

import sqlite3
from typing import Dict, List, Optional
from datetime import datetime
import os

DB_PATH = "neurobotanica.db"

class WholePlantAnalysisEngine:
    def __init__(self, db=None):
        self.db = db
        self.cache = {}

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

    def analyze_whole_plant(self, plant_id: str, target_condition: str, customer_tier: str = "computational_only") -> Dict:
        """
        Analyze whole-plant profile for optimal formulation.
        
        Args:
            plant_id: Plant compound ID (e.g., 'cbd').
            target_condition: Condition to treat.
            customer_tier: For TK access.
        
        Returns:
            Dict with ratios, mechanisms, and recommendations.
        """
        # Get plant compounds
        compounds = self._get_plant_compounds(plant_id)
        if not compounds:
            return {"error": "Plant not found"}
        
        # Compute ratios
        ratios = self._compute_cannabinoid_ratios(compounds)
        
        # TK enhancement
        tk_data = None
        if customer_tier == "tk_enhanced":
            tk_data = self._query_tk_formulations(plant_id, target_condition)
            if tk_data and not self._verify_consent(tk_data.get("consent_ids")):
                tk_data = None  # Fallback
        
        # Synergy integration
        synergy_score = self._integrate_synergy(compounds, target_condition)
        
        # Build formulation
        formulation = {
            "plant_id": plant_id,
            "target_condition": target_condition,
            "cannabinoid_ratios": ratios,
            "predicted_mechanisms": ["anti_inflammatory", "pain_relief"],
            "synergy_score": synergy_score,
            "tk_enhanced": tk_data is not None,
            "recommended_dose_mg": 25.0,  # Mock
            "processing_time_ms": 0
        }
        
        # Log formulation
        self._log_formulation(formulation, customer_tier)
        
        return formulation

    def _get_plant_compounds(self, plant_id: str) -> List[Dict]:
        """Get compounds for the plant."""
        query = "SELECT compound_id, compound_name, molecular_weight FROM neurobotanica_compounds WHERE compound_id LIKE ?"
        results = self._execute_query(query, (f"%{plant_id}%",))
        if self.db:
            return [{"compound_id": r.compound_id, "compound_name": r.compound_name, "molecular_weight": r.molecular_weight} for r in results]
        else:
            return [dict(r) for r in results]

    def _compute_cannabinoid_ratios(self, compounds: List[Dict]) -> Dict:
        """Compute optimal ratios (mock: equal distribution)."""
        total = len(compounds)
        return {c["compound_id"]: 1.0 / total for c in compounds}

    def _query_tk_formulations(self, plant_id: str, condition: str) -> Optional[Dict]:
        """Query TK-based formulations."""
        query = """
        SELECT consent_ids FROM neurobotanica_formulations
        WHERE ingredients LIKE ? AND target_condition = ? AND requires_consent = 1
        """
        results = self._execute_query(query, (f"%{plant_id}%", condition))
        if results:
            return dict(results[0])
        return None

    def _verify_consent(self, consent_ids: str) -> bool:
        """Verify TK consents."""
        if not consent_ids:
            return False
        ids = consent_ids.split(",")
        for cid in ids:
            query = "SELECT consent_status FROM omnipath_consent_artifacts WHERE consent_id = ?"
            results = self._execute_query(query, (cid.strip(),))
            if results:
                if self.db:
                    status = results[0].consent_status
                else:
                    status = results[0]["consent_status"]
                if status != "active":
                    self._trigger_policy_halt("tk_consent_denied", {"consent_ids": consent_ids})
                    return False
            else:
                return False
        return True

    def _integrate_synergy(self, compounds: List[Dict], condition: str) -> float:
        """Integrate synergy scores."""
        # Mock: average synergy
        return 0.8

    def _trigger_policy_halt(self, violation_type: str, details: Dict):
        """Trigger policy halt."""
        self._log_audit("policy_violation", details, False)
        print(f"Policy violation: {violation_type}")

    def _log_audit(self, event_type: str, event_data: Dict, result: bool):
        """Log to audit."""
        query = """
        INSERT INTO omnipath_audit_log (event_type, application, event_data, result, hipaa_deidentified)
        VALUES (?, 'neurobotanica', ?, ?, 1)
        """
        if self.db:
            stmt = self.db.prepare(query).bind(event_type, str(event_data), result, 1)
            stmt.run()
        else:
            conn = self._get_connection()
            if 'psycopg2' in str(type(conn)):
                query = query.replace('?', '%s')
                cursor = conn.cursor()
                cursor.execute(query, (event_type, str(event_data), result))
                conn.commit()
                cursor.close()
                conn.close()
            else:
                conn.execute(query, (event_type, str(event_data), result))
                conn.commit()
                conn.close()

    def _log_formulation(self, formulation: Dict, tier: str):
        """Log to formulations."""
        import uuid
        fid = str(uuid.uuid4())
        query = """
        INSERT INTO neurobotanica_formulations 
        (formulation_id, formulation_name, ingredients, target_condition, synergy_score, tk_enhanced, requires_consent, created_at)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """
        ingredients = ",".join(formulation["cannabinoid_ratios"].keys())
        if self.db:
            stmt = self.db.prepare(query).bind(fid, f"{formulation['plant_id']}_formulation", ingredients, formulation["target_condition"], formulation["synergy_score"], formulation["tk_enhanced"], 1 if tier == "tk_enhanced" else 0, datetime.now().isoformat())
            stmt.run()
        else:
            conn = self._get_connection()
            if 'psycopg2' in str(type(conn)):
                query = query.replace('?', '%s')
                cursor = conn.cursor()
                cursor.execute(query, (fid, f"{formulation['plant_id']}_formulation", ingredients, formulation["target_condition"], formulation["synergy_score"], formulation["tk_enhanced"], 1 if tier == "tk_enhanced" else 0, datetime.now()))
                conn.commit()
                cursor.close()
                conn.close()
            else:
                conn.execute(query, (fid, f"{formulation['plant_id']}_formulation", ingredients, formulation["target_condition"], formulation["synergy_score"], formulation["tk_enhanced"], 1 if tier == "tk_enhanced" else 0, datetime.now()))
                conn.commit()
                conn.close()

    def close(self):
        pass

# Example usage
if __name__ == "__main__":
    analyzer = WholePlantAnalysisEngine()
    result = analyzer.analyze_whole_plant("cbd", "anxiety", "tk_enhanced")
    print(result)
    analyzer.close()