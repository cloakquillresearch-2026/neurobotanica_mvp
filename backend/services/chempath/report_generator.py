"""
ChemPath Report Generator - Molecular Characterization Reports

Generates manufacturer-facing reports from ChemPath analysis results.
Supports Markdown output (PDF generation planned for future).
"""

from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from typing import Optional, List, Dict, Any

from .analyzer import ChemPathResponse, QCFlag, QCSeverity


class ReportFormat(Enum):
    """Supported report output formats."""
    MARKDOWN = "markdown"
    HTML = "html"
    JSON = "json"


@dataclass
class ReportSection:
    """A section of the report."""
    title: str
    content: str
    order: int


class ChemPathReportGenerator:
    """Generate molecular characterization reports.
    
    Report Sections:
    1. Identity & Structure
    2. Physicochemical Descriptors (2D)
    3. 3D Properties (if available)
    4. Certificate of Analysis Summary
    5. QC Flags & Data Quality
    6. Completeness Assessment
    """
    
    REPORT_TITLE = "Molecular Characterization Report"
    DISCLAIMER = "*For planning/compliance readiness only; lab confirmation required.*"
    
    def __init__(self):
        """Initialize report generator."""
        pass
    
    def generate_report(
        self,
        response: ChemPathResponse,
        format: ReportFormat = ReportFormat.MARKDOWN,
        include_disclaimer: bool = True
    ) -> str:
        """Generate full report from ChemPath analysis response."""
        if format == ReportFormat.MARKDOWN:
            return self._generate_markdown(response, include_disclaimer)
        elif format == ReportFormat.HTML:
            return self._generate_html(response, include_disclaimer)
        elif format == ReportFormat.JSON:
            return self._generate_json(response)
        else:
            raise ValueError(f"Unsupported format: {format}")
    
    def _generate_markdown(
        self,
        response: ChemPathResponse,
        include_disclaimer: bool
    ) -> str:
        """Generate Markdown report."""
        sections = []
        
        # Header
        sections.append(self._header_section(response))
        
        # Section 1: Identity & Structure
        sections.append(self._identity_section(response))
        
        # Section 2: 2D Descriptors
        sections.append(self._descriptors_2d_section(response))
        
        # Section 3: 3D Properties (if available)
        if response.computed_descriptors.descriptors_3d:
            sections.append(self._descriptors_3d_section(response))
        
        # Section 4: COA Summary (if available)
        if response.coa_qc_flags or self._has_coa_data(response):
            sections.append(self._coa_section(response))
        
        # Section 5: QC Flags
        sections.append(self._qc_section(response))
        
        # Section 6: Data Completeness
        sections.append(self._completeness_section(response))
        
        # Footer
        if include_disclaimer:
            sections.append(f"\n---\n{self.DISCLAIMER}")
        
        return "\n".join(sections)
    
    def _header_section(self, response: ChemPathResponse) -> str:
        """Generate report header."""
        timestamp = datetime.fromisoformat(response.analyzed_at).strftime("%Y-%m-%d %H:%M UTC")
        
        return f"""# {self.REPORT_TITLE}

**Compound:** {response.compound_name}  
**Generated:** {timestamp}  
**Job ID:** `{response.chempath_job_id}`  
**Status:** {response.status.upper()}
"""
    
    def _identity_section(self, response: ChemPathResponse) -> str:
        """Generate identity and structure section."""
        structure = response.normalized_structure
        
        section = """## 1. Identity & Structure

"""
        # Structure details
        section += f"**Canonical SMILES:** `{structure.canonical_smiles or 'N/A'}`\n\n"
        
        if structure.inchi_key:
            section += f"**InChIKey:** `{structure.inchi_key}`\n\n"
        
        if structure.inchi:
            # Truncate long InChI
            inchi_display = structure.inchi
            if len(inchi_display) > 80:
                inchi_display = inchi_display[:77] + "..."
            section += f"**InChI:** `{inchi_display}`\n\n"
        
        if structure.molecular_formula:
            section += f"**Molecular Formula:** {structure.molecular_formula}\n\n"
        
        # Normalization notes
        if structure.normalization_notes:
            section += "**Normalization Notes:**\n"
            for note in structure.normalization_notes:
                section += f"- {note}\n"
            section += "\n"
        
        # Validity status
        if not structure.is_valid:
            section += "⚠️ **Structure validation failed** - see QC flags below\n\n"
        
        return section
    
    def _descriptors_2d_section(self, response: ChemPathResponse) -> str:
        """Generate 2D descriptors section."""
        desc = response.computed_descriptors.descriptors_2d
        
        if not desc:
            return """## 2. Physicochemical Descriptors

*No 2D descriptors computed (structure validation failed)*

"""
        
        section = """## 2. Physicochemical Descriptors

### 2D Properties

| Property | Value |
|----------|-------|
"""
        # Key properties
        property_map = {
            "molecular_weight": ("Molecular Weight", "Da"),
            "exact_mass": ("Exact Mass", "Da"),
            "logp": ("LogP", ""),
            "tpsa": ("TPSA", "Å²"),
            "hbd": ("H-Bond Donors", ""),
            "hba": ("H-Bond Acceptors", ""),
            "rotatable_bonds": ("Rotatable Bonds", ""),
            "num_rings": ("Ring Count", ""),
            "num_aromatic_rings": ("Aromatic Rings", ""),
            "fraction_csp3": ("Fraction Csp³", ""),
            "num_stereocenters": ("Stereocenters", ""),
            "lipinski_violations": ("Lipinski Violations", ""),
            "qed": ("QED Score", ""),
            "labute_asa": ("Labute ASA", "Å²"),
            "heavy_atom_count": ("Heavy Atoms", "")
        }
        
        for key, (name, unit) in property_map.items():
            if key in desc:
                value = desc[key]
                if unit:
                    section += f"| {name} | {value} {unit} |\n"
                else:
                    section += f"| {name} | {value} |\n"
        
        # Drug-likeness assessment
        violations = desc.get("lipinski_violations", 0)
        qed = desc.get("qed", 0)
        
        section += "\n### Drug-Likeness Assessment\n\n"
        
        if violations == 0:
            section += "✅ **Passes Lipinski Rule of 5** (0 violations)\n"
        elif violations <= 1:
            section += f"⚠️ **Minor Lipinski violations** ({violations} violation)\n"
        else:
            section += f"❌ **Multiple Lipinski violations** ({violations} violations)\n"
        
        if qed >= 0.5:
            section += f"✅ **Good QED Score** ({qed})\n"
        elif qed > 0:
            section += f"⚠️ **Moderate QED Score** ({qed})\n"
        
        section += "\n"
        return section
    
    def _descriptors_3d_section(self, response: ChemPathResponse) -> str:
        """Generate 3D descriptors section."""
        desc_3d = response.computed_descriptors.descriptors_3d
        meta = response.computed_descriptors.conformer_meta
        
        section = """### 3D Properties

"""
        if meta and meta.get("error"):
            section += f"*3D generation failed: {meta['error']}*\n\n"
            return section
        
        if meta:
            section += f"**Method:** {meta.get('method', 'Unknown')}\n"
            section += f"**Optimization:** {meta.get('optimization', 'None')}\n\n"
        
        if desc_3d:
            section += "| Property | Value |\n"
            section += "|----------|-------|\n"
            
            property_names = {
                "radius_of_gyration": "Radius of Gyration",
                "inertial_shape_factor": "Inertial Shape Factor",
                "eccentricity": "Eccentricity",
                "asphericity": "Asphericity",
                "spherocity_index": "Spherocity Index",
                "npr1": "NPR1",
                "npr2": "NPR2"
            }
            
            for key, name in property_names.items():
                if key in desc_3d:
                    section += f"| {name} | {desc_3d[key]} |\n"
        
        section += "\n"
        return section
    
    def _coa_section(self, response: ChemPathResponse) -> str:
        """Generate COA summary section."""
        section = """## 3. Certificate of Analysis Summary

"""
        # Note: In full implementation, this would include actual COA data
        # MVP version focuses on QC flag summary
        
        coa_flags = [f for f in response.coa_qc_flags]
        
        if not coa_flags:
            section += "*No COA data provided for analysis*\n\n"
            return section
        
        # Categorize flags
        errors = [f for f in coa_flags if f.severity == QCSeverity.ERROR]
        warnings = [f for f in coa_flags if f.severity == QCSeverity.WARNING]
        infos = [f for f in coa_flags if f.severity == QCSeverity.INFO]
        
        if errors:
            section += f"### ❌ Errors ({len(errors)})\n"
            for flag in errors:
                section += f"- **{flag.code}**: {flag.message}\n"
            section += "\n"
        
        if warnings:
            section += f"### ⚠️ Warnings ({len(warnings)})\n"
            for flag in warnings:
                section += f"- **{flag.code}**: {flag.message}\n"
            section += "\n"
        
        if infos:
            section += f"### ℹ️ Information ({len(infos)})\n"
            for flag in infos:
                section += f"- **{flag.code}**: {flag.message}\n"
            section += "\n"
        
        if not errors and not warnings:
            section += "✅ **COA passes quality checks**\n\n"
        
        return section
    
    def _qc_section(self, response: ChemPathResponse) -> str:
        """Generate QC flags section."""
        section = """## 4. Structure QC Flags

"""
        structure_flags = response.structure_qc_flags
        
        if not structure_flags:
            section += "✅ **No structure quality issues detected**\n\n"
            return section
        
        # Categorize flags
        errors = [f for f in structure_flags if f.severity == QCSeverity.ERROR]
        warnings = [f for f in structure_flags if f.severity == QCSeverity.WARNING]
        infos = [f for f in structure_flags if f.severity == QCSeverity.INFO]
        
        section += "| Severity | Code | Message |\n"
        section += "|----------|------|---------|"
        
        for flag in errors:
            section += f"\n| ❌ Error | `{flag.code}` | {flag.message} |"
        
        for flag in warnings:
            section += f"\n| ⚠️ Warning | `{flag.code}` | {flag.message} |"
        
        for flag in infos:
            section += f"\n| ℹ️ Info | `{flag.code}` | {flag.message} |"
        
        section += "\n\n"
        return section
    
    def _completeness_section(self, response: ChemPathResponse) -> str:
        """Generate data completeness section."""
        score = response.data_completeness_score
        
        section = f"""## 5. Data Completeness

**Completeness Score:** {score}%

"""
        # Visual indicator
        filled = int(score / 10)
        empty = 10 - filled
        bar = "█" * filled + "░" * empty
        section += f"`{bar}` {score}%\n\n"
        
        # Assessment
        if score >= 80:
            section += "✅ **Excellent** - Data is comprehensive for regulatory planning\n"
        elif score >= 60:
            section += "⚠️ **Good** - Most key data present, some gaps exist\n"
        elif score >= 40:
            section += "⚠️ **Moderate** - Significant data gaps, recommend additional analysis\n"
        else:
            section += "❌ **Insufficient** - Major data gaps, additional characterization required\n"
        
        # Recommendations
        section += "\n### Recommendations\n\n"
        
        if score < 100:
            if not response.normalized_structure.inchi_key:
                section += "- Generate InChIKey for structure registry lookup\n"
            if not response.computed_descriptors.descriptors_3d:
                section += "- Compute 3D conformers for shape-based analysis\n"
            if not response.coa_qc_flags and score < 60:
                section += "- Provide Certificate of Analysis for compliance assessment\n"
        
        if score >= 80:
            section += "- Ready for ToxPath assessment\n"
            section += "- Consider RegPath pathway analysis\n"
        
        section += "\n"
        return section
    
    def _has_coa_data(self, response: ChemPathResponse) -> bool:
        """Check if response has COA data."""
        return len(response.coa_qc_flags) > 0
    
    def _generate_html(
        self,
        response: ChemPathResponse,
        include_disclaimer: bool
    ) -> str:
        """Generate HTML report (converts markdown to HTML)."""
        # Simple HTML wrapper around markdown
        markdown_content = self._generate_markdown(response, include_disclaimer)
        
        html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>ChemPath Report - {response.compound_name}</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; 
               max-width: 900px; margin: 0 auto; padding: 20px; }}
        h1 {{ color: #1a1a1a; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #2c3e50; margin-top: 30px; }}
        h3 {{ color: #34495e; }}
        table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #3498db; color: white; }}
        code {{ background-color: #f4f4f4; padding: 2px 6px; border-radius: 3px; }}
        pre {{ background-color: #f4f4f4; padding: 15px; border-radius: 5px; overflow-x: auto; }}
        .disclaimer {{ font-style: italic; color: #666; margin-top: 30px; padding-top: 20px; 
                      border-top: 1px solid #ddd; }}
    </style>
</head>
<body>
<pre>{markdown_content}</pre>
</body>
</html>"""
        
        return html
    
    def _generate_json(self, response: ChemPathResponse) -> str:
        """Generate JSON report."""
        import json
        
        report_data = {
            "report_type": "molecular_characterization",
            "generated_at": datetime.utcnow().isoformat(),
            "data": response.to_dict()
        }
        
        return json.dumps(report_data, indent=2)
    
    def generate_summary(self, response: ChemPathResponse) -> Dict[str, Any]:
        """Generate a brief summary of the analysis."""
        structure = response.normalized_structure
        desc = response.computed_descriptors.descriptors_2d
        
        # Count issues
        total_flags = len(response.structure_qc_flags) + len(response.coa_qc_flags)
        error_flags = len([f for f in response.structure_qc_flags + response.coa_qc_flags 
                          if f.severity == QCSeverity.ERROR])
        
        return {
            "compound_name": response.compound_name,
            "job_id": response.chempath_job_id,
            "status": response.status,
            "structure_valid": structure.is_valid,
            "molecular_weight": desc.get("molecular_weight"),
            "logp": desc.get("logp"),
            "lipinski_violations": desc.get("lipinski_violations", 0),
            "qed_score": desc.get("qed"),
            "has_3d_data": response.computed_descriptors.descriptors_3d is not None,
            "completeness_score": response.data_completeness_score,
            "total_qc_flags": total_flags,
            "error_count": error_flags,
            "ready_for_toxpath": response.data_completeness_score >= 60 and error_flags == 0
        }
