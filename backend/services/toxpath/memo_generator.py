"""
ToxPath Memo Generator - Assessment Report Generation

Generates formatted toxicology assessment memos from ToxPath results.
"""

from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from typing import Optional, Dict, Any
import json

from .assessor import ToxPathResponse, RiskTier, AlertSeverity


class MemoFormat(Enum):
    """Output format for memos."""
    MARKDOWN = "markdown"
    HTML = "html"
    JSON = "json"


@dataclass
class MemoConfig:
    """Configuration for memo generation."""
    include_testing_costs: bool = True
    include_timeline: bool = True
    include_consultation_section: bool = True
    include_assumptions: bool = True
    company_name: Optional[str] = None
    prepared_by: Optional[str] = None


class ToxPathMemoGenerator:
    """Generate formatted toxicology assessment memos."""
    
    MEMO_VERSION = "1.0.0"
    
    def __init__(self, config: Optional[MemoConfig] = None):
        self.config = config or MemoConfig()
    
    def generate(
        self, 
        assessment: ToxPathResponse,
        format: MemoFormat = MemoFormat.MARKDOWN
    ) -> str:
        """Generate memo in specified format."""
        if format == MemoFormat.MARKDOWN:
            return self._generate_markdown(assessment)
        elif format == MemoFormat.HTML:
            return self._generate_html(assessment)
        elif format == MemoFormat.JSON:
            return self._generate_json(assessment)
        else:
            raise ValueError(f"Unsupported format: {format}")
    
    def _generate_markdown(self, assessment: ToxPathResponse) -> str:
        """Generate Markdown memo."""
        lines = []
        
        # Header
        lines.append("# ToxPath Toxicity Risk Assessment Memo")
        lines.append("")
        lines.append(f"**Assessment ID:** `{assessment.toxpath_assessment_id}`")
        lines.append(f"**Generated:** {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')}")
        if self.config.company_name:
            lines.append(f"**Organization:** {self.config.company_name}")
        if self.config.prepared_by:
            lines.append(f"**Prepared by:** {self.config.prepared_by}")
        lines.append("")
        
        # Compound Information
        lines.append("---")
        lines.append("## Compound Information")
        lines.append("")
        lines.append(f"- **Compound Reference:** `{assessment.compound_ref}`")
        if assessment.compound_name:
            lines.append(f"- **Compound Name:** {assessment.compound_name}")
        lines.append(f"- **Administration Route:** {assessment.route.capitalize()}")
        lines.append("")
        
        # Properties snapshot
        props = assessment.properties_snapshot
        if props.get("is_valid"):
            lines.append("### Molecular Properties")
            lines.append("")
            lines.append("| Property | Value |")
            lines.append("|----------|-------|")
            lines.append(f"| Molecular Weight | {props.get('molecular_weight', 'N/A')} |")
            lines.append(f"| LogP | {props.get('logp', 'N/A')} |")
            lines.append(f"| TPSA | {props.get('tpsa', 'N/A')} |")
            lines.append(f"| H-Bond Donors | {props.get('hbd', 'N/A')} |")
            lines.append(f"| H-Bond Acceptors | {props.get('hba', 'N/A')} |")
            lines.append(f"| Rotatable Bonds | {props.get('rotatable_bonds', 'N/A')} |")
            lines.append("")
        
        # Risk Summary - prominent
        lines.append("---")
        lines.append("## ⚠️ Risk Assessment Summary")
        lines.append("")
        
        risk = assessment.risk_summary
        tier_emoji = self._get_tier_emoji(risk.overall_tier)
        lines.append(f"### {tier_emoji} Overall Risk Tier: **{risk.overall_tier.label}**")
        lines.append("")
        lines.append(f"**Testing Depth Recommendation:** {risk.overall_tier.testing_depth}")
        lines.append("")
        lines.append(f"**Rationale:** {risk.tier_rationale}")
        lines.append("")
        
        # Top Risks
        if risk.top_risks:
            lines.append("### Top Identified Risks")
            lines.append("")
            for risk_item in risk.top_risks:
                lines.append(f"- {risk_item}")
            lines.append("")
        
        # Key Unknowns
        if risk.key_unknowns:
            lines.append("### Key Unknowns")
            lines.append("")
            for unknown in risk.key_unknowns:
                lines.append(f"- {unknown}")
            lines.append("")
        
        # Route-specific concerns
        if risk.route_specific_concerns:
            lines.append("### Route-Specific Concerns")
            lines.append("")
            for concern in risk.route_specific_concerns:
                lines.append(f"- {concern}")
            lines.append("")
        
        # Structural Alerts
        if assessment.alerts:
            lines.append("---")
            lines.append("## Structural Alerts Detected")
            lines.append("")
            lines.append("| Alert | Severity | Mechanism | Affected Organs |")
            lines.append("|-------|----------|-----------|-----------------|")
            for alert in assessment.alerts:
                severity_badge = self._get_severity_badge(alert.severity)
                organs = ", ".join(alert.affected_organs) if alert.affected_organs else "N/A"
                mechanism = alert.mechanism or "N/A"
                lines.append(f"| {alert.name} | {severity_badge} | {mechanism} | {organs} |")
            lines.append("")
        else:
            lines.append("---")
            lines.append("## Structural Alerts")
            lines.append("")
            lines.append("✅ No significant structural alerts detected.")
            lines.append("")
        
        # Testing Plan
        lines.append("---")
        lines.append("## Recommended Testing Plan")
        lines.append("")
        
        if assessment.testing_plan:
            for step in assessment.testing_plan:
                req_badge = "🔴 Required" if step.required else "🟡 Recommended"
                glp_badge = " | GLP Required" if step.gmp_glp_required else ""
                
                lines.append(f"### {step.order}. {step.test_name}")
                lines.append("")
                lines.append(f"**Type:** {step.test_type.replace('_', ' ').title()} {req_badge}{glp_badge}")
                lines.append("")
                lines.append(f"**Rationale:** {step.rationale}")
                lines.append("")
                if self.config.include_testing_costs:
                    lines.append(f"**Estimated Cost:** {step.estimated_cost_range}")
                if self.config.include_timeline:
                    lines.append(f"**Timeline:** {step.timeline}")
                lines.append("")
        else:
            lines.append("No specific testing recommended at this time.")
            lines.append("")
        
        # Consultation Section
        if self.config.include_consultation_section:
            lines.append("---")
            lines.append("## Consultation Requirements")
            lines.append("")
            if assessment.consultation_required:
                lines.append("⚠️ **Expert Consultation Recommended**")
                lines.append("")
                for flag in assessment.consultation_flags:
                    lines.append(f"- {flag}")
            else:
                lines.append("✅ Standard development pathway may proceed without immediate consultation.")
            lines.append("")
        
        # Assumptions and Limitations
        if self.config.include_assumptions and risk.key_assumptions:
            lines.append("---")
            lines.append("## Assumptions and Limitations")
            lines.append("")
            for assumption in risk.key_assumptions:
                lines.append(f"- {assumption}")
            lines.append("")
        
        # Disclaimer
        lines.append("---")
        lines.append("## Disclaimer")
        lines.append("")
        lines.append("*This assessment is based on computational structural analysis and does not ")
        lines.append("replace formal toxicology studies. Results should be interpreted by qualified ")
        lines.append("toxicologists. Regulatory submissions require appropriate GLP studies. ")
        lines.append("This memo does not constitute medical or regulatory advice.*")
        lines.append("")
        lines.append("---")
        lines.append(f"*Generated by ToxPath v{self.MEMO_VERSION} | NeuroBotanica Platform*")
        
        return "\n".join(lines)
    
    def _generate_html(self, assessment: ToxPathResponse) -> str:
        """Generate HTML memo."""
        md_content = self._generate_markdown(assessment)
        
        # Basic HTML wrapper with styling
        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>ToxPath Assessment - {assessment.toxpath_assessment_id}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            max-width: 900px;
            margin: 0 auto;
            padding: 20px;
            line-height: 1.6;
            color: #333;
        }}
        h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        h3 {{ color: #7f8c8d; }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 15px 0;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 10px;
            text-align: left;
        }}
        th {{ background-color: #f5f6fa; }}
        code {{
            background-color: #f5f6fa;
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 0.9em;
        }}
        .risk-low {{ color: #27ae60; font-weight: bold; }}
        .risk-moderate {{ color: #f39c12; font-weight: bold; }}
        .risk-high {{ color: #e67e22; font-weight: bold; }}
        .risk-very-high {{ color: #e74c3c; font-weight: bold; }}
        hr {{ border: none; border-top: 1px solid #eee; margin: 30px 0; }}
        .disclaimer {{
            background-color: #f9f9f9;
            padding: 15px;
            border-left: 4px solid #95a5a6;
            font-style: italic;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <div class="memo-content">
        {self._markdown_to_html(md_content)}
    </div>
</body>
</html>"""
        return html
    
    def _markdown_to_html(self, md: str) -> str:
        """Basic markdown to HTML conversion."""
        import re
        
        html = md
        
        # Headers
        html = re.sub(r'^### (.+)$', r'<h3>\1</h3>', html, flags=re.MULTILINE)
        html = re.sub(r'^## (.+)$', r'<h2>\1</h2>', html, flags=re.MULTILINE)
        html = re.sub(r'^# (.+)$', r'<h1>\1</h1>', html, flags=re.MULTILINE)
        
        # Bold
        html = re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', html)
        
        # Italic
        html = re.sub(r'\*(.+?)\*', r'<em>\1</em>', html)
        
        # Code
        html = re.sub(r'`(.+?)`', r'<code>\1</code>', html)
        
        # Tables (basic)
        lines = html.split('\n')
        in_table = False
        new_lines = []
        for line in lines:
            if '|' in line and line.strip().startswith('|'):
                if not in_table:
                    new_lines.append('<table>')
                    in_table = True
                if '---' in line:
                    continue  # Skip separator
                cells = [c.strip() for c in line.split('|')[1:-1]]
                # First row after table start is header
                is_header_row = len(new_lines) > 0 and new_lines[-1] == '<table>'
                tag = 'th' if is_header_row else 'td'
                row = '<tr>' + ''.join(f'<{tag}>{c}</{tag}>' for c in cells) + '</tr>'
                new_lines.append(row)
            else:
                if in_table:
                    new_lines.append('</table>')
                    in_table = False
                new_lines.append(line)
        if in_table:
            new_lines.append('</table>')
        html = '\n'.join(new_lines)
        
        # Lists
        html = re.sub(r'^- (.+)$', r'<li>\1</li>', html, flags=re.MULTILINE)
        
        # Horizontal rules
        html = html.replace('---', '<hr>')
        
        # Paragraphs (simple)
        html = re.sub(r'\n\n', r'</p><p>', html)
        
        return f'<p>{html}</p>'
    
    def _generate_json(self, assessment: ToxPathResponse) -> str:
        """Generate JSON memo."""
        memo_data = {
            "memo_version": self.MEMO_VERSION,
            "generated_at": datetime.utcnow().isoformat(),
            "config": {
                "company_name": self.config.company_name,
                "prepared_by": self.config.prepared_by,
                "include_testing_costs": self.config.include_testing_costs,
                "include_timeline": self.config.include_timeline
            },
            "assessment": assessment.to_dict()
        }
        return json.dumps(memo_data, indent=2)
    
    def _get_tier_emoji(self, tier: RiskTier) -> str:
        """Get emoji for risk tier."""
        emojis = {
            RiskTier.LOW: "🟢",
            RiskTier.MODERATE: "🟡",
            RiskTier.HIGH: "🟠",
            RiskTier.VERY_HIGH: "🔴"
        }
        return emojis.get(tier, "⚪")
    
    def _get_severity_badge(self, severity: AlertSeverity) -> str:
        """Get badge for alert severity."""
        badges = {
            AlertSeverity.INFO: "ℹ️ Info",
            AlertSeverity.CAUTION: "⚠️ Caution",
            AlertSeverity.WARNING: "🔶 Warning",
            AlertSeverity.CRITICAL: "🔴 Critical"
        }
        return badges.get(severity, severity.value)
