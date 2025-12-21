"""
RegPath Memo Generator - Regulatory Strategy Memo Generation

Generates formatted regulatory strategy memos from RegPath results.
"""

from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from typing import Optional, Dict, Any
import json

from .strategist import (
    RegPathResponse, 
    RegulatoryPathway, 
    ChecklistCategory,
    MilestonePhase
)


class MemoFormat(Enum):
    """Output format for memos."""
    MARKDOWN = "markdown"
    HTML = "html"
    JSON = "json"


@dataclass
class MemoConfig:
    """Configuration for memo generation."""
    include_timeline_table: bool = True
    include_cost_estimates: bool = True
    include_checklist: bool = True
    include_assumptions: bool = True
    company_name: Optional[str] = None
    prepared_by: Optional[str] = None


class RegPathMemoGenerator:
    """Generate formatted regulatory strategy memos."""
    
    MEMO_VERSION = "1.0.0"
    
    def __init__(self, config: Optional[MemoConfig] = None):
        self.config = config or MemoConfig()
    
    def generate(
        self,
        strategy: RegPathResponse,
        format: MemoFormat = MemoFormat.MARKDOWN
    ) -> str:
        """Generate memo in specified format."""
        if format == MemoFormat.MARKDOWN:
            return self._generate_markdown(strategy)
        elif format == MemoFormat.HTML:
            return self._generate_html(strategy)
        elif format == MemoFormat.JSON:
            return self._generate_json(strategy)
        else:
            raise ValueError(f"Unsupported format: {format}")
    
    def _generate_markdown(self, strategy: RegPathResponse) -> str:
        """Generate Markdown memo."""
        lines = []
        profile = strategy.product_profile
        
        # Header
        lines.append("# Regulatory Strategy Memo")
        lines.append("")
        lines.append(f"**Product:** {profile.product_name}")
        lines.append(f"**Type:** {profile.product_type}")
        lines.append(f"**Route:** {profile.route.capitalize()}")
        lines.append(f"**Target Market:** {', '.join(profile.target_markets)}")
        lines.append(f"**Generated:** {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')}")
        lines.append(f"**Strategy ID:** `{strategy.regpath_strategy_id}`")
        if self.config.company_name:
            lines.append(f"**Organization:** {self.config.company_name}")
        if self.config.prepared_by:
            lines.append(f"**Prepared by:** {self.config.prepared_by}")
        lines.append("")
        
        # Product Snapshot
        lines.append("---")
        lines.append("## 1. Product Snapshot")
        lines.append("")
        lines.append(f"**Product Name:** {profile.product_name}")
        lines.append(f"**Product Type:** {profile.product_type}")
        lines.append(f"**Intended Use:** {profile.intended_use}")
        lines.append(f"**Route of Administration:** {profile.route}")
        if profile.therapeutic_area:
            lines.append(f"**Therapeutic Area:** {profile.therapeutic_area}")
        if profile.active_ingredients:
            lines.append(f"**Active Ingredients:** {', '.join(profile.active_ingredients)}")
        lines.append("")
        lines.append(f"**Novelty Assessment:** {strategy.novelty_assessment.capitalize()}")
        lines.append(f"**Evidence Strength:** {strategy.evidence_strength.capitalize()}")
        lines.append("")
        
        # Key Assumptions
        if self.config.include_assumptions and strategy.key_assumptions:
            lines.append("### Key Assumptions")
            lines.append("")
            for assumption in strategy.key_assumptions:
                lines.append(f"- {assumption}")
            lines.append("")
        
        # Pathway Recommendation
        lines.append("---")
        lines.append("## 2. Pathway Recommendation")
        lines.append("")
        
        primary = strategy.primary_pathway
        lines.append(f"### Primary: {primary.pathway.value} ({primary.pathway.full_name})")
        lines.append("")
        lines.append(f"**Confidence:** {primary.confidence.capitalize()}")
        lines.append("")
        lines.append("**Rationale:**")
        for point in primary.rationale:
            lines.append(f"- {point}")
        lines.append("")
        
        lines.append("**Key Requirements:**")
        for req in primary.key_requirements:
            lines.append(f"- {req}")
        lines.append("")
        
        if self.config.include_cost_estimates:
            lines.append(f"**Estimated Cost:** {primary.estimated_cost_range}")
            timeline = primary.estimated_timeline_months
            lines.append(f"**Estimated Timeline:** {timeline[0]}-{timeline[1]} months")
            lines.append("")
        
        # Fallback pathway
        if strategy.fallback_pathway:
            fallback = strategy.fallback_pathway
            lines.append(f"### Fallback: {fallback.pathway.value} ({fallback.pathway.full_name})")
            lines.append("")
            lines.append(f"**Confidence:** {fallback.confidence.capitalize()}")
            lines.append("")
            lines.append("**Rationale:**")
            for point in fallback.rationale:
                lines.append(f"- {point}")
            lines.append("")
            
            if self.config.include_cost_estimates:
                lines.append(f"**Estimated Cost:** {fallback.estimated_cost_range}")
                timeline = fallback.estimated_timeline_months
                lines.append(f"**Estimated Timeline:** {timeline[0]}-{timeline[1]} months")
                lines.append("")
        
        # Critical Gating Questions
        if strategy.gating_questions:
            lines.append("### Critical Gating Questions")
            lines.append("")
            for i, question in enumerate(strategy.gating_questions, 1):
                lines.append(f"{i}. {question}")
            lines.append("")
        
        # Readiness Checklist
        if self.config.include_checklist and strategy.readiness_checklist:
            lines.append("---")
            lines.append("## 3. Readiness Checklist")
            lines.append("")
            
            # Group by category
            categories = {}
            for item in strategy.readiness_checklist:
                cat = item.category.value
                if cat not in categories:
                    categories[cat] = []
                categories[cat].append(item)
            
            for cat_name, items in categories.items():
                lines.append(f"### {cat_name}")
                lines.append("")
                lines.append("| Item | Description | Priority | Status |")
                lines.append("|------|-------------|----------|--------|")
                for item in items:
                    status_emoji = self._get_status_emoji(item.status)
                    lines.append(f"| {item.item} | {item.description} | {item.priority.capitalize()} | {status_emoji} |")
                lines.append("")
        
        # Timeline
        if self.config.include_timeline_table and strategy.timeline:
            lines.append("---")
            lines.append("## 4. Estimated Timeline")
            lines.append("")
            lines.append("| Phase | Milestone | Duration | Start | End | Deliverables |")
            lines.append("|-------|-----------|----------|-------|-----|--------------|")
            for milestone in strategy.timeline:
                end_month = milestone.start_month + milestone.duration_months
                deliverables = ", ".join(milestone.deliverables) if milestone.deliverables else "-"
                lines.append(f"| {milestone.phase.value} | {milestone.milestone} | {milestone.duration_months} mo | M{milestone.start_month} | M{end_month} | {deliverables} |")
            lines.append("")
            
            # Calculate total timeline
            if strategy.timeline:
                total_months = max(m.start_month + m.duration_months for m in strategy.timeline)
                lines.append(f"**Total Estimated Duration:** {total_months} months ({total_months/12:.1f} years)")
                lines.append("")
        
        # Next Actions
        lines.append("---")
        lines.append("## 5. Next Actions")
        lines.append("")
        for i, action in enumerate(strategy.next_actions, 1):
            lines.append(f"{i}. {action}")
        lines.append("")
        
        # Disclaimer
        lines.append("---")
        lines.append("## Disclaimer")
        lines.append("")
        lines.append("*Timeline estimates are based on typical pathways and may vary significantly ")
        lines.append("based on product complexity, FDA feedback, and development outcomes. ")
        lines.append("Consult regulatory counsel for binding guidance. Cost estimates are indicative ")
        lines.append("and do not constitute a formal budget or quote.*")
        lines.append("")
        lines.append("---")
        lines.append(f"*Generated by RegPath v{self.MEMO_VERSION} | NeuroBotanica Platform*")
        
        return "\n".join(lines)
    
    def _generate_html(self, strategy: RegPathResponse) -> str:
        """Generate HTML memo."""
        md_content = self._generate_markdown(strategy)
        
        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Regulatory Strategy - {strategy.product_profile.product_name}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            max-width: 900px;
            margin: 0 auto;
            padding: 20px;
            line-height: 1.6;
            color: #333;
        }}
        h1 {{ color: #1a5f7a; border-bottom: 3px solid #1a5f7a; padding-bottom: 10px; }}
        h2 {{ color: #2c3e50; margin-top: 30px; border-left: 4px solid #3498db; padding-left: 15px; }}
        h3 {{ color: #34495e; }}
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
        th {{ background-color: #f5f6fa; font-weight: 600; }}
        tr:nth-child(even) {{ background-color: #fafafa; }}
        code {{
            background-color: #f5f6fa;
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 0.9em;
        }}
        .pathway-primary {{
            background-color: #e8f5e9;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
        }}
        .pathway-fallback {{
            background-color: #fff3e0;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
        }}
        hr {{ border: none; border-top: 1px solid #eee; margin: 30px 0; }}
        .disclaimer {{
            background-color: #f9f9f9;
            padding: 15px;
            border-left: 4px solid #95a5a6;
            font-style: italic;
            font-size: 0.9em;
        }}
        .next-action {{
            background-color: #e3f2fd;
            padding: 8px 15px;
            margin: 5px 0;
            border-radius: 4px;
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
        
        # Tables
        lines = html.split('\n')
        in_table = False
        new_lines = []
        for line in lines:
            if '|' in line and line.strip().startswith('|'):
                if not in_table:
                    new_lines.append('<table>')
                    in_table = True
                if '---' in line:
                    continue
                cells = [c.strip() for c in line.split('|')[1:-1]]
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
        html = re.sub(r'^(\d+)\. (.+)$', r'<li>\2</li>', html, flags=re.MULTILINE)
        
        # Horizontal rules
        html = html.replace('---', '<hr>')
        
        return html
    
    def _generate_json(self, strategy: RegPathResponse) -> str:
        """Generate JSON memo."""
        memo_data = {
            "memo_version": self.MEMO_VERSION,
            "generated_at": datetime.utcnow().isoformat(),
            "config": {
                "company_name": self.config.company_name,
                "prepared_by": self.config.prepared_by,
                "include_cost_estimates": self.config.include_cost_estimates,
                "include_timeline_table": self.config.include_timeline_table
            },
            "strategy": strategy.to_dict()
        }
        return json.dumps(memo_data, indent=2)
    
    def _get_status_emoji(self, status: str) -> str:
        """Get emoji for status."""
        emojis = {
            "not_started": "⬜ Not Started",
            "in_progress": "🔄 In Progress",
            "complete": "✅ Complete",
            "na": "➖ N/A"
        }
        return emojis.get(status, status)
