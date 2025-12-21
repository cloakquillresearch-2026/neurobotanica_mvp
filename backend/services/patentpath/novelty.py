"""
Novelty Scoring - PatentPath Lite
Assesses novelty of cannabinoid innovations against prior art.

Features:
- Similarity-based scoring
- Claim element analysis
- CPC code matching
- Multi-factor novelty assessment
"""
import asyncio
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
import logging
import re
from collections import Counter

from .prior_art import PriorArtSearcher, PriorArtResult, SearchQuery

logger = logging.getLogger(__name__)


class NoveltyLevel(Enum):
    """Novelty assessment levels."""
    HIGH = "high"           # Novel - significant patentability potential
    MODERATE = "moderate"   # Some prior art, but distinguishable
    LOW = "low"            # Substantial overlap with prior art
    BLOCKED = "blocked"    # Direct conflict with existing patents


@dataclass
class SimilarityMatch:
    """Single similarity match with prior art."""
    patent_number: str
    patent_title: str
    similarity_score: float
    matching_keywords: List[str]
    matching_cpc_codes: List[str]
    concern_level: str  # "high", "moderate", "low"
    notes: str
    url: str


@dataclass
class NoveltyReport:
    """Comprehensive novelty assessment report."""
    
    # Innovation details
    innovation_title: str
    innovation_description: str
    innovation_keywords: List[str]
    innovation_cpc_codes: List[str]
    
    # Scores (0-1)
    overall_novelty_score: float
    keyword_novelty_score: float
    structural_novelty_score: float
    claim_novelty_score: float
    
    # Assessment
    novelty_level: NoveltyLevel
    patentability_potential: str  # "Strong", "Moderate", "Weak"
    
    # Prior art analysis
    prior_art_count: int
    high_relevance_count: int
    top_matches: List[SimilarityMatch]
    
    # Recommendations
    recommendations: List[str]
    differentiation_opportunities: List[str]
    
    # Metadata
    generated_at: str = field(default_factory=lambda: datetime.now().isoformat())
    search_sources: List[str] = field(default_factory=list)
    
    def to_dict(self) -> Dict:
        return {
            "innovation": {
                "title": self.innovation_title,
                "description": self.innovation_description[:500] if self.innovation_description else "",
                "keywords": self.innovation_keywords,
                "cpc_codes": self.innovation_cpc_codes
            },
            "scores": {
                "overall_novelty": round(self.overall_novelty_score, 3),
                "keyword_novelty": round(self.keyword_novelty_score, 3),
                "structural_novelty": round(self.structural_novelty_score, 3),
                "claim_novelty": round(self.claim_novelty_score, 3)
            },
            "assessment": {
                "novelty_level": self.novelty_level.value,
                "patentability_potential": self.patentability_potential
            },
            "prior_art": {
                "total_found": self.prior_art_count,
                "high_relevance_count": self.high_relevance_count,
                "top_matches": [
                    {
                        "patent_number": m.patent_number,
                        "title": m.patent_title,
                        "similarity": round(m.similarity_score, 3),
                        "concern_level": m.concern_level,
                        "url": m.url
                    }
                    for m in self.top_matches[:5]
                ]
            },
            "recommendations": self.recommendations,
            "differentiation_opportunities": self.differentiation_opportunities,
            "metadata": {
                "generated_at": self.generated_at,
                "sources": self.search_sources
            }
        }
    
    def get_summary(self) -> str:
        """Generate a human-readable summary."""
        return f"""
Novelty Assessment: {self.innovation_title}
====================================
Overall Novelty Score: {self.overall_novelty_score:.1%}
Novelty Level: {self.novelty_level.value.upper()}
Patentability Potential: {self.patentability_potential}

Prior Art Analysis:
- Total patents found: {self.prior_art_count}
- High relevance matches: {self.high_relevance_count}

Top Recommendations:
{chr(10).join(f"• {r}" for r in self.recommendations[:3])}
""".strip()


class NoveltyScorer:
    """Novelty assessment engine for cannabinoid innovations."""
    
    # Keywords indicating different innovation types
    INNOVATION_KEYWORDS = {
        "composition": ["formulation", "composition", "mixture", "combination", "blend"],
        "method": ["method", "process", "treatment", "administering", "delivering"],
        "compound": ["compound", "derivative", "analog", "dimer", "conjugate"],
        "device": ["device", "apparatus", "system", "inhaler", "vaporizer"],
        "use": ["use", "application", "treatment", "therapy"]
    }
    
    # CPC codes by technology area
    TECHNOLOGY_CPC_MAP = {
        "cannabinoid_compound": ["A61K31/352", "C07D311/80"],
        "formulation": ["A61K9/00", "A61K47/00"],
        "extraction": ["B01D11/02", "A61K36/185"],
        "delivery": ["A61K9/00", "A61M15/00"],
        "therapeutic_use": ["A61P25/00", "A61P29/00", "A61P35/00"]
    }
    
    def __init__(self, searcher: Optional[PriorArtSearcher] = None):
        """Initialize novelty scorer.
        
        Args:
            searcher: PriorArtSearcher instance (creates new one if None)
        """
        self.searcher = searcher or PriorArtSearcher()
    
    async def assess_novelty(
        self,
        title: str,
        description: str,
        keywords: Optional[List[str]] = None,
        claims: Optional[List[str]] = None,
        smiles: Optional[str] = None,
        cpc_codes: Optional[List[str]] = None
    ) -> NoveltyReport:
        """Assess novelty of an innovation.
        
        Args:
            title: Innovation title
            description: Detailed description
            keywords: Key terms (extracted if not provided)
            claims: Draft claims (optional)
            smiles: Chemical structure (optional)
            cpc_codes: Relevant CPC codes (inferred if not provided)
            
        Returns:
            NoveltyReport with comprehensive assessment
        """
        # Extract/validate keywords
        if not keywords:
            keywords = self._extract_keywords(title, description)
        
        # Infer CPC codes
        if not cpc_codes:
            cpc_codes = self._infer_cpc_codes(title, description, keywords)
        
        # Search prior art
        prior_art = await self._search_prior_art(keywords, cpc_codes)
        
        # Calculate similarity scores
        matches = self._calculate_similarities(
            title, description, keywords, claims, prior_art
        )
        
        # Calculate novelty scores
        keyword_score = self._calculate_keyword_novelty(keywords, prior_art)
        structural_score = self._calculate_structural_novelty(cpc_codes, prior_art)
        claim_score = self._calculate_claim_novelty(claims, prior_art) if claims else 0.7
        
        # Overall score (weighted average)
        overall_score = (
            keyword_score * 0.35 +
            structural_score * 0.35 +
            claim_score * 0.30
        )
        
        # Determine novelty level
        novelty_level = self._determine_novelty_level(overall_score, matches)
        
        # Generate recommendations
        recommendations = self._generate_recommendations(
            novelty_level, matches, keywords, cpc_codes
        )
        
        differentiation = self._identify_differentiation(
            description, prior_art, matches
        )
        
        # Patentability assessment
        patentability = self._assess_patentability(
            overall_score, novelty_level, len(prior_art)
        )
        
        return NoveltyReport(
            innovation_title=title,
            innovation_description=description,
            innovation_keywords=keywords,
            innovation_cpc_codes=cpc_codes,
            overall_novelty_score=overall_score,
            keyword_novelty_score=keyword_score,
            structural_novelty_score=structural_score,
            claim_novelty_score=claim_score,
            novelty_level=novelty_level,
            patentability_potential=patentability,
            prior_art_count=len(prior_art),
            high_relevance_count=sum(1 for m in matches if m.similarity_score > 0.6),
            top_matches=sorted(matches, key=lambda x: x.similarity_score, reverse=True)[:10],
            recommendations=recommendations,
            differentiation_opportunities=differentiation,
            search_sources=["USPTO PatentsView"]
        )
    
    def _extract_keywords(self, title: str, description: str) -> List[str]:
        """Extract key terms from title and description."""
        text = f"{title} {description}".lower()
        
        # Cannabis-related keywords to look for
        cannabis_terms = [
            "cannabinoid", "thc", "cbd", "cbg", "cbn", "terpene", "cannabis",
            "tetrahydrocannabinol", "cannabidiol", "cannabigerol", "cannabinol",
            "dimer", "linker", "receptor", "cb1", "cb2", "endocannabinoid",
            "formulation", "delivery", "extraction", "synthesis"
        ]
        
        keywords = []
        for term in cannabis_terms:
            if term in text:
                keywords.append(term)
        
        # Extract noun phrases (simplified)
        words = re.findall(r'\b[a-z]{4,}\b', text)
        word_freq = Counter(words)
        
        # Add frequent words (likely important)
        for word, freq in word_freq.most_common(10):
            if word not in keywords and freq >= 2:
                keywords.append(word)
        
        return keywords[:15]  # Limit
    
    def _infer_cpc_codes(
        self,
        title: str,
        description: str,
        keywords: List[str]
    ) -> List[str]:
        """Infer relevant CPC codes from content."""
        text = f"{title} {description}".lower()
        codes = set()
        
        # Always include cannabinoid base codes
        codes.add("A61K31/352")  # Cannabinoids
        
        # Check for specific technology areas
        if any(w in text for w in ["formulation", "composition", "mixture"]):
            codes.update(self.TECHNOLOGY_CPC_MAP["formulation"])
        
        if any(w in text for w in ["extract", "isolation", "purif"]):
            codes.update(self.TECHNOLOGY_CPC_MAP["extraction"])
        
        if any(w in text for w in ["deliver", "release", "transdermal"]):
            codes.update(self.TECHNOLOGY_CPC_MAP["delivery"])
        
        if any(w in text for w in ["dimer", "derivative", "analog", "synthe"]):
            codes.update(self.TECHNOLOGY_CPC_MAP["cannabinoid_compound"])
        
        if any(w in text for w in ["treatment", "therapy", "patient"]):
            codes.update(self.TECHNOLOGY_CPC_MAP["therapeutic_use"])
        
        return list(codes)[:8]
    
    async def _search_prior_art(
        self,
        keywords: List[str],
        cpc_codes: List[str]
    ) -> List[PriorArtResult]:
        """Search for relevant prior art."""
        query = SearchQuery(
            keywords=keywords[:5],  # Top keywords
            cpc_codes=cpc_codes,
            max_results=50
        )
        
        return await self.searcher.search(query=query)
    
    def _calculate_similarities(
        self,
        title: str,
        description: str,
        keywords: List[str],
        claims: Optional[List[str]],
        prior_art: List[PriorArtResult]
    ) -> List[SimilarityMatch]:
        """Calculate similarity scores with prior art."""
        matches = []
        
        innovation_text = f"{title} {description}".lower()
        if claims:
            innovation_text += " " + " ".join(claims).lower()
        
        for pa in prior_art:
            pa_text = f"{pa.title} {pa.abstract}".lower()
            
            # Keyword matching
            matching_keywords = [
                kw for kw in keywords
                if kw.lower() in pa_text
            ]
            
            # CPC matching
            matching_cpcs = []
            for code in pa.classifications:
                for innovation_code in keywords:
                    if code.startswith(innovation_code[:7]):
                        matching_cpcs.append(code)
            
            # Calculate similarity score
            keyword_overlap = len(matching_keywords) / max(len(keywords), 1)
            relevance_factor = pa.relevance_score
            
            # Jaccard-like similarity on words
            innovation_words = set(re.findall(r'\b\w+\b', innovation_text))
            pa_words = set(re.findall(r'\b\w+\b', pa_text))
            
            if innovation_words and pa_words:
                jaccard = len(innovation_words & pa_words) / len(innovation_words | pa_words)
            else:
                jaccard = 0
            
            similarity = (keyword_overlap * 0.4 + relevance_factor * 0.3 + jaccard * 0.3)
            
            # Determine concern level
            if similarity > 0.7:
                concern = "high"
            elif similarity > 0.4:
                concern = "moderate"
            else:
                concern = "low"
            
            # Notes
            notes = ""
            if matching_keywords:
                notes = f"Matches keywords: {', '.join(matching_keywords[:3])}"
            
            matches.append(SimilarityMatch(
                patent_number=pa.patent_number,
                patent_title=pa.title,
                similarity_score=similarity,
                matching_keywords=matching_keywords,
                matching_cpc_codes=matching_cpcs,
                concern_level=concern,
                notes=notes,
                url=pa.url
            ))
        
        return matches
    
    def _calculate_keyword_novelty(
        self,
        keywords: List[str],
        prior_art: List[PriorArtResult]
    ) -> float:
        """Calculate keyword-based novelty score."""
        if not keywords:
            return 0.5
        
        # Count how many prior art docs contain each keyword
        keyword_counts = {kw: 0 for kw in keywords}
        
        for pa in prior_art:
            pa_text = f"{pa.title} {pa.abstract}".lower()
            for kw in keywords:
                if kw.lower() in pa_text:
                    keyword_counts[kw] += 1
        
        # Calculate novelty (lower counts = higher novelty)
        total_pa = max(len(prior_art), 1)
        keyword_novelties = []
        
        for kw, count in keyword_counts.items():
            # Inverse frequency = novelty
            novelty = 1 - (count / total_pa)
            keyword_novelties.append(novelty)
        
        return sum(keyword_novelties) / max(len(keyword_novelties), 1)
    
    def _calculate_structural_novelty(
        self,
        cpc_codes: List[str],
        prior_art: List[PriorArtResult]
    ) -> float:
        """Calculate structural/classification novelty."""
        if not cpc_codes:
            return 0.5
        
        # Check how common our CPC codes are in prior art
        cpc_counts = {code: 0 for code in cpc_codes}
        
        for pa in prior_art:
            for code in cpc_codes:
                for pa_code in pa.classifications:
                    if pa_code.startswith(code[:7]):
                        cpc_counts[code] += 1
                        break
        
        total_pa = max(len(prior_art), 1)
        novelties = []
        
        for code, count in cpc_counts.items():
            novelty = 1 - (count / total_pa)
            novelties.append(novelty)
        
        return sum(novelties) / max(len(novelties), 1)
    
    def _calculate_claim_novelty(
        self,
        claims: Optional[List[str]],
        prior_art: List[PriorArtResult]
    ) -> float:
        """Calculate claim-based novelty (simplified)."""
        if not claims:
            return 0.7  # Default moderate
        
        # Extract claim elements
        claim_text = " ".join(claims).lower()
        
        # Check overlap with prior art
        max_overlap = 0
        for pa in prior_art:
            pa_text = f"{pa.title} {pa.abstract}".lower()
            
            # Simple word overlap
            claim_words = set(re.findall(r'\b\w{4,}\b', claim_text))
            pa_words = set(re.findall(r'\b\w{4,}\b', pa_text))
            
            if claim_words:
                overlap = len(claim_words & pa_words) / len(claim_words)
                max_overlap = max(max_overlap, overlap)
        
        return 1 - max_overlap
    
    def _determine_novelty_level(
        self,
        overall_score: float,
        matches: List[SimilarityMatch]
    ) -> NoveltyLevel:
        """Determine novelty level from score and matches."""
        high_concern_count = sum(1 for m in matches if m.concern_level == "high")
        
        if high_concern_count >= 3 or overall_score < 0.3:
            return NoveltyLevel.BLOCKED
        elif high_concern_count >= 1 or overall_score < 0.5:
            return NoveltyLevel.LOW
        elif overall_score < 0.7:
            return NoveltyLevel.MODERATE
        else:
            return NoveltyLevel.HIGH
    
    def _assess_patentability(
        self,
        score: float,
        level: NoveltyLevel,
        prior_art_count: int
    ) -> str:
        """Assess overall patentability potential."""
        if level == NoveltyLevel.HIGH:
            return "Strong"
        elif level == NoveltyLevel.MODERATE:
            return "Moderate"
        elif level == NoveltyLevel.LOW:
            return "Weak"
        else:
            return "Unlikely"
    
    def _generate_recommendations(
        self,
        level: NoveltyLevel,
        matches: List[SimilarityMatch],
        keywords: List[str],
        cpc_codes: List[str]
    ) -> List[str]:
        """Generate actionable recommendations."""
        recommendations = []
        
        if level == NoveltyLevel.BLOCKED:
            recommendations.append(
                "Consider pivoting innovation approach - significant prior art overlap detected"
            )
            recommendations.append(
                "Review high-concern patents for potential licensing opportunities"
            )
            
        elif level == NoveltyLevel.LOW:
            recommendations.append(
                "Focus claims on specific differentiating features"
            )
            recommendations.append(
                "Consider narrowing scope to avoid prior art conflicts"
            )
            
        elif level == NoveltyLevel.MODERATE:
            recommendations.append(
                "Strengthen claims by emphasizing unique technical effects"
            )
            recommendations.append(
                "Consider continuation-in-part for additional differentiation"
            )
            
        else:  # HIGH
            recommendations.append(
                "Strong novelty position - proceed with patent application"
            )
            recommendations.append(
                "Consider broader claim scope to maximize protection"
            )
        
        # Add specific recommendations based on matches
        high_matches = [m for m in matches if m.concern_level == "high"]
        if high_matches:
            recommendations.append(
                f"Review {len(high_matches)} high-relevance patents: " +
                ", ".join(m.patent_number for m in high_matches[:3])
            )
        
        return recommendations
    
    def _identify_differentiation(
        self,
        description: str,
        prior_art: List[PriorArtResult],
        matches: List[SimilarityMatch]
    ) -> List[str]:
        """Identify potential differentiation opportunities."""
        opportunities = []
        
        desc_lower = description.lower()
        
        # Check for differentiating features
        differentiators = {
            "dimer": "Novel dimer linkage strategy could provide differentiation",
            "linker": "Specific linker chemistry may distinguish from prior art",
            "ratio": "Specific component ratios could support patentability",
            "synerg": "Synergistic effects may provide non-obvious improvement",
            "nano": "Nanoformulation approach could differentiate",
            "targeted": "Targeted delivery mechanism offers differentiation",
            "sustained": "Sustained release profile may be patentable",
            "enhanced": "Enhanced bioavailability claims could distinguish"
        }
        
        for keyword, opportunity in differentiators.items():
            if keyword in desc_lower:
                opportunities.append(opportunity)
        
        # Look for gaps in prior art
        if not any("dimer" in m.matching_keywords for m in matches):
            opportunities.append(
                "Dimer-based approach appears underexplored in prior art"
            )
        
        if not any("terpene" in m.matching_keywords for m in matches):
            opportunities.append(
                "Terpene synergy effects may offer novel claims"
            )
        
        return opportunities[:5]


async def quick_novelty_check(
    title: str,
    description: str,
    keywords: Optional[List[str]] = None
) -> Dict:
    """Quick novelty assessment (convenience function).
    
    Args:
        title: Innovation title
        description: Brief description
        keywords: Optional keywords
        
    Returns:
        Simplified novelty assessment dict
    """
    scorer = NoveltyScorer()
    report = await scorer.assess_novelty(title, description, keywords)
    
    return {
        "novelty_score": round(report.overall_novelty_score, 2),
        "novelty_level": report.novelty_level.value,
        "patentability": report.patentability_potential,
        "prior_art_found": report.prior_art_count,
        "high_concern_matches": report.high_relevance_count,
        "recommendation": report.recommendations[0] if report.recommendations else None
    }
