"""
Prior Art Search - PatentPath Lite
Searches USPTO PatentsView API and EPO Open Patent Services for prior art.

Supports:
- Keyword/phrase search
- CPC/IPC classification filtering
- Date range filtering
- Cannabinoid-specific search templates
"""
import asyncio
import aiohttp
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from datetime import datetime, date
from enum import Enum
import logging
import json
import re
from urllib.parse import urlencode, quote

logger = logging.getLogger(__name__)


class PatentSource(Enum):
    """Patent database sources."""
    USPTO = "uspto"
    EPO = "epo"
    ALL = "all"


@dataclass
class PriorArtResult:
    """Single prior art patent/publication result."""
    patent_number: str
    title: str
    abstract: str
    publication_date: str
    inventors: List[str]
    assignees: List[str]
    classifications: List[str]  # CPC/IPC codes
    claims_count: int
    source: PatentSource
    relevance_score: float  # 0-1 relevance to query
    url: str
    
    # Additional metadata
    filing_date: Optional[str] = None
    priority_date: Optional[str] = None
    family_size: int = 0
    citation_count: int = 0
    
    def to_dict(self) -> Dict:
        return {
            "patent_number": self.patent_number,
            "title": self.title,
            "abstract": self.abstract[:500] + "..." if len(self.abstract) > 500 else self.abstract,
            "publication_date": self.publication_date,
            "inventors": self.inventors,
            "assignees": self.assignees,
            "classifications": self.classifications,
            "claims_count": self.claims_count,
            "source": self.source.value,
            "relevance_score": round(self.relevance_score, 3),
            "url": self.url,
            "filing_date": self.filing_date,
            "family_size": self.family_size,
            "citation_count": self.citation_count
        }


@dataclass
class SearchQuery:
    """Structured patent search query."""
    keywords: List[str]
    cpc_codes: List[str] = field(default_factory=list)
    ipc_codes: List[str] = field(default_factory=list)
    date_from: Optional[str] = None  # YYYY-MM-DD
    date_to: Optional[str] = None
    assignee: Optional[str] = None
    inventor: Optional[str] = None
    max_results: int = 50
    
    def to_query_string(self) -> str:
        """Convert to search string."""
        parts = []
        
        if self.keywords:
            parts.append(" AND ".join(f'"{kw}"' for kw in self.keywords))
        
        return " ".join(parts)


class USPTOClient:
    """USPTO PatentsView API client.
    
    API Documentation: https://patentsview.org/apis/api-endpoints
    Free API, no authentication required.
    """
    
    BASE_URL = "https://api.patentsview.org/patents/query"
    
    # Cannabinoid-related CPC codes
    CANNABIS_CPC_CODES = [
        "A61K31/352",  # Cannabinoids
        "A61K36/185",  # Cannabis (plant extracts)
        "C07D311/80",  # Benzopyrans (THC structure)
        "A61P25/00",   # CNS disorders
        "A61K31/05",   # Phenols (CBD)
    ]
    
    async def search(
        self,
        query: SearchQuery,
        timeout: int = 30
    ) -> List[PriorArtResult]:
        """Search USPTO PatentsView API."""
        results = []
        
        # Build query
        api_query = self._build_query(query)
        
        params = {
            "q": json.dumps(api_query),
            "f": json.dumps([
                "patent_number",
                "patent_title",
                "patent_abstract",
                "patent_date",
                "inventor_first_name",
                "inventor_last_name",
                "assignee_organization",
                "cpc_subgroup_id",
                "num_claims"
            ]),
            "o": json.dumps({"per_page": query.max_results})
        }
        
        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(
                    self.BASE_URL,
                    params=params,
                    timeout=aiohttp.ClientTimeout(total=timeout)
                ) as response:
                    if response.status == 200:
                        data = await response.json()
                        results = self._parse_results(data, query)
                    else:
                        logger.warning(f"USPTO API returned {response.status}")
                        
        except asyncio.TimeoutError:
            logger.error("USPTO API timeout")
        except Exception as e:
            logger.error(f"USPTO API error: {e}")
        
        return results
    
    def _build_query(self, query: SearchQuery) -> Dict:
        """Build PatentsView API query object."""
        conditions = []
        
        # Keyword search in title/abstract
        if query.keywords:
            keyword_conditions = []
            for kw in query.keywords:
                keyword_conditions.append({
                    "_or": [
                        {"_text_all": {"patent_title": kw}},
                        {"_text_all": {"patent_abstract": kw}}
                    ]
                })
            conditions.append({"_and": keyword_conditions})
        
        # CPC classification
        if query.cpc_codes:
            cpc_conditions = [
                {"_begins": {"cpc_subgroup_id": code}}
                for code in query.cpc_codes
            ]
            conditions.append({"_or": cpc_conditions})
        
        # Date range
        if query.date_from:
            conditions.append({
                "_gte": {"patent_date": query.date_from}
            })
        if query.date_to:
            conditions.append({
                "_lte": {"patent_date": query.date_to}
            })
        
        # Assignee
        if query.assignee:
            conditions.append({
                "_contains": {"assignee_organization": query.assignee}
            })
        
        if len(conditions) == 1:
            return conditions[0]
        elif len(conditions) > 1:
            return {"_and": conditions}
        else:
            # Default: cannabis-related patents
            return {
                "_or": [
                    {"_text_all": {"patent_abstract": "cannabinoid"}},
                    {"_text_all": {"patent_abstract": "cannabis"}},
                    {"_begins": {"cpc_subgroup_id": "A61K31/352"}}
                ]
            }
    
    def _parse_results(
        self,
        data: Dict,
        query: SearchQuery
    ) -> List[PriorArtResult]:
        """Parse USPTO API response."""
        results = []
        patents = data.get("patents", [])
        
        for patent in patents:
            try:
                # Parse inventors
                inventors = []
                for inv in patent.get("inventors", []):
                    name = f"{inv.get('inventor_first_name', '')} {inv.get('inventor_last_name', '')}".strip()
                    if name:
                        inventors.append(name)
                
                # Parse assignees
                assignees = [
                    a.get("assignee_organization", "")
                    for a in patent.get("assignees", [])
                    if a.get("assignee_organization")
                ]
                
                # Parse CPC codes
                classifications = list(set(
                    c.get("cpc_subgroup_id", "")
                    for c in patent.get("cpcs", [])
                    if c.get("cpc_subgroup_id")
                ))
                
                patent_number = patent.get("patent_number", "")
                
                result = PriorArtResult(
                    patent_number=patent_number,
                    title=patent.get("patent_title", ""),
                    abstract=patent.get("patent_abstract", ""),
                    publication_date=patent.get("patent_date", ""),
                    inventors=inventors[:5],  # Limit
                    assignees=assignees[:3],
                    classifications=classifications[:10],
                    claims_count=patent.get("num_claims", 0),
                    source=PatentSource.USPTO,
                    relevance_score=self._calculate_relevance(patent, query),
                    url=f"https://patents.google.com/patent/US{patent_number}"
                )
                results.append(result)
                
            except Exception as e:
                logger.warning(f"Error parsing patent: {e}")
                continue
        
        # Sort by relevance
        results.sort(key=lambda x: x.relevance_score, reverse=True)
        return results
    
    def _calculate_relevance(self, patent: Dict, query: SearchQuery) -> float:
        """Calculate relevance score based on keyword matches."""
        score = 0.0
        
        title = patent.get("patent_title", "").lower()
        abstract = patent.get("patent_abstract", "").lower()
        
        for keyword in query.keywords:
            kw = keyword.lower()
            # Title match (higher weight)
            if kw in title:
                score += 0.3
            # Abstract match
            if kw in abstract:
                score += 0.2
        
        # CPC match bonus
        patent_cpcs = [c.get("cpc_subgroup_id", "") for c in patent.get("cpcs", [])]
        for cpc in query.cpc_codes:
            for patent_cpc in patent_cpcs:
                if patent_cpc.startswith(cpc[:7]):  # Match at subclass level
                    score += 0.15
                    break
        
        return min(score, 1.0)


class EPOClient:
    """EPO Open Patent Services client (fallback/simulation).
    
    Note: Full EPO OPS requires registration. This provides
    a simulation for development.
    """
    
    BASE_URL = "https://ops.epo.org/3.2/rest-services"
    
    async def search(
        self,
        query: SearchQuery,
        timeout: int = 30
    ) -> List[PriorArtResult]:
        """Search EPO (simulated for MVP)."""
        # For MVP, return empty list
        # Full implementation requires EPO API registration
        logger.info("EPO search - requires API registration for production")
        return []


class PriorArtSearcher:
    """Main prior art search interface."""
    
    # Pre-defined search templates for cannabinoid research
    SEARCH_TEMPLATES = {
        "terpene_formulation": SearchQuery(
            keywords=["terpene", "cannabinoid", "formulation", "composition"],
            cpc_codes=["A61K31/352", "A61K36/185", "A61K47/00"],
            max_results=30
        ),
        "dimer_synthesis": SearchQuery(
            keywords=["cannabinoid", "dimer", "linked", "synthesis"],
            cpc_codes=["C07D311/80", "A61K31/352"],
            max_results=30
        ),
        "therapeutic_method": SearchQuery(
            keywords=["cannabinoid", "treatment", "method", "therapeutic"],
            cpc_codes=["A61K31/352", "A61P25/00", "A61P29/00"],
            max_results=30
        ),
        "extraction_process": SearchQuery(
            keywords=["cannabis", "extraction", "process", "method"],
            cpc_codes=["A61K36/185", "B01D11/02"],
            max_results=30
        ),
        "delivery_system": SearchQuery(
            keywords=["cannabinoid", "delivery", "controlled", "release"],
            cpc_codes=["A61K9/00", "A61K31/352"],
            max_results=30
        )
    }
    
    def __init__(self, sources: Optional[List[PatentSource]] = None):
        """Initialize searcher.
        
        Args:
            sources: Patent sources to search (default: USPTO only for MVP)
        """
        self.sources = sources or [PatentSource.USPTO]
        self.uspto = USPTOClient()
        self.epo = EPOClient()
    
    async def search(
        self,
        query: Optional[SearchQuery] = None,
        template: Optional[str] = None,
        keywords: Optional[List[str]] = None,
        cpc_codes: Optional[List[str]] = None,
        date_from: Optional[str] = None,
        date_to: Optional[str] = None,
        max_results: int = 50
    ) -> List[PriorArtResult]:
        """Search for prior art.
        
        Args:
            query: Pre-built SearchQuery object
            template: Name of search template to use
            keywords: Search keywords (alternative to query/template)
            cpc_codes: CPC codes to filter
            date_from: Start date YYYY-MM-DD
            date_to: End date YYYY-MM-DD
            max_results: Maximum results per source
            
        Returns:
            List of PriorArtResult sorted by relevance
        """
        # Build query
        if query is None:
            if template and template in self.SEARCH_TEMPLATES:
                query = self.SEARCH_TEMPLATES[template]
            elif keywords:
                query = SearchQuery(
                    keywords=keywords,
                    cpc_codes=cpc_codes or [],
                    date_from=date_from,
                    date_to=date_to,
                    max_results=max_results
                )
            else:
                raise ValueError("Must provide query, template, or keywords")
        
        all_results = []
        
        # Search each source
        tasks = []
        for source in self.sources:
            if source in [PatentSource.USPTO, PatentSource.ALL]:
                tasks.append(self.uspto.search(query))
            if source in [PatentSource.EPO, PatentSource.ALL]:
                tasks.append(self.epo.search(query))
        
        if tasks:
            results_lists = await asyncio.gather(*tasks, return_exceptions=True)
            
            for results in results_lists:
                if isinstance(results, list):
                    all_results.extend(results)
                elif isinstance(results, Exception):
                    logger.error(f"Search error: {results}")
        
        # Deduplicate and sort
        seen = set()
        unique_results = []
        for r in sorted(all_results, key=lambda x: x.relevance_score, reverse=True):
            if r.patent_number not in seen:
                seen.add(r.patent_number)
                unique_results.append(r)
        
        return unique_results[:max_results]
    
    async def search_for_compound(
        self,
        compound_name: str,
        smiles: Optional[str] = None,
        include_related: bool = True
    ) -> List[PriorArtResult]:
        """Search prior art for a specific compound.
        
        Args:
            compound_name: Name of compound (e.g., "THC", "CBD")
            smiles: Optional SMILES string
            include_related: Include structurally related compounds
            
        Returns:
            Prior art results related to the compound
        """
        keywords = [compound_name.lower()]
        
        # Add common synonyms
        synonyms = {
            "thc": ["tetrahydrocannabinol", "delta-9-thc", "Δ9-THC", "dronabinol"],
            "cbd": ["cannabidiol", "epidiolex"],
            "cbg": ["cannabigerol"],
            "cbn": ["cannabinol"],
            "thcv": ["tetrahydrocannabivarin"],
            "cbdv": ["cannabidivarin"]
        }
        
        name_lower = compound_name.lower()
        if name_lower in synonyms:
            keywords.extend(synonyms[name_lower])
        
        query = SearchQuery(
            keywords=keywords,
            cpc_codes=USPTOClient.CANNABIS_CPC_CODES[:3],
            max_results=30
        )
        
        return await self.search(query=query)
    
    async def search_for_application(
        self,
        therapeutic_area: str,
        formulation_type: Optional[str] = None
    ) -> List[PriorArtResult]:
        """Search prior art for therapeutic application.
        
        Args:
            therapeutic_area: Disease/condition (e.g., "pain", "epilepsy")
            formulation_type: Optional formulation type
            
        Returns:
            Prior art results
        """
        keywords = ["cannabinoid", therapeutic_area.lower()]
        
        if formulation_type:
            keywords.append(formulation_type.lower())
        
        # Map therapeutic areas to CPC codes
        therapeutic_cpc = {
            "pain": ["A61P29/00"],
            "epilepsy": ["A61P25/08"],
            "anxiety": ["A61P25/22"],
            "cancer": ["A61P35/00"],
            "inflammation": ["A61P29/00"],
            "nausea": ["A61P1/08"],
            "glaucoma": ["A61P27/06"]
        }
        
        cpc_codes = therapeutic_cpc.get(
            therapeutic_area.lower(),
            ["A61P25/00"]  # Default: CNS
        )
        
        query = SearchQuery(
            keywords=keywords,
            cpc_codes=cpc_codes + ["A61K31/352"],
            max_results=30
        )
        
        return await self.search(query=query)
    
    def get_template_names(self) -> List[str]:
        """Get available search template names."""
        return list(self.SEARCH_TEMPLATES.keys())
    
    def get_template(self, name: str) -> Optional[SearchQuery]:
        """Get a search template by name."""
        return self.SEARCH_TEMPLATES.get(name)
