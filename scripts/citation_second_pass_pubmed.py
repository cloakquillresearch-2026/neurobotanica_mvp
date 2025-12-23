"""Second-pass citation hardening using PubMed (NCBI eUtils).

Goal
----
Fill in unambiguous PubMed identifiers (PMID, DOI) into existing `citation` strings
without changing the JSON schema.

Safety
------
This script is intentionally conservative:
- It only updates a study when it finds exactly ONE PubMed record AND
  the returned PubMed title matches the study_title after normalization.
- If the title match is not exact, it skips.
- It never deletes existing citation text; it appends identifiers when missing.

Usage
-----
Dry run (report only):
  python scripts/citation_second_pass_pubmed.py

Apply updates:
  python scripts/citation_second_pass_pubmed.py --apply

Limit to specific conditions:
  python scripts/citation_second_pass_pubmed.py --apply --conditions TOURETTE_SYNDROME GLAUCOMA PARKINSONS

Notes
-----
NCBI recommends providing a tool name and email for API usage. Set:
  NCBI_TOOL and NCBI_EMAIL environment variables if you want.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import requests


DATA_DIR = Path("data") / "norml_extraction"
DEFAULT_GLOB = "*_studies.json"


def _normalize_title(title: str) -> str:
    title = title.strip().lower()
    title = re.sub(r"\s+", " ", title)
    # Remove punctuation that commonly varies in PubMed titles.
    title = re.sub(r"[\.,;:!\?\(\)\[\]\{\}\"'’]", "", title)
    title = title.replace("–", "-").replace("—", "-")
    title = re.sub(r"\s*-\s*", "-", title)
    return title


def _citation_has_pmid(citation: str) -> bool:
    return bool(re.search(r"\bPMID\s*:\s*\d+\b", citation, flags=re.IGNORECASE))


def _citation_has_doi(citation: str) -> bool:
    return bool(re.search(r"\bdoi\s*:\s*10\.\d{4,9}/\S+", citation, flags=re.IGNORECASE))


def _safe_get_str(obj: Dict[str, Any], key: str) -> str:
    value = obj.get(key)
    return value if isinstance(value, str) else ""


@dataclass
class PubMedHit:
    pmid: str
    title: str
    doi: Optional[str]


def _eutils_params() -> Dict[str, str]:
    params: Dict[str, str] = {
        "tool": os.environ.get("NCBI_TOOL", "neurobotanica_citation_second_pass"),
        "email": os.environ.get("NCBI_EMAIL", ""),
    }
    # Remove empty keys (NCBI accepts missing email).
    return {k: v for k, v in params.items() if v}


def _pubmed_esearch(term: str, timeout_s: float = 20.0) -> List[str]:
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pubmed",
        "retmode": "json",
        "retmax": "20",
        "term": term,
        **_eutils_params(),
    }
    resp = requests.get(base, params=params, timeout=timeout_s)
    resp.raise_for_status()
    data = resp.json()
    return data.get("esearchresult", {}).get("idlist", []) or []


def _pubmed_esummary_metadata(pmids: List[str], timeout_s: float = 20.0) -> Dict[str, Dict[str, str]]:
    if not pmids:
        return {}
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "pubmed",
        "retmode": "json",
        "id": ",".join(pmids),
        **_eutils_params(),
    }
    resp = requests.get(base, params=params, timeout=timeout_s)
    resp.raise_for_status()
    data = resp.json()
    result = data.get("result", {}) if isinstance(data, dict) else {}
    metadata: Dict[str, Dict[str, str]] = {}
    for pmid in pmids:
        entry = result.get(pmid)
        if not isinstance(entry, dict):
            continue
        title = (entry.get("title") or "").strip().rstrip(".")
        journal = (entry.get("fulljournalname") or entry.get("source") or "").strip()
        pubdate = (entry.get("pubdate") or entry.get("sortpubdate") or "").strip()
        year_match = re.search(r"(19|20)\d{2}", pubdate)
        year = year_match.group(0) if year_match else ""
        metadata[pmid] = {
            "title": title,
            "journal": journal,
            "year": year,
        }
    return metadata


def _pubmed_efetch_xml(pmid: str, timeout_s: float = 20.0) -> str:
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "pubmed",
        "retmode": "xml",
        "id": pmid,
        **_eutils_params(),
    }
    resp = requests.get(base, params=params, timeout=timeout_s)
    resp.raise_for_status()
    return resp.text


def _parse_efetch_for_title_and_doi(xml_text: str) -> Tuple[str, Optional[str]]:
    # Minimal XML extraction; avoids extra dependencies.
    # Title appears in <ArticleTitle> ... </ArticleTitle>
    title_match = re.search(r"<ArticleTitle>(.*?)</ArticleTitle>", xml_text, flags=re.DOTALL)
    title = "" if not title_match else re.sub(r"\s+", " ", title_match.group(1)).strip()

    # DOI appears in <ArticleId IdType=\"doi\">...</ArticleId>
    doi_match = re.search(r"<ArticleId\s+IdType=\"doi\">(.*?)</ArticleId>", xml_text)
    doi = None if not doi_match else doi_match.group(1).strip()

    # Strip basic XML entities if present.
    title = title.replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">")
    return title, doi


def _build_pubmed_query(study_title: str, year: str) -> str:
    # Title-first query; year is optional if missing.
    # Use [Title] field to keep specificity high.
    # Quote the title to reduce tokenization ambiguity.
    safe_title = study_title.replace('"', " ")
    # Use [ti] (title) tag for eUtils/PubMed fielding.
    title_term = f"\"{safe_title}\"[ti]"
    year_term = f"{year}[dp]" if year else ""
    return " AND ".join([t for t in [title_term, year_term] if t])


def _build_pubmed_metadata_queries(first_author: str, journal: str, year: str) -> List[str]:
    queries: List[str] = []
    if first_author and journal and year:
        safe_journal = journal.replace('"', " ")
        queries.append(f"{first_author}[Author] AND \"{safe_journal}\"[Jour] AND {year}[dp]")
    if first_author and year:
        queries.append(f"{first_author}[Author] AND {year}[dp]")
    if first_author and journal:
        safe_journal = journal.replace('"', " ")
        queries.append(f"{first_author}[Author] AND \"{safe_journal}\"[Jour]")
    return queries


def _first_author_surname(study: Dict[str, Any]) -> str:
    authors = _safe_get_str(study, "authors")
    if not authors:
        return ""
    # Heuristic: take the first comma-separated token, then the first word.
    first = authors.split(",", 1)[0].strip()
    if not first:
        return ""
    first_word = first.split()[0].strip()
    return re.sub(r"[^A-Za-z\-]", "", first_word)


def _normalize_journal(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def _extract_keywords(text: str) -> Set[str]:
    tokens = re.split(r"[^a-z0-9]+", text.lower())
    return {t for t in tokens if len(t) >= 4 and t not in {"study", "trial", "phase", "pilot", "open", "label"}}


def _select_candidate_by_score(
    meta: Dict[str, Dict[str, str]],
    study_title: str,
    expected_journal: str,
    expected_year: str,
) -> Optional[str]:
    keywords = _extract_keywords(study_title)
    scored: List[Tuple[str, int, bool, bool]] = []
    for pmid, info in meta.items():
        cand_title = info.get("title", "")
        cand_keywords = _extract_keywords(cand_title)
        overlap = len(keywords & cand_keywords) if keywords else 0
        journal_match = bool(expected_journal) and _normalize_journal(info.get("journal", "")) == expected_journal
        year_match = bool(expected_year) and info.get("year", "") == expected_year
        scored.append((pmid, overlap, journal_match, year_match))

    if not scored:
        return None

    def priority(entry: Tuple[str, int, bool, bool]) -> Tuple[int, int, int]:
        _, overlap, journal_match, year_match = entry
        return (1 if journal_match else 0, 1 if year_match else 0, overlap)

    scored.sort(key=priority, reverse=True)
    top = scored[0]
    top_overlap = top[1]
    top_journal = top[2]
    top_year = top[3]

    min_required = 2
    if top_journal or top_year:
        min_required = 1

    if top_overlap < min_required:
        return None

    if len(scored) > 1 and priority(scored[1]) == priority(top):
        return None

    return top[0]


def _extract_year(study: Dict[str, Any]) -> str:
    year = _safe_get_str(study, "year")
    if year and re.fullmatch(r"\d{4}", year):
        return year

    citation = _safe_get_str(study, "citation")
    m = re.search(r"\b(19\d{2}|20\d{2})\b", citation)
    return m.group(1) if m else ""


def _append_identifiers(citation: str, pmid: Optional[str], doi: Optional[str]) -> str:
    updated = citation.strip()

    # Normalize trailing punctuation/spacing.
    updated = re.sub(r"\s+", " ", updated)
    updated = updated.rstrip(" .;")

    parts: List[str] = [updated] if updated else []

    if pmid and not _citation_has_pmid(updated):
        parts.append(f"PMID: {pmid}")
    if doi and not _citation_has_doi(updated):
        parts.append(f"doi: {doi}")

    return "; ".join(parts) + "." if parts else updated


def _load_condition_file(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _write_condition_file(path: Path, payload: Dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
        f.write("\n")


def _iter_condition_files() -> Iterable[Path]:
    yield from sorted(DATA_DIR.glob(DEFAULT_GLOB))


def _condition_name_from_payload(payload: Dict[str, Any], fallback: str) -> str:
    cond = payload.get("condition")
    return cond if isinstance(cond, str) and cond.strip() else fallback


def _fallback_condition_from_filename(path: Path) -> str:
    # e.g. tourette_syndrome_studies.json -> TOURETTE_SYNDROME
    name = path.name.replace("_studies.json", "").upper()
    return name


def main() -> int:
    parser = argparse.ArgumentParser(description="Second-pass citation hardening via PubMed")
    parser.add_argument("--apply", action="store_true", help="Write updates back to files")
    parser.add_argument(
        "--conditions",
        nargs="*",
        default=None,
        help="Optional list of condition names to limit (e.g., TOURETTE_SYNDROME GLAUCOMA PARKINSONS)",
    )
    parser.add_argument("--sleep", type=float, default=0.34, help="Seconds to sleep between NCBI requests")
    parser.add_argument(
        "--debug",
        type=int,
        default=0,
        help="Print up to N sample queries and their hit counts (dry-run style visibility)",
    )
    args = parser.parse_args()

    allowed_conditions = set([c.upper() for c in args.conditions]) if args.conditions else None

    total_candidates = 0
    total_updates = 0
    total_skipped_ambiguous = 0
    total_skipped_mismatch = 0
    total_resolved_by_title = 0
    total_resolved_by_metadata = 0
    debug_remaining = max(0, int(args.debug))

    for path in _iter_condition_files():
        payload = _load_condition_file(path)
        condition = _condition_name_from_payload(payload, _fallback_condition_from_filename(path))
        condition_upper = condition.upper()

        if allowed_conditions is not None and condition_upper not in allowed_conditions:
            continue

        studies = payload.get("studies")
        if not isinstance(studies, list):
            continue

        file_updates = 0

        for study in studies:
            if not isinstance(study, dict):
                continue

            citation = _safe_get_str(study, "citation")
            study_title = _safe_get_str(study, "study_title") or _safe_get_str(study, "title")
            if not study_title:
                continue

            # Skip if already has both identifiers.
            if _citation_has_pmid(citation) and _citation_has_doi(citation):
                continue

            total_candidates += 1
            year = _extract_year(study)
            surname = _first_author_surname(study)
            journal = _safe_get_str(study, "journal")
            expected_journal = _normalize_journal(journal) if journal else ""

            attempts: List[Tuple[str, str]] = []
            title_query = _build_pubmed_query(study_title, year)
            if surname:
                title_query = f"({title_query}) AND {surname}[Author]"
            if title_query:
                attempts.append(("title", title_query))

            for metadata_query in _build_pubmed_metadata_queries(surname, journal, year):
                attempts.append(("metadata", metadata_query))

            pmid: Optional[str] = None
            matched_meta: Dict[str, str] = {}
            pmid_source = ""

            for mode, query in attempts:
                try:
                    ids = _pubmed_esearch(query)
                except Exception:
                    continue

                time.sleep(args.sleep)

                if debug_remaining > 0:
                    sid = _safe_get_str(study, "study_id")
                    print(f"DEBUG query_hits={len(ids)} study_id={sid} mode={mode} query={query}")
                    debug_remaining -= 1

                if not ids:
                    continue

                shortlist = ids[:10]
                try:
                    meta = _pubmed_esummary_metadata(shortlist)
                except Exception:
                    meta = {}

                time.sleep(args.sleep)

                if mode == "title":
                    matches = [
                        candidate
                        for candidate, candidate_meta in meta.items()
                        if _normalize_title(candidate_meta.get("title", "")) == _normalize_title(study_title)
                    ]

                    if len(matches) == 1:
                        pmid = matches[0]
                        matched_meta = meta.get(pmid, {})
                        pmid_source = "title"
                        total_resolved_by_title += 1
                        break
                    elif len(matches) > 1:
                        total_skipped_ambiguous += 1
                        continue
                    else:
                        continue
                else:
                    strict_matches = []
                    for candidate, candidate_meta in meta.items():
                        cand_journal = _normalize_journal(candidate_meta.get("journal", ""))
                        cand_year = candidate_meta.get("year", "")
                        if expected_journal and expected_journal != cand_journal:
                            continue
                        if year and cand_year and cand_year != year:
                            continue
                        strict_matches.append(candidate)

                    if len(strict_matches) == 1:
                        pmid = strict_matches[0]
                        matched_meta = meta.get(pmid, {})
                        pmid_source = "metadata_strict"
                        total_resolved_by_metadata += 1
                        break
                    elif len(strict_matches) > 1:
                        total_skipped_ambiguous += 1
                        continue

                    scored_candidate = _select_candidate_by_score(meta, study_title, expected_journal, year)
                    if scored_candidate:
                        pmid = scored_candidate
                        matched_meta = meta.get(pmid, {})
                        pmid_source = "metadata_score"
                        total_resolved_by_metadata += 1
                        break
                    else:
                        continue

            if not pmid:
                total_skipped_ambiguous += 1
                continue

            try:
                xml_text = _pubmed_efetch_xml(pmid)
            except Exception:
                continue

            time.sleep(args.sleep)

            pub_title, doi = _parse_efetch_for_title_and_doi(xml_text)
            if not pub_title:
                continue

            if _normalize_title(pub_title) != _normalize_title(study_title):
                if not pmid_source.startswith("metadata"):
                    total_skipped_mismatch += 1
                    continue

            updated_citation = _append_identifiers(citation, pmid=pmid, doi=doi)
            if updated_citation != citation:
                study["citation"] = updated_citation
                file_updates += 1
                total_updates += 1

        if file_updates and args.apply:
            _write_condition_file(path, payload)

        if file_updates:
            print(f"{path.as_posix()}: updated {file_updates} citations")

    print("\nSummary")
    print("-------")
    print(f"candidates_checked: {total_candidates}")
    print(f"citations_updated:  {total_updates}")
    print(f"skipped_ambiguous: {total_skipped_ambiguous}")
    print(f"skipped_mismatch:  {total_skipped_mismatch}")
    print(f"resolved_by_title: {total_resolved_by_title}")
    print(f"resolved_by_metadata: {total_resolved_by_metadata}")
    print(f"mode: {'APPLY' if args.apply else 'DRY_RUN'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
