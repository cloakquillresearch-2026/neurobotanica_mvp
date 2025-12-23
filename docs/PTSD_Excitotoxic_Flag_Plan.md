# PTSD Excitotoxic Flag Plan

**Last updated:** 2025-12-23

## Objective
Translate the glaucoma mechanistic evidence (Nucci et al. 2008, study ID `GLAUCOMA_MECHANISTIC_002`) into a reusable metadata flag that documents excitotoxic shielding for PTSD practices and enforces that requirement during validation.

## Data Sources
- `data/norml_extraction/glaucoma_studies.json#GLAUCOMA_MECHANISTIC_002` — retinal ganglion cell protection from glutamate excitotoxicity (40% improved survival).
- `data/norml_extraction/ptsd_studies.json` — literature set used during PTSD review and validation packs.

## Schema Additions
- `neuroprotective_metadata` (per TK practice)
  - `excitotoxic_flag` (bool) — true once a practice has vetted mechanistic evidence.
  - `excitotoxic_evidence_source` (string) — NORML or internal study ID providing the mechanistic citation.
  - `excitotoxic_reference_description` (string) — <= 240 char human summary.
  - `translation_notes` (string) — how the mechanism maps into PTSD neural circuitry.
  - `last_verified` (ISO date) — audit trail for regulatory submissions.

These keys were added to `data/processed/tk_practice_template.json` and populated for `TK_0013` (Cannabis for PTSD and trauma) inside `data/processed/tk_practices.json`.

## Validation Guardrail
`scripts/validate_dataset.py` now includes `validate_ptsd_excitotoxic_flag()`:
- Scans all TK practices whose `indications` include `ptsd` (case insensitive).
- Fails validation if `neuroprotective_metadata.excitotoxic_flag` or `excitotoxic_evidence_source` is missing.
- Surfaces violations inside the exported `excitotoxic_checks` block so auditors see provenance.

## Workflow Impact
1. Curators attach the excitotoxic metadata whenever PTSD-related practices are created or updated.
2. Dataset validation must pass the excitotoxic check before correlations are released to modeling teams.
3. Future PTSD practices can reference additional mechanistic sources (e.g., hippocampal or amygdala excitotoxic models) by adding more evidence fields; validator logic already enforces at least one citation.

This plan keeps NeuroBotanica aligned with the patent claims around excitotoxic mitigation while ensuring every PTSD correlation run is backed by explicit, reviewable metadata.
