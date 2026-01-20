# Poor Correlation Triage – 23 Dec 2025

## Snapshot
- **Total poor correlations:** 114 (all TK → Genomic)
- **Dosage mix:** 52 high / 41 medium / 21 low
- **Top TK practices triggering low scores:**
  1. Hemp seed for lactation support (13)
  2. Cannabis for nausea and vomiting (9)
  3. Cannabis salve for burns (9)
  4. Hemp for hypertension (7)
  5. Hemp seed for appetite stimulation (6)
  6. Cannabis smoke for asthma (6)
  7. Hemp protein for muscle recovery (6)
  8. Hemp fiber tea for diabetes (6)
  9. Hemp seed oil for skin conditions (5)
  10. Cannabis for chemotherapy side effects (5)
- **Common genomic targets:** PPAR-γ (11), TRPV1 (9), AMPK/mTOR/5-HT3 (6 each), ghrelin & D2 receptors (5 each).
- **Confidence floor:** worst offenders fall between 0.37–0.45 despite hardened sources.
- **Hardened sources:** 100% of the poor set already include NORML evidence—root issue is low model confidence, not missing citations.

## Root Causes
1. **Dosage fan-out on weak base scores** – The generator emits low/medium/high variants even when the baseline `overall_confidence` < 0.45, so high-dose rows inherit penalties (−0.03 modifier) and sink into "poor".
2. **Practice novelty vs. mechanistic coverage** – Women’s health (lactation) and metabolic hemp practices have sparse mechanistic descriptors (tissue_confidence = 0), so bridge vectors underweight them for PPAR, oxytocin, prolactin targets.
3. **Target over-assignment** – PPAR-γ / TRPV1 / AMPK are currently assigned to a broad set of practices without differentiating tissue context, diluting confidences.

## Recommended Remediations
1. **Dosage guardrail**
   - Skip generating medium/high dose variants when the base hypothesis confidence < 0.48.
   - Alternatively clamp `confidence_modifier` to 0 for weak hypotheses so we do not create artificial low rows.
2. **Practice-specific evidence floor**
   - Tag high-risk practices (`TK_0027` lactation, `TK_0016` hypertension, etc.) with a `requires_additional_mechanistic_evidence` flag in `tk_practices.json` and suppress their correlations until tissue/expression metadata is enriched.
3. **Target curation**
   - Revisit mapping heuristics for PPAR-γ/AMPK/TRPV1 and restrict to practices whose indication sets overlap with validated NORML conditions (arthritis, chronic pain, neuropathy). This should remove ~35 low-signal entries immediately.
4. **Tissue confidence floor**
   - In `TKGenomicCorrelator`, require `tissue_confidence >= 0.2` before assigning receptor/pathway targets, otherwise downgrade to `requires_validation` instead of surfacing as training data.
5. **Follow-up QA loop**
   - After the above fixes, rerun `scripts/generate_correlations.py` + `scripts/post_harden_training_correlations.py` and monitor `data/processed/poor_correlation_triage.json` to ensure the poor bucket drops below 5% of total correlations.

## Artefacts
- Machine-readable stats: `data/processed/poor_correlation_triage.json`
- Generator summary helper: `scripts/analyze_poor_correlations.py`
