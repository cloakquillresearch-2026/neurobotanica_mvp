"""Generate black-and-white patent figures for NeuroBotanica provisional.

Outputs:
- docs/figures/Figure_01_System_Architecture.svg (+ .png)
- ... through Figure_08_Continuous_Learning.svg (+ .png)

Design goals:
- Patent-friendly: monochrome line drawings, no color dependencies.
- Vector-first: SVG for crisp scaling; PNG for easy embedding.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle


@dataclass(frozen=True)
class Box:
    x: float
    y: float
    w: float
    h: float
    text: str


def _setup_ax(figsize=(11, 8.5)):
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    return fig, ax


def _add_box(ax, box: Box, fontsize=10, lw=1.5):
    rect = Rectangle(
        (box.x, box.y),
        box.w,
        box.h,
        linewidth=lw,
        edgecolor="black",
        facecolor="white",
    )
    ax.add_patch(rect)
    ax.text(
        box.x + box.w / 2,
        box.y + box.h / 2,
        box.text,
        ha="center",
        va="center",
        fontsize=fontsize,
        family="DejaVu Sans",
        wrap=True,
    )


def _arrow(ax, x1, y1, x2, y2, lw=1.3, style="-|>"):
    arr = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle=style,
        mutation_scale=12,
        linewidth=lw,
        color="black",
        shrinkA=6,
        shrinkB=6,
    )
    ax.add_patch(arr)


def _title(ax, title: str):
    ax.text(
        0.5,
        0.98,
        title,
        ha="center",
        va="top",
        fontsize=13,
        family="DejaVu Sans",
        fontweight="bold",
    )


def _save(fig, out_base: Path):
    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure_01(out_dir: Path):
    fig, ax = _setup_ax()
    _title(ax, "Figure 1. High-Level System Architecture")

    infra = Box(0.05, 0.62, 0.90, 0.28, "CONSENT-GATED FOUNDATIONAL INFRASTRUCTURE")
    proc = Box(0.05, 0.33, 0.90, 0.24, "ANALYTICAL PROCESSING FRAMEWORK")
    pred = Box(0.05, 0.08, 0.90, 0.20, "DIMERIC PREDICTION METHODOLOGY")

    for b in (infra, proc, pred):
        _add_box(ax, b, fontsize=11)

    # Infrastructure submodules
    infra_boxes = [
        Box(0.08, 0.66, 0.21, 0.18, "Policy\nOrchestration"),
        Box(0.31, 0.66, 0.21, 0.18, "Consent\nCompilation"),
        Box(0.54, 0.66, 0.21, 0.18, "Attribution\nRegistry"),
        Box(0.77, 0.66, 0.16, 0.18, "Benefit-\nSharing\nExecution"),
    ]
    for b in infra_boxes:
        _add_box(ax, b, fontsize=9)

    # Processing submodules
    proc_boxes = [
        Box(0.08, 0.38, 0.27, 0.14, "Data Ingestion\n(studies, evidence, structures)"),
        Box(0.38, 0.38, 0.27, 0.14, "Molecular Processing\n(conformers, dimers)"),
        Box(0.68, 0.38, 0.25, 0.14, "Regulatory Outputs\n(FDA Schedule III docs)"),
    ]
    for b in proc_boxes:
        _add_box(ax, b, fontsize=9)

    # Prediction submodules
    pred_boxes = [
        Box(0.08, 0.11, 0.27, 0.12, "Reactive Site\nIdentification"),
        Box(0.38, 0.11, 0.27, 0.12, "Geometric Compatibility\n+ Oxidative Modeling"),
        Box(0.68, 0.11, 0.25, 0.12, "Formation Probability\n+ Ranked Outputs"),
    ]
    for b in pred_boxes:
        _add_box(ax, b, fontsize=9)

    # Arrows between layers
    _arrow(ax, 0.50, 0.62, 0.50, 0.57)
    _arrow(ax, 0.50, 0.33, 0.50, 0.28)

    out_base = out_dir / "Figure_01_System_Architecture"
    _save(fig, out_base)


def figure_02(out_dir: Path):
    fig, ax = _setup_ax()
    _title(ax, "Figure 2. Dimeric Structure Prediction Workflow")

    steps = [
        Box(0.07, 0.82, 0.86, 0.10, "Input: Parent cannabinoid (SMILES/InChI/SDF) + Oxidative aging parameters"),
        Box(0.07, 0.69, 0.86, 0.10, "Validate structure + normalize representations"),
        Box(0.07, 0.56, 0.86, 0.10, "Graph-based reactive site identification"),
        Box(0.07, 0.43, 0.86, 0.10, "3D conformer generation + geometric compatibility assessment"),
        Box(0.07, 0.30, 0.86, 0.10, "Oxidative environment modeling (temp, O2, light, duration)"),
        Box(0.07, 0.17, 0.86, 0.10, "ML refinement + uncertainty quantification"),
        Box(0.07, 0.04, 0.86, 0.10, "Output: Predicted dimers + formation probability scores + ranked candidates"),
    ]
    for b in steps:
        _add_box(ax, b, fontsize=10)

    for i in range(len(steps) - 1):
        _arrow(ax, 0.50, steps[i].y, 0.50, steps[i + 1].y + steps[i + 1].h)

    out_base = out_dir / "Figure_02_Dimer_Prediction_Workflow"
    _save(fig, out_base)


def figure_03(out_dir: Path):
    fig, ax = _setup_ax(figsize=(11, 6.5))
    _title(ax, "Figure 3. Traditional Knowledge Consent Validation Integration")

    tk = Box(0.05, 0.55, 0.28, 0.25, "Traditional Knowledge\nPreparation Method\n(Inputs)")
    detector = Box(0.36, 0.60, 0.28, 0.15, "TK Detection\n+ Metadata")
    consent = Box(0.67, 0.60, 0.28, 0.15, "Consent Compilation\nService")

    decision_allow = Box(0.67, 0.35, 0.13, 0.15, "ALLOW")
    decision_deny = Box(0.82, 0.35, 0.13, 0.15, "DENY")

    analysis = Box(0.36, 0.22, 0.59, 0.18, "Authorized Analysis\n(Dimer prediction / linkage / docs)")

    for b in (tk, detector, consent, decision_allow, decision_deny, analysis):
        _add_box(ax, b, fontsize=10)

    _arrow(ax, tk.x + tk.w, tk.y + tk.h / 2, detector.x, detector.y + detector.h / 2)
    _arrow(ax, detector.x + detector.w, detector.y + detector.h / 2, consent.x, consent.y + consent.h / 2)

    _arrow(ax, consent.x + consent.w / 2, consent.y, decision_allow.x + decision_allow.w / 2, decision_allow.y + decision_allow.h)
    _arrow(ax, consent.x + consent.w / 2, consent.y, decision_deny.x + decision_deny.w / 2, decision_deny.y + decision_deny.h)

    _arrow(ax, decision_allow.x + decision_allow.w / 2, decision_allow.y, analysis.x + analysis.w * 0.35, analysis.y + analysis.h)

    ax.text(0.86, 0.18, "If DENY:\nReturn proof of denial\n+ halt propagation", ha="center", va="center", fontsize=9)

    out_base = out_dir / "Figure_03_Consent_Validation"
    _save(fig, out_base)


def figure_04(out_dir: Path):
    fig, ax = _setup_ax()
    _title(ax, "Figure 4. Traditional Knowledge Linkage Methodology")

    b1 = Box(0.05, 0.75, 0.40, 0.16, "Preparation method description\n(text / structured fields)")
    b2 = Box(0.55, 0.75, 0.40, 0.16, "Parameter extraction\n(aging duration, temp, oxygen, light, curing)")

    b3 = Box(0.05, 0.50, 0.90, 0.16, "Map extracted parameters to oxidative model inputs")
    b4 = Box(0.05, 0.28, 0.90, 0.16, "Run dimeric formation models + generate predicted enrichment")
    b5 = Box(0.05, 0.06, 0.90, 0.16, "Output: linked methods → prioritized dimer candidates + attribution metadata")

    for b in (b1, b2, b3, b4, b5):
        _add_box(ax, b, fontsize=10)

    _arrow(ax, b1.x + b1.w, b1.y + b1.h / 2, b2.x, b2.y + b2.h / 2)
    _arrow(ax, 0.50, b2.y, 0.50, b3.y + b3.h)
    _arrow(ax, 0.50, b3.y, 0.50, b4.y + b4.h)
    _arrow(ax, 0.50, b4.y, 0.50, b5.y + b5.h)

    out_base = out_dir / "Figure_04_TK_Linkage"
    _save(fig, out_base)


def figure_05(out_dir: Path):
    fig, ax = _setup_ax()
    _title(ax, "Figure 5. Automated Benefit-Sharing Execution")

    manifest = Box(0.05, 0.76, 0.90, 0.14, "Cryptographically signed attribution manifest\n(contributions, weights, identities)")
    scoring = Box(0.05, 0.58, 0.90, 0.14, "Proportional contribution scoring + payment calculation")
    contract = Box(0.05, 0.40, 0.90, 0.14, "Smart contract generation + execution")
    ledger = Box(0.05, 0.22, 0.90, 0.14, "Immutable ledger update (audit trail + receipts)")
    dist = Box(0.05, 0.04, 0.90, 0.14, "Distribution to source communities\n(automated transfers + verification)")

    for b in (manifest, scoring, contract, ledger, dist):
        _add_box(ax, b, fontsize=10)

    for a, b in ((manifest, scoring), (scoring, contract), (contract, ledger), (ledger, dist)):
        _arrow(ax, 0.50, a.y, 0.50, b.y + b.h)

    out_base = out_dir / "Figure_05_Benefit_Sharing"
    _save(fig, out_base)


def figure_06(out_dir: Path):
    fig, ax = _setup_ax()
    _title(ax, "Figure 6. FDA Schedule III Compliance Documentation Generation")

    inp = Box(0.05, 0.78, 0.90, 0.14, "Inputs: predicted dimer structure + evidence + analytics")
    cmc = Box(0.05, 0.55, 0.28, 0.18, "CMC Package\n(identity, purity, stability, specs)")
    pharm = Box(0.36, 0.55, 0.28, 0.18, "Pharmacology\n(MoA, targets, PK/PD)")
    clin = Box(0.67, 0.55, 0.28, 0.18, "Clinical Evidence\n(mapping + citations)")

    out = Box(0.05, 0.25, 0.90, 0.20, "Output: submission-ready documentation\n(modular templates for Schedule III FDA pathways)")

    for b in (inp, cmc, pharm, clin, out):
        _add_box(ax, b, fontsize=10)

    _arrow(ax, 0.50, inp.y, 0.20, cmc.y + cmc.h)
    _arrow(ax, 0.50, inp.y, 0.50, pharm.y + pharm.h)
    _arrow(ax, 0.50, inp.y, 0.80, clin.y + clin.h)

    _arrow(ax, cmc.x + cmc.w / 2, cmc.y, out.x + out.w * 0.25, out.y + out.h)
    _arrow(ax, pharm.x + pharm.w / 2, pharm.y, out.x + out.w * 0.50, out.y + out.h)
    _arrow(ax, clin.x + clin.w / 2, clin.y, out.x + out.w * 0.75, out.y + out.h)

    out_base = out_dir / "Figure_06_FDA_Documentation"
    _save(fig, out_base)


def figure_07(out_dir: Path):
    fig, ax = _setup_ax()
    _title(ax, "Figure 7. Modular Analytical Processing Architecture")

    # Left: disclosed methodology
    disclosed = Box(0.05, 0.18, 0.42, 0.68, "Disclosed Methodology Layer\n(Patent-protected workflows)")
    _add_box(ax, disclosed, fontsize=11)

    disclosed_boxes = [
        Box(0.08, 0.70, 0.36, 0.12, "Workflow Orchestration"),
        Box(0.08, 0.54, 0.36, 0.12, "Consent / Attribution Hooks"),
        Box(0.08, 0.38, 0.36, 0.12, "Prediction Pipeline Interfaces"),
        Box(0.08, 0.22, 0.36, 0.12, "Regulatory Output Interfaces"),
    ]
    for b in disclosed_boxes:
        _add_box(ax, b, fontsize=9)

    # Right: proprietary implementations
    proprietary = Box(0.53, 0.18, 0.42, 0.68, "Implementation Layer\n(Trade secrets accessed via APIs)")
    _add_box(ax, proprietary, fontsize=11)

    prop_boxes = [
        Box(0.56, 0.70, 0.36, 0.12, "Formation Probability Engine"),
        Box(0.56, 0.54, 0.36, 0.12, "ML Models + Calibration"),
        Box(0.56, 0.38, 0.36, 0.12, "Template / Report Generator"),
        Box(0.56, 0.22, 0.36, 0.12, "Optimization + Scaling"),
    ]
    for b in prop_boxes:
        _add_box(ax, b, fontsize=9)

    # API boundary arrows
    for y in (0.76, 0.60, 0.44, 0.28):
        _arrow(ax, 0.47, y, 0.53, y)
        ax.text(0.50, y + 0.03, "API", ha="center", va="bottom", fontsize=9)

    out_base = out_dir / "Figure_07_Modular_Architecture"
    _save(fig, out_base)


def figure_08(out_dir: Path):
    fig, ax = _setup_ax()
    _title(ax, "Figure 8. Continuous Learning and Model Refinement")

    pred = Box(0.08, 0.70, 0.84, 0.16, "Generate predictions + rank candidates\n(probability + uncertainty)")
    synth = Box(0.08, 0.50, 0.84, 0.16, "Experimental validation\n(synthesis + characterization)")
    data = Box(0.08, 0.30, 0.84, 0.16, "Validation dataset ingestion\n(ground truth outcomes)")
    retrain = Box(0.08, 0.10, 0.84, 0.16, "Automated retraining trigger\n(update models; preserve stability)")

    for b in (pred, synth, data, retrain):
        _add_box(ax, b, fontsize=10)

    _arrow(ax, 0.50, pred.y, 0.50, synth.y + synth.h)
    _arrow(ax, 0.50, synth.y, 0.50, data.y + data.h)
    _arrow(ax, 0.50, data.y, 0.50, retrain.y + retrain.h)

    # Loop-back
    _arrow(ax, 0.92, 0.18, 0.92, 0.78, style="-|>")
    _arrow(ax, 0.92, 0.78, 0.80, 0.78, style="-|>")
    ax.text(0.92, 0.50, "Feedback\nloop", ha="center", va="center", rotation=90, fontsize=9)

    # Active learning annotation
    ax.text(0.50, 0.02, "Active learning selects next experiments based on uncertainty / value-of-information", ha="center", va="bottom", fontsize=9)

    out_base = out_dir / "Figure_08_Continuous_Learning"
    _save(fig, out_base)


def main():
    repo_root = Path(__file__).resolve().parents[1]
    out_dir = repo_root / "docs" / "figures"

    figure_01(out_dir)
    figure_02(out_dir)
    figure_03(out_dir)
    figure_04(out_dir)
    figure_05(out_dir)
    figure_06(out_dir)
    figure_07(out_dir)
    figure_08(out_dir)

    outputs = sorted(out_dir.glob("Figure_*.svg"))
    print(f"Generated {len(outputs)} SVG figures in: {out_dir}")


if __name__ == "__main__":
    main()
