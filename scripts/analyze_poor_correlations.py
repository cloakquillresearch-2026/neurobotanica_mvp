import json
from collections import Counter
from pathlib import Path


def summarize_poor_correlations():
    payload = json.loads(Path("data/processed/training_correlations.json").read_text(encoding="utf-8"))
    correlations = payload.get("correlations", [])

    poor = [c for c in correlations if (c.get("quality") or "").lower() == "poor"]

    by_direction = Counter(c.get("direction") for c in poor)
    by_practice = Counter(
        (c.get("tk_practice_name") or c.get("tk_practice_predicted") or "UNKNOWN") for c in poor
    )
    by_target = Counter(c.get("genomic_target") for c in poor)
    by_dosage = Counter(c.get("dosage") for c in poor)
    low_conf = sorted(poor, key=lambda c: c.get("confidence", 0))[:15]

    return {
        "total_poor": len(poor),
        "by_direction": by_direction,
        "by_practice": by_practice,
        "by_target": by_target,
        "by_dosage": by_dosage,
        "low_confidence_examples": [
            {
                "correlation_id": c["correlation_id"],
                "confidence": round(float(c.get("confidence", 0)), 4),
                "tk_practice": c.get("tk_practice_name") or c.get("tk_practice_predicted"),
                "genomic_target": c.get("genomic_target"),
                "dosage": c.get("dosage"),
            }
            for c in low_conf
        ],
    }


def main():
    summary = summarize_poor_correlations()

    print(f"Poor correlations: {summary['total_poor']}")
    print("By direction:")
    for direction, count in summary["by_direction"].most_common():
        print(f"  {direction}: {count}")

    print("\nTop TK practices driving poor correlations:")
    for practice, count in summary["by_practice"].most_common(10):
        print(f"  {practice}: {count}")

    print("\nTop genomic targets represented:")
    for target, count in summary["by_target"].most_common(10):
        print(f"  {target}: {count}")

    print("\nDosage profile distribution:")
    for dosage, count in summary["by_dosage"].most_common():
        print(f"  {dosage}: {count}")

    print("\nLowest confidence poor correlations:")
    for corr in summary["low_confidence_examples"]:
        print(
            f"  {corr['correlation_id']} | conf={corr['confidence']:.3f} | "
            f"practice={corr['tk_practice']} | target={corr['genomic_target']} | "
            f"dosage={corr['dosage']}"
        )

    serializable = {
        "total_poor": summary["total_poor"],
        "by_direction": dict(summary["by_direction"]),
        "by_practice": dict(summary["by_practice"]),
        "by_target": dict(summary["by_target"]),
        "by_dosage": dict(summary["by_dosage"]),
        "low_confidence_examples": summary["low_confidence_examples"],
    }
    Path("data/processed/poor_correlation_triage.json").write_text(
        json.dumps(serializable, indent=2),
        encoding="utf-8",
    )
    print("\nWrote data/processed/poor_correlation_triage.json")


if __name__ == "__main__":
    main()
