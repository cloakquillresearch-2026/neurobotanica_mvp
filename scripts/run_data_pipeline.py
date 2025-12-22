"""
GenomePath Data Pipeline - Master Orchestrator
==============================================

Runs the complete data pipeline from extraction to validation.

Usage:
    python scripts/run_data_pipeline.py

Pipeline Steps:
    1. Extract genomic targets from NORML studies
    2. Build TK practice dataset from published sources
    3. Generate TK↔Genomic correlations using trade secret algorithms
    4. Validate all datasets for quality and compliance

Outputs:
    data/processed/genomic_targets.json
    data/processed/tk_practices.json
    data/processed/training_correlations.json
    data/processed/validation_report.json
"""

import sys
import subprocess
from pathlib import Path
from datetime import datetime


def run_script(script_name: str, description: str) -> bool:
    """Run a pipeline script and return success status."""
    
    print()
    print("=" * 70)
    print(f"STEP: {description}")
    print("=" * 70)
    print()
    
    script_path = Path(__file__).parent / script_name
    
    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            check=True,
            capture_output=False,
            text=True
        )
        
        print()
        print(f"✅ {description} - COMPLETE")
        return True
        
    except subprocess.CalledProcessError as e:
        print()
        print(f"❌ {description} - FAILED")
        print(f"Error: {e}")
        return False
    except Exception as e:
        print()
        print(f"❌ {description} - ERROR")
        print(f"Unexpected error: {e}")
        return False


def main():
    """Run complete data pipeline."""
    
    start_time = datetime.now()
    
    print()
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 68 + "║")
    print("║" + "  GenomePath TS-GP-001 Data Pipeline".center(68) + "║")
    print("║" + "  Complete TK↔Genomic Training Data Generation".center(68) + "║")
    print("║" + " " * 68 + "║")
    print("╚" + "=" * 68 + "╝")
    print()
    print(f"Started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Pipeline steps
    steps = [
        ("extract_genomic_targets.py", "Extract Genomic Targets from NORML Studies"),
        ("build_tk_dataset.py", "Build TK Practice Dataset"),
        ("generate_correlations.py", "Generate TK↔Genomic Correlations"),
        ("validate_dataset.py", "Validate Complete Dataset")
    ]
    
    results = []
    for script, description in steps:
        success = run_script(script, description)
        results.append((description, success))
        
        if not success:
            print()
            print("=" * 70)
            print("⚠️  Pipeline stopped due to error")
            print("=" * 70)
            break
    
    # Print final summary
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    
    print()
    print()
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 68 + "║")
    print("║" + "  PIPELINE SUMMARY".center(68) + "║")
    print("║" + " " * 68 + "║")
    print("╚" + "=" * 68 + "╝")
    print()
    
    for description, success in results:
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"{status}  {description}")
    
    print()
    print(f"Completed: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Duration: {duration:.1f} seconds")
    print()
    
    # Final status
    all_success = all(success for _, success in results)
    
    if all_success:
        print("╔" + "=" * 68 + "╗")
        print("║" + " " * 68 + "║")
        print("║" + "  🎉 PIPELINE COMPLETE - ALL STEPS PASSED".center(76) + "║")
        print("║" + " " * 68 + "║")
        print("╚" + "=" * 68 + "╝")
        print()
        print("Generated Datasets:")
        print("  📊 data/processed/genomic_targets.json")
        print("  📚 data/processed/tk_practices.json")
        print("  🔗 data/processed/training_correlations.json")
        print("  ✅ data/processed/validation_report.json")
        print()
        print("Next Steps:")
        print("  1. Review validation_report.json for any warnings")
        print("  2. Manually validate top 50 correlations against literature")
        print("  3. Use tk_practice_template.json to add more TK practices")
        print("  4. Re-run pipeline after adding new data")
        print()
        return 0
    else:
        print("╔" + "=" * 68 + "╗")
        print("║" + " " * 68 + "║")
        print("║" + "  ❌ PIPELINE FAILED - SEE ERRORS ABOVE".center(76) + "║")
        print("║" + " " * 68 + "║")
        print("╚" + "=" * 68 + "╝")
        print()
        return 1


if __name__ == "__main__":
    exit(main())
