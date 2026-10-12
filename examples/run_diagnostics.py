"""
run_diagnostics.py - is our scoring function worth trusting?

Runs the six checks on our own configuration, then tries to build a wet-lab
proposal sheet.

    python examples/run_diagnostics.py

On the config as shipped the checks FAIL and the sheet refuses to issue
anything. That is the intended demonstration, not a broken example. We left
it this way because a tool that claims to catch a problem should be shown
catching one.

You need yeast-GEM.xml in the working directory:
    wget https://github.com/SysBioChalmers/yeast-GEM/raw/main/model/yeast-GEM.xml
"""
import itertools
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from preflight import (
    load_config, CobraEvaluator, run_diagnostics, build_proposal_sheet,
)

CONFIG = os.path.join(os.path.dirname(__file__), "..", "configs",
                       "yeast9_methyl_myristate.yaml")
TARGET = {"C14:0": 0.6, "C16:0": 0.2, "C18:0": 0.1, "C18:1": 0.1}


def main():
    config = load_config(CONFIG)
    evaluator = CobraEvaluator(config)

    report = run_diagnostics(evaluator, TARGET, max_order=2, cap=20)
    print(report.to_text())

    candidates = [
        c for order in (1, 2)
        for c in itertools.combinations(config.candidate_genes, order)
        if evaluator.is_valid(frozenset(c), TARGET)
    ][:20]

    sheet = build_proposal_sheet(evaluator, TARGET, candidates, top_n=5,
                                  diagnostics_report=report)
    print()
    print(sheet.to_text())

    out = os.path.join(os.path.dirname(__file__), "proposals_yeast9.csv")
    sheet.to_csv(out)
    print(f"\nWrote {out}")
    if sheet.gated:
        print("\nNo proposals were issued, by design. Fix the blocking checks "
              "above and re-run.")


if __name__ == "__main__":
    main()
