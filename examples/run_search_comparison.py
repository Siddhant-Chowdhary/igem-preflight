"""
run_search_comparison.py - which search should you actually use?

Replaces the RL-based example we started with. Both that example and the
score it ran on turned out to be unsound: the score cannot separate
candidates (run_diagnostics.py) and the RL agent loses to simpler strategies.

Measured over 20 random seeds on the synthetic landscape in benchmark.py,
up to 3 knockouts:

  genes  strategy      finds exact optimum   mean shortfall   evaluations
  10     brute force   100% (by definition)  0%                  175
  10     greedy         65%                  1.14%                23
  10     RL             10%                  10.13%            1,204

  20     brute force   100%                  0%                1,350
  20     greedy         40%                  1.21%                44

  30     brute force   100%                  0%                4,525
  30     greedy         50%                  0.63%                71

On our real evaluator (about 0.30 s per FBA solve), 10 genes at order 3:

  brute force    175 evals =  52 s   exact
  greedy          23 evals =   7 s   misses the best answer ~35% of the time
  RL           1,204 evals = 361 s   misses it ~90% of the time

So: use brute force at this scale. Saving 45 seconds is not worth a
one-in-three chance of handing the wet lab the wrong strain. Greedy earns its
place once enumeration stops being affordable, which is tens of genes at
order 4 and up, and that is not where we are.

Use the RL agent only if you can show on your own landscape that it beats
both. On ours it beat neither.
"""
import itertools
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from preflight import (
    load_config, CobraEvaluator, brute_force_best, greedy_search,
    run_diagnostics, combination_count,
)

CONFIG = os.path.join(os.path.dirname(__file__), "..", "configs",
                       "yeast9_methyl_myristate.yaml")
TARGET = {"C14:0": 0.6, "C16:0": 0.2, "C18:0": 0.1, "C18:1": 0.1}


def main():
    config = load_config(CONFIG)
    n = len(config.candidate_genes)
    total = combination_count(n, config.max_knockout_order)
    print(f"{n} candidate genes, order <= {config.max_knockout_order} "
          f"-> {total} combinations")
    print(f"Exhaustive enumeration is therefore ~{total * 0.3:.0f} s of FBA. "
          "Cheap, and exact.\n")

    evaluator = CobraEvaluator(config)

    # Check the score is worth searching over BEFORE searching over it.
    report = run_diagnostics(evaluator, TARGET, max_order=2, cap=20)
    print(f"Scoring setup usable: {report.usable}")
    if not report.usable:
        print("\nThe score cannot discriminate candidates, so neither strategy")
        print("below is meaningful on it. Running them anyway only to show the")
        print("mechanics; do not act on the output.\n")
        for r in report.blocking_failures:
            print(f"  {r.check_id} {r.name}: {r.headline}")
        print()

    print("Brute force (exact, recommended at this scale):")
    bf_combo, bf_score = brute_force_best(config, evaluator, TARGET)
    print(f"  {sorted(bf_combo)}  score={bf_score:.6f}")

    print("\nGreedy (approximate, cheaper -- for when enumeration is too big):")
    g_combo, g_score = greedy_search(config, evaluator, TARGET)
    print(f"  {sorted(g_combo)}  score={g_score:.6f}")

    gap = bf_score - g_score
    print(f"\nGreedy shortfall vs exact: {gap:.3e}")
    print("\nSee the module docstring for why brute force is the "
          "recommendation here.")


if __name__ == "__main__":
    main()
