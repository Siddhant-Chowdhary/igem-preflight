"""
run_benchmark.py - was the reinforcement learning worth it?

Reproduces our search-strategy comparison. Runs entirely on a synthetic
landscape, so it needs no model file and finishes in a few minutes.

    python examples/run_benchmark.py

Our answer was no. Greedy hill-climbing reached 99.4-99.9% of the exhaustive
optimum at every epistasis level we tried, while the linear Q-learning agent
reached 80-85% using far more evaluations. Random search on the same budget
matched greedy, which told us the landscape has many near-optimal answers.

If you get a different answer on your own landscape, that is the interesting
case and we would like to hear about it.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

from preflight import (
    run_scaling_benchmark, benchmark_table, crossover_summary,
)
from preflight.benchmark import (
    SyntheticLandscape, search_exhaustive, search_greedy, search_random,
    search_rl,
)

# Measure this on your own evaluator and change it before quoting any
# projected time. 0.30 s is this project's 3-stage lexicographic Yeast9 cost.
FBA_SECONDS_PER_EVAL = 0.30


def scaling():
    print("=" * 74)
    print("SCALING: how cost and quality move with the candidate gene count")
    print("=" * 74)
    rows = run_scaling_benchmark(gene_counts=(10, 20, 30, 40, 50),
                                 max_order=3, n_episodes=300, seed=0)
    print()
    print(benchmark_table(rows, FBA_SECONDS_PER_EVAL))
    print()
    print(crossover_summary(rows, FBA_SECONDS_PER_EVAL))
    return rows


def epistasis_sweep(n_genes=30, seeds=5):
    print()
    print("=" * 74)
    print("EPISTASIS: does gene interaction break greedy hill-climbing?")
    print("=" * 74)
    print(f"{'epistasis':>10} {'greedy':>16} {'rl':>16} {'random':>16}"
          "   (% of exhaustive optimum, mean +- sd)")
    for eps in (0.0, 0.5, 1.0, 2.0, 4.0):
        g, r, q = [], [], []
        for seed in range(seeds):
            L = SyntheticLandscape(n_genes=n_genes, seed=seed,
                                   epistasis_strength=eps,
                                   epistasis_pairs=200)
            t = L.random_target(np.random.default_rng(seed + 100))
            L.reset_counter(); _, ex = search_exhaustive(L, t, 3)
            L.reset_counter(); _, gs = search_greedy(L, t, 3)
            L.reset_counter()
            _, rs = search_rl(L, L.random_target, t, 3, n_episodes=300,
                              seed=seed)
            budget = L.n_evals
            L.reset_counter()
            _, qs = search_random(L, t, 3, budget=budget, seed=seed)
            g.append(gs / ex); r.append(rs / ex); q.append(qs / ex)
        print(f"{eps:>10.1f} "
              f"{np.mean(g)*100:>11.1f}+-{np.std(g)*100:>3.0f}% "
              f"{np.mean(r)*100:>11.1f}+-{np.std(r)*100:>3.0f}% "
              f"{np.mean(q)*100:>11.1f}+-{np.std(q)*100:>3.0f}%", flush=True)
    print("\nRandom search matching the RL evaluation budget scoring as well")
    print("as greedy is itself a finding: this landscape has many")
    print("near-optimal solutions, and where that is true, search")
    print("sophistication buys very little.")


if __name__ == "__main__":
    scaling()
    epistasis_sweep()
