"""
Tests for the synthetic landscape used in the search benchmark.

Mostly sanity: the landscape has to be deterministic for a given seed,
different across seeds, bounded in [0,1], and honest about how many
evaluations a search used. If the evaluation counter is wrong then every
cost comparison in the benchmark is wrong too, so it gets its own test.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

from preflight.benchmark import (
    SyntheticLandscape, search_exhaustive, search_greedy, combination_count,
)


def test_landscape_is_deterministic():
    a = SyntheticLandscape(n_genes=8, seed=3)
    b = SyntheticLandscape(n_genes=8, seed=3)
    t = a.random_target(np.random.default_rng(0))
    assert a.evaluate(("G00", "G01"), t) == b.evaluate(("G00", "G01"), t)


def test_different_seeds_differ():
    a = SyntheticLandscape(n_genes=8, seed=1)
    b = SyntheticLandscape(n_genes=8, seed=2)
    t = a.random_target(np.random.default_rng(0))
    assert a.evaluate(("G00",), t) != b.evaluate(("G00",), t)


def test_scores_in_unit_range():
    L = SyntheticLandscape(n_genes=10, seed=0)
    t = L.random_target(np.random.default_rng(0))
    for combo in [(), ("G00",), ("G00", "G01"), ("G00", "G01", "G02")]:
        s = L.evaluate(combo, t)
        assert 0.0 <= s <= 1.0, f"{combo} -> {s}"


def test_eval_counter_counts_unique_only():
    L = SyntheticLandscape(n_genes=6, seed=0)
    t = L.random_target(np.random.default_rng(0))
    L.reset_counter()
    L.evaluate(("G00",), t)
    L.evaluate(("G00",), t)
    assert L.n_evals == 1


def test_combination_count_matches_enumeration():
    L = SyntheticLandscape(n_genes=7, seed=0)
    t = L.random_target(np.random.default_rng(0))
    L.reset_counter()
    search_exhaustive(L, t, 3)
    assert L.n_evals == combination_count(7, 3)


def test_greedy_never_beats_exhaustive():
    for seed in range(3):
        L = SyntheticLandscape(n_genes=12, seed=seed)
        t = L.random_target(np.random.default_rng(seed))
        _, ex = search_exhaustive(L, t, 3)
        _, gr = search_greedy(L, t, 3)
        assert gr <= ex + 1e-12


def test_greedy_is_cheap():
    L = SyntheticLandscape(n_genes=20, seed=0)
    t = L.random_target(np.random.default_rng(0))
    L.reset_counter()
    search_greedy(L, t, 3)
    assert L.n_evals < combination_count(20, 3) / 10


if __name__ == "__main__":
    tests = [v for k, v in list(globals().items()) if k.startswith("test_")]
    p = f = 0
    for t in tests:
        try:
            t(); p += 1
        except Exception as e:  # noqa: BLE001
            f += 1; print(f"FAILED {t.__name__}: {e}")
    print(f"\n{p} passed, {f} failed")
