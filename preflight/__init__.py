"""
preflight - checks whether a metabolic model's scoring function is measuring
anything, before anyone acts on what it says

We built this for an iGEM 2026 dry-lab project engineering yeast toward a
C14:0-rich fatty-acid profile, then generalised it so any team with a
COBRApy-compatible model can use it by writing a config file.

The short version of why: our scoring function ranked gene knockouts, brute
force found a clear winner out of 175 combinations, and the winner beat doing
nothing by 0.010%. No GC-MS could measure that. The search was fine. The
score was not, and nothing about the output said so.

Quick start:

    from preflight import load_config, CobraEvaluator, run_diagnostics

    config = load_config("configs/yeast9_methyl_myristate.yaml")
    evaluator = CobraEvaluator(config)
    target = {"C14:0": 0.6, "C16:0": 0.2, "C18:0": 0.1, "C18:1": 0.1}

    report = run_diagnostics(evaluator, target)
    print(report.to_text())
    print(report.usable)      # False means any ranking is an artefact

The checks only need an object with evaluate(combo, target) -> float. They do
not need COBRApy, our config system or our evaluator. Checks that need more
than that skip instead of crashing, so this works on any scoring function.

See README.md for results and INSTALL.md for setup.
"""

from .config import SearchConfig, load_config
from .filters import Rule, load_rules, is_valid_combo
from .agent import LinearQAgent, build_features, feature_dim
from .search import (
    run_episode, train_agent, brute_force_best, greedy_search,
    evaluate_agent_vs_bruteforce, sample_random_target_profile,
)
from .cobra_evaluator import CobraEvaluator, cosine_similarity, growth_penalty
from .diagnostics import (
    run_diagnostics, DiagnosticReport, CheckResult, SOLVER_NOISE_FLOOR,
    DETECTABLE_RELATIVE_MARGIN, MEASURABLE_RELATIVE_MARGIN,
)
from .proposals import build_proposal_sheet, ProposalSheet, Proposal
from .benchmark import (
    SyntheticLandscape, run_scaling_benchmark, benchmark_table,
    crossover_summary, combination_count,
)

__all__ = [
    "SearchConfig", "load_config",
    "Rule", "load_rules", "is_valid_combo",
    "LinearQAgent", "build_features", "feature_dim",
    "run_episode", "train_agent", "brute_force_best", "greedy_search",
    "evaluate_agent_vs_bruteforce", "sample_random_target_profile",
    "CobraEvaluator", "cosine_similarity", "growth_penalty",
    "run_diagnostics", "DiagnosticReport", "CheckResult",
    "SOLVER_NOISE_FLOOR", "DETECTABLE_RELATIVE_MARGIN",
    "MEASURABLE_RELATIVE_MARGIN",
    "build_proposal_sheet", "ProposalSheet", "Proposal",
    "SyntheticLandscape", "run_scaling_benchmark", "benchmark_table",
    "crossover_summary", "combination_count",
]

__version__ = "1.0.0"
