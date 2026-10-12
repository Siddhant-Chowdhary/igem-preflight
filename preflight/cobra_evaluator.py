"""
cobra_evaluator.py - scores a knockout combination against a target profile
using FBA, for any COBRApy-compatible model

Nothing organism-specific lives here. Which model, which genes, which
reactions, how the growth penalty is shaped - all of it comes from a
SearchConfig (see config.py). To point this at a different organism you write
a new YAML file and change nothing in this directory.

WHY THE OBJECTIVE HAS THREE STAGES

We started with the obvious thing: maximise a weighted sum of the target
exchange reactions, with the weights set to the ratio we wanted. It does not
work. A linear program only ever returns a vertex, so the solver gives you
the single best-paying product and ignores the rest no matter what weights
you set. We confirmed this on our own model at several weight ratios and got
the same all-or-nothing answer every time.

What works is a lexicographic solve:

  Stage 1  maximise total flux across all target reactions, call it T_max
  Stage 2  require at least alpha * T_max total, then minimise the L1
           deviation from the target RATIO. This is the stage that actually
           chases a blend instead of one product.
  Stage 3  lock that ratio match in place, then maximise growth.

Stage 3 exists because we got it wrong the first time. We removed the growth
constraint and let the LP decide, which gave it no reason to keep any growth
at all, so every combination collapsed to the minimum and the growth penalty
below became a constant that could not distinguish anything. Locking the
ratio match first and then asking how much growth is still affordable gives a
number that actually varies between combinations.

THE GROWTH PENALTY

A hard pass/fail growth cutoff gives a search no gradient near the boundary,
so we used a smooth Hill-type curve instead. The shape is borrowed from the
toxicity terms in the dynamic metabolic-control literature, adapted for
steady-state FBA:

    g = achieved_growth / wt_growth              clipped to [0, 1]
    deficit d = 1 - g
    penalty(d) = d^n / (K^n + d^n)
    score = cosine_similarity(predicted, target) * (1 - penalty)

There is still a small hard floor (config: hard_min_growth_fraction) to throw
out genuinely dead strains. A smooth curve should not paper over a genotype
that cannot grow at all.
"""

import warnings
warnings.filterwarnings("ignore")
import numpy as np
from cobra.io import read_sbml_model

from .config import SearchConfig
from .filters import is_valid_combo


def cosine_similarity(predicted: dict, target: dict) -> float:
    fas = list(target.keys())
    p = np.array([predicted.get(k, 0.0) for k in fas])
    t = np.array([target.get(k, 0.0) for k in fas])
    if np.linalg.norm(p) < 1e-12 or np.linalg.norm(t) < 1e-12:
        return 0.0
    return float(np.dot(p, t) / (np.linalg.norm(p) * np.linalg.norm(t)))


def growth_penalty(achieved_growth: float, wt_growth: float, K: float, n: float) -> float:
    """Smooth Hill-type multiplier in [0, 1]. 1.0 = no growth lost."""
    if wt_growth <= 0:
        return 0.0
    g = max(0.0, min(1.0, achieved_growth / wt_growth))
    d = 1.0 - g
    penalty = (d ** n) / (K ** n + d ** n)
    return 1.0 - penalty


class CobraEvaluator:
    """
    Loads the model once (expensive) and exposes an `evaluate()` method
    matching the interface the search engine expects:
        evaluate(knockout_combo: tuple[str, ...], target_profile: dict) -> float
    Results are cached by (combo, target_profile) since FBA solves are the
    expensive step and the same combo/target pair is often re-queried
    during RL training and validation.
    """

    def __init__(self, config: SearchConfig):
        self.config = config
        print(f"Loading model for '{config.name}' ({config.model_path})...")
        self._model = read_sbml_model(config.model_path)
        wt_sol = self._model.optimize()
        self._wt_growth = wt_sol.objective_value
        print(f"Loaded. WT growth rate: {self._wt_growth:.6f}")
        self._cache: dict = {}
        self._parts: dict = {}   # key -> (raw cosine, achieved growth)

    def is_valid(self, combo: frozenset, target_profile: dict, open_gates: set | None = None) -> bool:
        return is_valid_combo(combo, target_profile, self.config.filter_rules, open_gates)

    def diagnose(self, knockout_combo, target_profile: dict) -> dict:
        """Break a score into its parts: raw cosine match, growth as a fraction
        of wild-type, growth multiplier, and final score. For debugging/reporting."""
        score = self.evaluate(knockout_combo, target_profile)
        key = (tuple(sorted(knockout_combo)), tuple(sorted(target_profile.items())))
        raw, g = self._parts.get(key, (0.0, 0.0))
        return {"cosine": raw, "growth_fraction": g / self._wt_growth if self._wt_growth else 0.0,
                "score": score}

    def evaluate(self, knockout_combo, target_profile: dict) -> float:
        key = (tuple(sorted(knockout_combo)), tuple(sorted(target_profile.items())))
        if key in self._cache:
            return self._cache[key]

        cfg = self.config
        with self._model as model:
            for gene_common in knockout_combo:
                model_id = cfg.gene_to_model_id.get(gene_common)
                if model_id is None:
                    continue
                try:
                    model.genes.get_by_id(model_id).knock_out()
                except KeyError:
                    continue

            biomass_rxn = model.reactions.get_by_id(cfg.biomass_reaction_id)
            biomass_rxn.lower_bound = cfg.hard_min_growth_fraction * self._wt_growth

            target_rxns = {
                name: model.reactions.get_by_id(rid)
                for name, rid in cfg.exchange_reactions.items() if name in target_profile
            }
            if not target_rxns:
                self._cache[key] = 0.0
                return 0.0

            # Stage 1: max achievable total output across target reactions
            model.objective = {rxn: 1.0 for rxn in target_rxns.values()}
            model.objective.direction = "max"
            T_max = model.slim_optimize()
            if T_max is None or T_max < 1e-9:
                self._cache[key] = 0.0
                return 0.0

            # Stage 2: floor total output, minimize L1 deviation from target ratio
            total_con = model.problem.Constraint(
                sum(rxn.flux_expression for rxn in target_rxns.values()),
                lb=cfg.ratio_match_alpha * T_max, name="min_total_output")
            model.add_cons_vars([total_con])

            dev_vars = [model.problem.Variable(f"dev_{name}", lb=0) for name in target_rxns]
            model.add_cons_vars(dev_vars)
            for name, rxn in target_rxns.items():
                d = model.variables[f"dev_{name}"]
                target_flux = target_profile[name] * T_max
                model.add_cons_vars([model.problem.Constraint(rxn.flux_expression - d, ub=target_flux)])
                model.add_cons_vars([model.problem.Constraint(rxn.flux_expression + d, lb=target_flux)])

            model.objective = sum(model.variables[f"dev_{name}"] for name in target_rxns)
            model.objective.direction = "min"
            min_dev = model.slim_optimize()
            if min_dev is None:
                self._cache[key] = 0.0
                return 0.0

            # Stage 3: lock the ratio match, then maximize growth
            dev_lock = model.problem.Constraint(
                sum(model.variables[f"dev_{name}"] for name in target_rxns),
                ub=min_dev, name="lock_ratio_match")
            model.add_cons_vars([dev_lock])
            model.objective = biomass_rxn
            model.objective.direction = "max"
            sol = model.optimize()

            if sol.status != "optimal":
                self._cache[key] = 0.0
                return 0.0

            predicted_profile = {
                name: max(0.0, sol.fluxes.get(rid, 0.0))
                for name, rid in cfg.exchange_reactions.items()
            }
            achieved_growth = sol.fluxes.get(cfg.biomass_reaction_id, 0.0)

        raw_score = cosine_similarity(predicted_profile, target_profile)
        penalty_mult = growth_penalty(achieved_growth, self._wt_growth,
                                       cfg.growth_penalty_K, cfg.growth_penalty_n)
        score = raw_score * penalty_mult
        self._parts[key] = (raw_score, achieved_growth)
        self._cache[key] = score
        return score
