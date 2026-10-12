"""
diagnostics.py - six checks on a scoring function

We wrote this after our own pipeline handed us a ranked list of knockouts
that turned out to mean nothing. The search had worked. The arithmetic was
right. The best combination beat doing nothing by 0.010%, which no GC-MS we
have access to could ever measure.

There are good test suites for genome-scale models already. memote (Lieven
et al. 2020, doi:10.1038/s41587-020-0446-y) checks annotations, mass balance,
stoichiometric consistency and energy-generating cycles. We ran it and we
recommend it. What none of them check is the scoring function you write on
top of the model, which is where our problem was.

So these are six checks on the scoring setup. Each one returns PASS, WARN or
FAIL with a sentence explaining what it found and what to do about it.

HOW MUCH OF THIS IS NEW

We looked before claiming anything. We went through memote, FBCModelTests,
FROG, COBRApy's own validators, FVA, flux sampling, loopless FBA, pFBA,
production envelopes, and the strain-design tools OptKnock, OptGene,
OptForce, RobustKnock, cameo, MEWpy, OptFlux and StrainDesign. Of the six
checks here, four already exist somewhere:

  D1  null-intervention parity   We could not find this automated anywhere.
                                 StrainDesign's docs mention in passing that
                                 the wild type can come back as a solution,
                                 but nothing tests for it.
  D2  landscape discrimination   Also could not find this anywhere. No tool
                                 we looked at reports whether its own score
                                 can separate candidates.
  D3  scalarization degeneracy   Not new. This is FVA (Mahadevan & Schilling
                                 2003) and alternate-optima analysis (Reed &
                                 Palsson 2004), and the underlying result
                                 about weighted sums only reaching vertices
                                 is Das & Dennis 1997
                                 (doi:10.1007/BF01197559). We automate it and
                                 explain it. We did not discover it.
  D4  free production            Not new, and memote already ships this as
                                 find_reactions_with_unbounded_flux_default_
                                 condition. Ours just looks at your specific
                                 target reactions so it lands in the same
                                 report. Run memote as well.
  D5  growth confound            Not new (Schuetz et al. 2007). Careful here:
                                 growth-coupled production is usually the
                                 point of a strain design (Klamt & Mahadevan
                                 2015), so this check never returns FAIL. It
                                 only asks whether the coupling is what you
                                 meant.
  D6  parameter sensitivity      Ordinary sensitivity analysis.

What we think is actually ours: D1, D2, and the fact that the proposal step
refuses to run when they fail (see proposals.py). Please do not describe the
whole suite as novel, because four of six are not.

WHY THE THRESHOLDS ARE WHAT THEY ARE

Our first version compared score differences against solver tolerance, about
one part in a million. It passed the broken pipeline described above, because
2e-5 is bigger than 1e-6. The check was useless.

Solver tolerance is the wrong yardstick. The question is not whether a
difference is real arithmetic. It is whether anyone could ever measure it.
GC-MS fatty-acid composition varies by roughly 1-5% between biological
replicates, so a model that ranks two designs 0.01% apart is making a
distinction no experiment can test. The thresholds below are therefore
relative and tied to measurement resolution, not to the solver. If your assay
is tighter than ours, tighten them.
"""

from dataclasses import dataclass, field
import itertools
from typing import Any

import numpy as np

# Below this, two scores are indistinguishable from solver noise rather than
# genuinely different. LP solvers typically converge to ~1e-9 but FBA scores
# built from several solves accumulate more error than that.
SOLVER_NOISE_FLOOR = 1e-6

# THE THRESHOLD THAT ACTUALLY MATTERS, AND WHY IT IS NOT THE SOLVER'S
#
# Our first version of this module compared score differences against
# SOLVER_NOISE_FLOOR, and it let our own broken pipeline through with a
# PASS-equivalent: the best knockout beat doing nothing by 2e-5, which is
# comfortably above 1e-6, so the check said "fine" about a landscape we
# already knew was flat.
#
# Solver tolerance is the wrong standard. The question a strain-design score
# has to answer is not "is this difference real arithmetic?" but "could the
# wet lab ever tell these two strains apart?" GC-MS fatty-acid composition,
# reported as a fraction of total peak area, carries roughly 1-5% relative
# error between biological replicates. A model that ranks two designs 0.01%
# apart is making a distinction no experiment can test, so for the purpose of
# proposing bench work the landscape is flat even though the arithmetic is
# sound.
#
# These are therefore RELATIVE thresholds tied to measurement resolution, and
# they are deliberately conservative. Override them for your own assay:
# a method with tighter replicate error earns a smaller threshold.
MEASURABLE_RELATIVE_MARGIN = 0.05   # 5% -- a difference worth building a strain for
DETECTABLE_RELATIVE_MARGIN = 0.01   # 1% -- the optimistic edge of GC-MS resolution

PASS, WARN, FAIL, SKIP = "PASS", "WARN", "FAIL", "SKIP"


@dataclass
class CheckResult:
    check_id: str
    name: str
    status: str
    headline: str           # one plain-language sentence
    detail: str = ""        # what was measured
    fix: str = ""           # what to do about it
    evidence: dict = field(default_factory=dict)
    citations: list = field(default_factory=list)
    novelty: str = "known"  # "new" or "known" -- see module docstring

    @property
    def blocking(self) -> bool:
        """A FAIL on a blocking check means downstream rankings are unusable."""
        return self.status == FAIL and self.check_id in ("D1", "D2")


@dataclass
class DiagnosticReport:
    results: list
    target_profile: dict
    n_candidates_tested: int

    @property
    def blocking_failures(self) -> list:
        return [r for r in self.results if r.blocking]

    @property
    def usable(self) -> bool:
        """True when the scoring setup can support a ranking of candidates."""
        return not self.blocking_failures

    def counts(self) -> dict:
        out = {PASS: 0, WARN: 0, FAIL: 0, SKIP: 0}
        for r in self.results:
            out[r.status] += 1
        return out

    def to_text(self) -> str:
        lines = []
        w = 74
        lines.append("=" * w)
        lines.append("SCORING-SETUP DIAGNOSTIC REPORT")
        lines.append("=" * w)
        tp = ", ".join(f"{k}={v}" for k, v in self.target_profile.items())
        lines.append(f"Target profile    : {tp}")
        lines.append(f"Candidates tested : {self.n_candidates_tested}")
        c = self.counts()
        lines.append(f"Result            : {c[PASS]} pass, {c[WARN]} warn, "
                     f"{c[FAIL]} fail, {c[SKIP]} skipped")
        lines.append("")
        verdict = ("USABLE -- rankings from this scoring setup are meaningful"
                   if self.usable else
                   "NOT USABLE -- a blocking check failed; rankings are "
                   "meaningless")
        lines.append(f"VERDICT: {verdict}")
        lines.append("=" * w)
        for r in self.results:
            tag = "[new]" if r.novelty == "new" else ""
            lines.append("")
            lines.append(f"{r.status:5} {r.check_id}  {r.name} {tag}".rstrip())
            lines.append(f"      {r.headline}")
            if r.detail:
                for ln in r.detail.split("\n"):
                    lines.append(f"      {ln}")
            if r.status in (WARN, FAIL) and r.fix:
                lines.append(f"      FIX: {r.fix}")
            if r.citations:
                lines.append(f"      refs: {'; '.join(r.citations)}")
        lines.append("")
        lines.append("=" * w)
        lines.append("Run memote (doi:10.1038/s41587-020-0446-y) as well -- it")
        lines.append("checks the MODEL; this suite checks the SCORING SETUP.")
        lines.append("=" * w)
        return "\n".join(lines)

    def to_markdown(self) -> str:
        icon = {PASS: "PASS", WARN: "WARN", FAIL: "FAIL", SKIP: "n/a"}
        lines = ["# Scoring-setup diagnostic report", ""]
        verdict = ("**Usable.** Rankings from this scoring setup are meaningful."
                   if self.usable else
                   "**Not usable.** A blocking check failed, so any ranking "
                   "this setup produces is meaningless.")
        lines += [verdict, "",
                  f"Target profile: `{self.target_profile}`  ",
                  f"Candidates tested: {self.n_candidates_tested}", "",
                  "| Check | Result | Finding |", "| --- | --- | --- |"]
        for r in self.results:
            head = r.headline.replace("|", "/")
            lines.append(f"| {r.check_id} {r.name} | {icon[r.status]} | {head} |")
        lines += ["", "## Detail", ""]
        for r in self.results:
            lines.append(f"### {r.check_id} {r.name} — {icon[r.status]}")
            lines.append(r.headline)
            if r.detail:
                lines += ["", "```", r.detail, "```"]
            if r.status in (WARN, FAIL) and r.fix:
                lines += ["", f"**Fix:** {r.fix}"]
            if r.citations:
                lines += ["", f"*References: {'; '.join(r.citations)}*"]
            lines.append("")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _default_candidates(evaluator, target_profile, max_order=2, cap=40):
    """Candidate combos to probe, taken from the evaluator's own config and
    filters so we test exactly what the search would actually consider."""
    cfg = getattr(evaluator, "config", None)
    if cfg is None:
        raise ValueError(
            "evaluator has no .config; pass candidate_combos= explicitly")
    genes = cfg.candidate_genes
    order_cap = min(max_order, cfg.max_knockout_order)
    out = []
    for order in range(1, order_cap + 1):
        for combo in itertools.combinations(genes, order):
            if hasattr(evaluator, "is_valid") and not evaluator.is_valid(
                    frozenset(combo), target_profile):
                continue
            out.append(combo)
    if len(out) > cap:
        step = len(out) / cap
        out = [out[int(i * step)] for i in range(cap)]
    return out


def _score_all(evaluator, combos, target_profile):
    return {c: float(evaluator.evaluate(c, target_profile)) for c in combos}


# ---------------------------------------------------------------------------
# D1 -- null-intervention parity                                       [NEW]
# ---------------------------------------------------------------------------

def check_null_intervention_parity(evaluator, target_profile, scores) -> CheckResult:
    """Does doing nothing score as well as the best intervention?

    If the unedited strain scores as highly as the best knockout, the score
    is not measuring the effect of the intervention. Every ranking built on
    it is then an artefact. This is the single most informative check in the
    suite and the cheapest to run.
    """
    null_score = float(evaluator.evaluate((), target_profile))
    if not scores:
        return CheckResult("D1", "null-intervention parity", SKIP,
                           "No candidates to compare against.", novelty="new")
    best_combo = max(scores, key=scores.get)
    best = scores[best_combo]
    margin = best - null_score
    rel = margin / abs(null_score) if abs(null_score) > 1e-12 else float("inf")

    ev = {"null_score": null_score, "best_score": best,
          "best_combo": list(best_combo), "margin": margin,
          "relative_margin": rel}
    detail = (f"no-intervention score : {null_score:.6f}\n"
              f"best candidate        : {list(best_combo)} -> {best:.6f}\n"
              f"margin                : {margin:.3e}  "
              f"({rel*100:.3f}% of baseline)")
    cites = ["StrainDesign docs note the wild type can be returned as a "
             "solution (doi:10.1093/bioinformatics/btac632)"]

    fix = ("The objective is essentially satisfiable without any edit. "
           "Constrain the route the product must take (e.g. require flux "
           "through the synthesis pathway), bound the product exchange "
           "realistically, or score product yield at fixed growth instead of "
           "profile shape.")

    if margin <= SOLVER_NOISE_FLOOR:
        return CheckResult(
            "D1", "null-intervention parity", FAIL,
            "Doing nothing scores as well as the best intervention, so the "
            "score is not measuring the intervention at all.",
            detail, fix, ev, cites, novelty="new")
    if rel < DETECTABLE_RELATIVE_MARGIN:
        return CheckResult(
            "D1", "null-intervention parity", FAIL,
            f"The best intervention beats doing nothing by only "
            f"{rel*100:.3f}%, which is below what the assay could resolve, so "
            "the ranking is untestable even though the arithmetic is sound.",
            detail, fix, ev, cites, novelty="new")
    if rel < MEASURABLE_RELATIVE_MARGIN:
        return CheckResult(
            "D1", "null-intervention parity", WARN,
            f"The best intervention beats doing nothing by {rel*100:.2f}%, "
            "near the edge of replicate measurement error.",
            detail, "Treat the ranking as provisional and check the margin "
            "against your assay's own replicate error before committing bench "
            "time.", ev, cites, novelty="new")
    return CheckResult(
        "D1", "null-intervention parity", PASS,
        f"The best intervention beats doing nothing by {rel*100:.1f}%, a "
        "margin an experiment could resolve.",
        detail, evidence=ev, citations=cites, novelty="new")


# ---------------------------------------------------------------------------
# D2 -- landscape discrimination                                       [NEW]
# ---------------------------------------------------------------------------

def check_landscape_discrimination(evaluator, target_profile, scores) -> CheckResult:
    """Can the score tell candidates apart at all?

    Reports the spread of scores across candidates and how many genuinely
    distinguishable tiers exist once solver noise is accounted for. A search
    over a landscape with one tier is a random number generator with extra
    steps.
    """
    if len(scores) < 2:
        return CheckResult("D2", "landscape discrimination", SKIP,
                           "Need at least two candidates.", novelty="new")
    vals = np.array(sorted(scores.values()))
    spread = float(vals[-1] - vals[0])
    std = float(np.std(vals))
    scale = float(np.median(np.abs(vals)))
    rel_spread = spread / scale if scale > 1e-12 else float("inf")

    # Two tier counts, because they answer different questions: what the
    # arithmetic can separate, and what an experiment could separate.
    def count_tiers(min_gap):
        tiers, last = 1, vals[0]
        for v in vals[1:]:
            if v - last > min_gap:
                tiers += 1
                last = v
        return tiers

    tiers_numeric = count_tiers(SOLVER_NOISE_FLOOR)
    tiers_measurable = count_tiers(DETECTABLE_RELATIVE_MARGIN * scale)

    ev = {"n": len(vals), "min": float(vals[0]), "max": float(vals[-1]),
          "spread": spread, "relative_spread": rel_spread, "std": std,
          "tiers_numeric": tiers_numeric,
          "tiers_measurable": tiers_measurable}
    detail = (f"candidates              : {len(vals)}\n"
              f"score range             : {vals[0]:.6f} .. {vals[-1]:.6f}\n"
              f"spread                  : {spread:.3e} "
              f"({rel_spread*100:.3f}% of typical score)\n"
              f"std deviation           : {std:.3e}\n"
              f"tiers the solver can see: {tiers_numeric} of {len(vals)}\n"
              f"tiers an assay can see  : {tiers_measurable} of {len(vals)}")
    cites = ["Bernstein et al. 2021 (doi:10.1186/s13059-021-02289-z)",
             "Kaste & Shachar-Hill 2024 (arXiv:2303.12651)"]

    if spread <= SOLVER_NOISE_FLOOR:
        return CheckResult(
            "D2", "landscape discrimination", FAIL,
            "Every candidate scores identically, so the score cannot rank "
            "anything.",
            detail, "The score is flat. Find what it is actually responding "
            "to (run D5), then redefine it so the intervention changes it.",
            ev, cites, novelty="new")
    if tiers_measurable <= 1:
        return CheckResult(
            "D2", "landscape discrimination", FAIL,
            f"The solver separates these candidates into {tiers_numeric} "
            "levels, but all of them fall inside one measurement-resolution "
            "band, so no experiment could confirm the ordering.",
            detail, "The ranking is numerically real and experimentally "
            "meaningless. Redefine the score so differences exceed assay "
            "error, or report a single undifferentiated set.",
            ev, cites, novelty="new")
    if tiers_measurable <= 2:
        return CheckResult(
            "D2", "landscape discrimination", WARN,
            f"Only {tiers_measurable} measurement-resolvable levels across "
            f"{len(vals)} candidates.",
            detail, "A ranking this coarse can order tiers but not candidates "
            "within a tier. Report tiers, not a ranked list.",
            ev, cites, novelty="new")
    return CheckResult(
        "D2", "landscape discrimination", PASS,
        f"The score separates candidates into {tiers_measurable} levels an "
        "experiment could resolve.",
        detail, evidence=ev, citations=cites, novelty="new")


# ---------------------------------------------------------------------------
# D3 -- scalarization degeneracy                                     [KNOWN]
# ---------------------------------------------------------------------------

def check_scalarization_degeneracy(evaluator, target_profile,
                                    n_weight_samples=5) -> CheckResult:
    """Does a weighted-sum objective over the target products collapse?

    NOTE ON WHAT THIS TESTS: this check builds a weighted-sum objective on
    the model directly and sweeps the weights. It is therefore a test of the
    weighted-sum FORMULATION, not of whatever objective your evaluator
    currently uses. A FAIL means "if you scalarize these products with
    weights, the weights will not control the blend" -- which is worth
    knowing even if you have already moved to goal programming, because it
    documents why you had to.

    This is FVA / alternate-optima analysis applied to a multi-product
    objective: a linear program only ever returns a vertex, so maximizing a
    weighted sum of product exchanges gives you the single best-paying
    product at almost every weight ratio, not a blend. Known result, not our
    discovery -- see the citations.
    """
    model = getattr(evaluator, "_model", None)
    cfg = getattr(evaluator, "config", None)
    cites = ["Das & Dennis 1997 (doi:10.1007/BF01197559)",
             "Mahadevan & Schilling 2003 (doi:10.1016/j.ymben.2003.09.002)",
             "Reed & Palsson 2004 (doi:10.1101/gr.2546004)"]
    if model is None or cfg is None:
        return CheckResult("D3", "scalarization degeneracy", SKIP,
                           "Needs a COBRA model on the evaluator.",
                           citations=cites)

    products = [p for p in cfg.exchange_reactions if p in target_profile]
    if len(products) < 2:
        return CheckResult("D3", "scalarization degeneracy", SKIP,
                           "Only one target product; nothing to scalarize.",
                           citations=cites)

    rng = np.random.default_rng(0)
    supports = []
    rows = []
    try:
        for i in range(n_weight_samples):
            w = (np.ones(len(products)) / len(products) if i == 0
                 else rng.dirichlet(np.ones(len(products))))
            with model as m:
                rxns = {p: m.reactions.get_by_id(cfg.exchange_reactions[p])
                        for p in products}
                m.objective = {rxns[p]: float(w[j])
                               for j, p in enumerate(products)}
                m.objective.direction = "max"
                sol = m.optimize()
                if sol.status != "optimal":
                    continue
                fl = {p: max(0.0, sol.fluxes.get(
                    cfg.exchange_reactions[p], 0.0)) for p in products}
            total = sum(fl.values())
            active = frozenset(p for p, v in fl.items()
                               if total > 0 and v / total > 0.01)
            supports.append(active)
            rows.append((np.round(w, 2).tolist(), len(active),
                         sorted(active)))
    except Exception as exc:  # noqa: BLE001 -- diagnostics must not crash a run
        return CheckResult("D3", "scalarization degeneracy", SKIP,
                           f"Could not run the weight sweep: {exc}",
                           citations=cites)

    if not supports:
        return CheckResult("D3", "scalarization degeneracy", SKIP,
                           "No optimal solutions in the weight sweep.",
                           citations=cites)

    sizes = [len(s) for s in supports]
    detail = "\n".join(
        f"weights {r[0]} -> {r[1]} product(s) active: {r[2]}" for r in rows)
    ev = {"n_samples": len(supports), "support_sizes": sizes,
          "products": products}

    if max(sizes) <= 1:
        return CheckResult(
            "D3", "scalarization degeneracy", FAIL,
            "A weighted sum over these products cannot control the blend: at "
            "every weight tested the solver produced just one product. (This "
            "tests the weighted-sum formulation itself, not your current "
            "objective.)",
            detail, "Do not use a weighted sum for a ratio target. Use "
            "lexicographic / goal programming: maximize total output first, "
            "then minimize deviation from the target ratio. Or inspect the "
            "trade-off surface with COBRApy's production_envelope.",
            ev, cites)
    if len(set(supports)) == 1 and len(supports) > 1:
        return CheckResult(
            "D3", "scalarization degeneracy", WARN,
            "The same set of products stays active at every weight, so the "
            "weights have little influence.",
            detail, "Confirm the weights are doing what you expect; consider "
            "a goal-programming formulation.", ev, cites)
    return CheckResult(
        "D3", "scalarization degeneracy", PASS,
        "Different weights produce different product blends.",
        detail, evidence=ev, citations=cites)


# ---------------------------------------------------------------------------
# D4 -- free production / unbounded export                           [KNOWN]
# ---------------------------------------------------------------------------

def check_free_production(evaluator, target_profile) -> CheckResult:
    """Can the target be produced without consuming anything?

    memote already ships this class of check
    (find_reactions_with_unbounded_flux_default_condition,
    detect_energy_generating_cycles) and you should run memote too. This is a
    convenience check aimed only at the specific target reactions in your
    config, so it appears in the same report as the rest.
    """
    model = getattr(evaluator, "_model", None)
    cfg = getattr(evaluator, "config", None)
    cites = ["memote (doi:10.1038/s41587-020-0446-y) -- run it as well",
             "Fritzemeier et al. 2017 (doi:10.1371/journal.pcbi.1005494)",
             "Thiele & Palsson 2010 (doi:10.1038/nprot.2009.203)"]
    if model is None or cfg is None:
        return CheckResult("D4", "free production", SKIP,
                           "Needs a COBRA model on the evaluator.",
                           citations=cites)

    findings, rows = [], []
    try:
        for name, rid in cfg.exchange_reactions.items():
            if name not in target_profile:
                continue
            rxn = model.reactions.get_by_id(rid)
            ub, lb = rxn.upper_bound, rxn.lower_bound
            huge = abs(ub) >= 1000
            rows.append(f"{name:8} {rid:10} bounds [{lb:g}, {ub:g}]"
                        f"{'   <- effectively unbounded' if huge else ''}")
            if huge:
                findings.append(f"{name} ({rid}) upper bound {ub:g}")
    except Exception as exc:  # noqa: BLE001
        return CheckResult("D4", "free production", SKIP,
                           f"Could not inspect bounds: {exc}", citations=cites)

    detail = "\n".join(rows)
    ev = {"effectively_unbounded": findings}
    if findings:
        return CheckResult(
            "D4", "free production", WARN,
            f"{len(findings)} target exchange reaction(s) are effectively "
            "unbounded, so the model may produce the target more freely than "
            "a cell could.",
            detail, "Set product exchange upper bounds from measured or "
            "literature secretion rates, and run memote's unbounded-flux and "
            "energy-cycle checks on the model.", ev, cites)
    return CheckResult(
        "D4", "free production", PASS,
        "Target exchange reactions carry finite, plausible bounds.",
        detail, evidence=ev, citations=cites)


# ---------------------------------------------------------------------------
# D5 -- growth confound                                              [KNOWN]
# ---------------------------------------------------------------------------

def check_growth_confound(evaluator, target_profile, scores) -> CheckResult:
    """Is the score tracking product formation, or only retained growth?

    CAREFUL: growth-coupled production is the INTENDED mechanism of most
    strain designs (Klamt & Mahadevan 2015; von Kamp & Klamt 2017), so
    correlation between score and growth is not by itself a defect. This
    check therefore never returns FAIL. What it flags is the specific case
    where the product term contributes no variance at all -- the score moves
    only with growth, so the product half of the objective is doing no work.
    """
    cites = ["Schuetz et al. 2007 (doi:10.1038/msb4100162)",
             "Feist & Palsson 2010 (doi:10.1016/j.mib.2010.03.003)",
             "Klamt & Mahadevan 2015 (doi:10.1016/j.ymben.2015.05.006) -- "
             "growth coupling is often intended"]
    if not hasattr(evaluator, "diagnose"):
        return CheckResult("D5", "growth confound", SKIP,
                           "Evaluator has no diagnose() method to decompose "
                           "the score.", citations=cites)
    if len(scores) < 3:
        return CheckResult("D5", "growth confound", SKIP,
                           "Need at least three candidates.", citations=cites)

    prod, growth, total = [], [], []
    try:
        for combo in scores:
            d = evaluator.diagnose(combo, target_profile)
            prod.append(d.get("cosine", 0.0))
            growth.append(d.get("growth_fraction", 0.0))
            total.append(d.get("score", scores[combo]))
    except Exception as exc:  # noqa: BLE001
        return CheckResult("D5", "growth confound", SKIP,
                           f"Could not decompose scores: {exc}",
                           citations=cites)

    prod, growth, total = np.array(prod), np.array(growth), np.array(total)
    pv, gv = float(np.std(prod)), float(np.std(growth))
    # Judge each term's variation against its own scale, so "the product term
    # technically moves by 1e-5" does not count as the term doing work.
    p_scale = max(float(np.median(np.abs(prod))), 1e-12)
    g_scale = max(float(np.median(np.abs(growth))), 1e-12)
    p_inert = (pv / p_scale) < DETECTABLE_RELATIVE_MARGIN
    g_inert = (gv / g_scale) < DETECTABLE_RELATIVE_MARGIN

    def corr(a, b):
        if np.std(a) <= 1e-12 or np.std(b) <= 1e-12:
            return float("nan")
        return float(np.corrcoef(a, b)[0, 1])

    cg, cp = corr(total, growth), corr(total, prod)
    ev = {"product_match_std": pv, "growth_std": gv,
          "corr_score_growth": cg, "corr_score_product": cp}
    detail = (f"product-match (cosine) spread : {pv:.3e} "
              f"({pv/p_scale*100:.3f}% of its scale"
              f"{', inert' if p_inert else ''})\n"
              f"growth-fraction spread        : {gv:.3e} "
              f"({gv/g_scale*100:.3f}% of its scale"
              f"{', inert' if g_inert else ''})\n"
              f"corr(score, growth)           : {cg:.3f}\n"
              f"corr(score, product match)    : {cp:.3f}")
    ev |= {"product_inert": bool(p_inert), "growth_inert": bool(g_inert)}

    if p_inert and not g_inert:
        return CheckResult(
            "D5", "growth confound", WARN,
            "The product-match term is identical for every candidate, so the "
            "score is driven entirely by retained growth.",
            detail, "Your score is a growth measurement wearing a "
            "product-match label. If growth coupling is what you intended, "
            "say so explicitly and score it directly; if not, constrain the "
            "product route so the match term varies.", ev, cites)
    if p_inert and g_inert:
        return CheckResult(
            "D5", "growth confound", WARN,
            "Neither the product-match nor the growth term varies "
            "meaningfully across candidates.",
            detail, "Nothing in the score responds to the intervention. See "
            "D1 and D2.", ev, cites)
    if not np.isnan(cg) and abs(cg) > 0.95 and (np.isnan(cp) or abs(cp) < 0.3):
        return CheckResult(
            "D5", "growth confound", WARN,
            "The score follows growth almost exactly and barely follows the "
            "product match.",
            detail, "Check this is the growth-coupled behaviour you intended "
            "rather than the product term being inert.", ev, cites)
    if g_inert and not p_inert:
        return CheckResult(
            "D5", "growth confound", PASS,
            "The score follows the product-match term; growth is constant "
            "across these candidates, so it is not confounding the ranking.",
            detail, evidence=ev, citations=cites)
    return CheckResult(
        "D5", "growth confound", PASS,
        "Both the product-match and growth terms vary and contribute to the "
        "score.", detail, evidence=ev, citations=cites)


# ---------------------------------------------------------------------------
# D6 -- parameter sensitivity                                        [KNOWN]
# ---------------------------------------------------------------------------

def check_parameter_sensitivity(evaluator, target_profile, scores,
                                 probe_combos=None) -> CheckResult:
    """Do the scoring parameters actually change anything?

    Ordinary sensitivity analysis. If scores are unchanged when a key
    parameter moves, that parameter's constraint is not binding, which often
    means it is not doing the job you think it is.
    """
    import dataclasses
    cfg = getattr(evaluator, "config", None)
    cites = ["standard sensitivity analysis; "
             "Bernstein et al. 2021 (doi:10.1186/s13059-021-02289-z)"]
    if cfg is None or not dataclasses.is_dataclass(cfg):
        return CheckResult("D6", "parameter sensitivity", SKIP,
                           "Needs a dataclass config to perturb.",
                           citations=cites)
    # NOTE: this rebuilds the evaluator once per parameter value, which for a
    # COBRA evaluator means reloading the model each time. It is the slowest
    # check in the suite by a wide margin. Pass fewer probe_combos, or skip
    # this check, if that cost matters.
    probe_combos = probe_combos or list(scores)[:4]
    if not probe_combos:
        return CheckResult("D6", "parameter sensitivity", SKIP,
                           "No candidates to probe.", citations=cites)

    EvalCls = type(evaluator)
    sweeps = {"ratio_match_alpha": [0.2, 0.8],
              "hard_min_growth_fraction": [0.05, 0.3]}
    rows, inert = [], []
    for param, values in sweeps.items():
        if not hasattr(cfg, param):
            continue
        base = [float(evaluator.evaluate(c, target_profile))
                for c in probe_combos]
        deltas = []
        for v in values:
            try:
                alt = EvalCls(dataclasses.replace(cfg, **{param: v}))
                got = [float(alt.evaluate(c, target_profile))
                       for c in probe_combos]
            except Exception as exc:  # noqa: BLE001
                rows.append(f"{param}={v}: could not evaluate ({exc})")
                continue
            d = float(np.max(np.abs(np.array(got) - np.array(base))))
            deltas.append(d)
            rows.append(f"{param}={v:<5} max score change {d:.3e}")
        if deltas and max(deltas) <= SOLVER_NOISE_FLOOR:
            inert.append(param)

    if not rows:
        return CheckResult("D6", "parameter sensitivity", SKIP,
                           "No recognised parameters to sweep.",
                           citations=cites)
    detail = "\n".join(rows)
    ev = {"inert_parameters": inert}
    if inert:
        return CheckResult(
            "D6", "parameter sensitivity", WARN,
            f"Changing {', '.join(inert)} does not move any score.",
            detail, "That constraint is not binding. Either it is redundant, "
            "or the solution is pinned somewhere else and the parameter can "
            "never take effect.", ev, cites)
    return CheckResult(
        "D6", "parameter sensitivity", PASS,
        "Scores respond to changes in the scoring parameters.",
        detail, evidence=ev, citations=cites)


# ---------------------------------------------------------------------------
# runner
# ---------------------------------------------------------------------------

def run_diagnostics(evaluator, target_profile, candidate_combos=None,
                     max_order=2, cap=40, verbose=True) -> DiagnosticReport:
    """Run every check and return a report.

    `evaluator` needs `.evaluate(combo, target_profile) -> float`. Checks use
    more when it is there: `.config`, `.is_valid()`, `.diagnose()`, `._model`.
    Anything missing makes the relevant check SKIP rather than fail.
    """
    if candidate_combos is None:
        candidate_combos = _default_candidates(
            evaluator, target_profile, max_order=max_order, cap=cap)
    if verbose:
        print(f"[diagnostics] scoring {len(candidate_combos)} candidates ...")
    scores = _score_all(evaluator, candidate_combos, target_profile)

    results = [
        check_null_intervention_parity(evaluator, target_profile, scores),
        check_landscape_discrimination(evaluator, target_profile, scores),
        check_scalarization_degeneracy(evaluator, target_profile),
        check_free_production(evaluator, target_profile),
        check_growth_confound(evaluator, target_profile, scores),
        check_parameter_sensitivity(evaluator, target_profile, scores),
    ]
    return DiagnosticReport(results=results, target_profile=dict(target_profile),
                            n_candidates_tested=len(candidate_combos))
