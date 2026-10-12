"""
proposals.py - wet-lab proposal sheets that refuse to print when the score
is not good enough to justify them

Every strain-design tool we looked at (OptKnock, OptGene, OptForce,
RobustKnock, cameo, MEWpy, OptFlux, StrainDesign) will give you a ranked list
of designs whether or not its own scoring can support one. None of them check
first and none of them refuse.

This one refuses. It runs diagnostics.py, and if a blocking check fails it
produces an empty sheet explaining what is wrong instead of a ranking that
looks authoritative and is not. If we are going to ask the wet lab for three
weeks of cloning, we should be able to say why.

WHAT WE ASSUMED ABOUT THE WET LAB

We assumed a standard yeast workflow: knockouts built by marker-based
homologous recombination, fatty-acid profile read by GC-MS of derivatised
total lipid extract, reported as a fraction of total peak area. We had not
confirmed this when we wrote the sheets, so the "how to test" and "effort"
columns are our estimate and not a commitment from anyone. If the real
workflow is different, those columns need rewriting.

CONFIDENCE TIERS

  A  The score separates this design from the alternatives by more than the
     assay could resolve, and the direction of the effect matches what we
     would expect from what the genes do.
  B  The score separates it, but we cannot explain why from first
     principles. Still worth testing, but a negative result would not tell
     us much about the model.
  C  Inside the noise. Listed so the sheet is complete. Do not spend bench
     time on these because of us.

Tier A is still a model prediction. None of this has been checked against
wet-lab data.
"""

from dataclasses import dataclass, field
import csv

from .diagnostics import (
    run_diagnostics, DiagnosticReport,
    DETECTABLE_RELATIVE_MARGIN, MEASURABLE_RELATIVE_MARGIN,
)


@dataclass
class Proposal:
    rank: int
    genes: tuple
    score: float
    margin_over_null: float
    tier: str
    rationale: str
    expected_signature: str
    falsified_if: str
    effort: str

    def as_row(self) -> dict:
        return {
            "rank": self.rank,
            "knockouts": " ".join(self.genes) if self.genes else "(none)",
            "n_genes": len(self.genes),
            "score": f"{self.score:.6f}",
            "margin_over_no_edit": f"{self.margin_over_null:.3e}",
            "confidence_tier": self.tier,
            "why": self.rationale,
            "expected_GCMS_signature": self.expected_signature,
            "model_falsified_if": self.falsified_if,
            "construction_effort": self.effort,
        }


@dataclass
class ProposalSheet:
    target_profile: dict
    diagnostics: DiagnosticReport
    proposals: list = field(default_factory=list)
    gated: bool = False
    gate_reason: str = ""
    assumptions: list = field(default_factory=list)

    def to_text(self) -> str:
        w = 78
        L = ["=" * w, "WET-LAB PROPOSAL SHEET", "=" * w,
             f"Target profile: " +
             ", ".join(f"{k}={v}" for k, v in self.target_profile.items()), ""]
        if self.gated:
            L += ["*** NO PROPOSALS ISSUED ***", "",
                  "The scoring setup did not pass its own diagnostics, so any",
                  "ranking it produced would be an artefact. Rather than hand",
                  "the wet lab a list we cannot stand behind, this sheet is",
                  "deliberately empty.", "", "Reason:", ""]
            for ln in self.gate_reason.split("\n"):
                L.append(f"  {ln}")
            L += ["", "Fix the blocking items in the diagnostic report, then",
                  "re-run. The diagnostic report is attached below.", ""]
        else:
            L += [f"{len(self.proposals)} proposal(s), highest confidence "
                  "first.", ""]
            for p in self.proposals:
                L += [f"--- Rank {p.rank}  [tier {p.tier}]  "
                      f"{' '.join(p.genes) if p.genes else '(none)'}",
                      f"    score            : {p.score:.6f}  "
                      f"(margin over no-edit {p.margin_over_null:.3e})",
                      f"    why              : {p.rationale}",
                      f"    expect to see    : {p.expected_signature}",
                      f"    model wrong if   : {p.falsified_if}",
                      f"    build effort     : {p.effort}", ""]
        if self.assumptions:
            L += ["-" * w, "ASSUMPTIONS THE DRY LAB MADE (please correct)", ""]
            L += [f"  - {a}" for a in self.assumptions]
            L.append("")
        L += ["-" * w, "DIAGNOSTIC REPORT", ""]
        L.append(self.diagnostics.to_text())
        return "\n".join(L)

    def to_csv(self, path: str) -> str:
        cols = ["rank", "knockouts", "n_genes", "score", "margin_over_no_edit",
                "confidence_tier", "why", "expected_GCMS_signature",
                "model_falsified_if", "construction_effort"]
        with open(path, "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=cols)
            wr.writeheader()
            if self.gated:
                wr.writerow({c: "" for c in cols} | {
                    "knockouts": "NO PROPOSALS ISSUED",
                    "why": self.gate_reason.replace("\n", " ")})
            for p in self.proposals:
                wr.writerow(p.as_row())
        return path


# ---------------------------------------------------------------------------
# rationale / expectation text
# ---------------------------------------------------------------------------

# Plain-language roles, so a sheet reads the same to a wet-lab member as to a
# modeller. Keys are gene names; anything absent falls back to a neutral line.
GENE_ROLES = {
    "POX1": ("acyl-CoA oxidase, first committed step of peroxisomal "
             "beta-oxidation", "less fatty-acid breakdown"),
    "PXA1": ("peroxisomal ABC transporter subunit importing fatty acids for "
             "breakdown", "less fatty acid entering breakdown"),
    "PXA2": ("the other half of that transporter",
             "less fatty acid entering breakdown"),
    "FAA1": ("acyl-CoA synthetase that re-activates free fatty acids",
             "more free fatty acid left unactivated and available for export"),
    "FAA4": ("second acyl-CoA synthetase with longer-chain preference",
             "more free fatty acid left unactivated"),
    "DGA1": ("diacylglycerol acyltransferase, makes storage triacylglycerol",
             "less carbon diverted into storage lipid"),
    "LRO1": ("phospholipid-diacylglycerol acyltransferase, also makes storage "
             "triacylglycerol", "less carbon diverted into storage lipid"),
    "ARE1": ("sterol acyltransferase making steryl esters",
             "less carbon into sterol-linked storage"),
    "ARE2": ("the main sterol acyltransferase",
             "less carbon into sterol-linked storage"),
    "OLE1": ("delta-9 desaturase, the only route to unsaturated fatty acids",
             "loss of unsaturated species, which is lethal without supplement"),
}


def _rationale(genes) -> str:
    if not genes:
        return "No edit; included as the baseline comparison."
    parts = []
    for g in genes:
        role, effect = GENE_ROLES.get(
            g, (f"{g} (role not annotated in this config)", "an unknown effect"))
        parts.append(f"{g} is {role}, so deleting it should give {effect}")
    return "; ".join(parts) + "."


def _top_products(target_profile, n=1):
    """Highest-weighted products in the target, or [] for an empty target."""
    return [k for k, _ in
            sorted(target_profile.items(), key=lambda kv: -kv[1])[:n]]


def _signature(genes, target_profile) -> str:
    if not genes:
        return "The unedited profile, by definition."
    top = _top_products(target_profile, 2)
    if not top:
        return ("No target products were specified, so no expected signature "
                "can be stated.")
    want = " and ".join(top)
    return (f"A shift in the GC-MS profile toward {want} relative to the "
            "parent strain, measured as a fraction of total fatty-acid peak "
            "area, not absolute titre.")


def _falsified(genes, target_profile) -> str:
    if not genes:
        return "Not applicable."
    top = _top_products(target_profile, 1)
    if not top:
        return ("No target products were specified, so no falsification "
                "criterion can be stated.")
    top = top[0]
    return (f"The {top} fraction does not increase relative to the parent, or "
            "total fatty-acid output collapses while the ratio shifts, which "
            "would mean the model traded away yield that it should have kept.")


def _effort(genes) -> str:
    n = len(genes)
    if n == 0:
        return "None; this is the parent strain."
    if n == 1:
        return "One marker-based deletion; standard."
    if n == 2:
        return "Two sequential deletions; needs two markers or marker recycling."
    return (f"{n} sequential deletions; consider whether the marker budget and "
            "timeline justify it before starting.")


# ---------------------------------------------------------------------------
# builder
# ---------------------------------------------------------------------------

def build_proposal_sheet(evaluator, target_profile, candidate_combos,
                          top_n=5, diagnostics_report=None,
                          extra_assumptions=None,
                          mechanism_expected=None) -> ProposalSheet:
    """Score candidates, run diagnostics, and either rank or refuse.

    `candidate_combos` is whatever your search produced (or an enumeration).
    `mechanism_expected` is an optional set of gene names whose direction of
    effect you are willing to defend mechanistically; combos drawn entirely
    from it can reach tier A. Leave it None to let every discriminated combo
    be tier B, which is the conservative default.
    """
    report = diagnostics_report or run_diagnostics(
        evaluator, target_profile,
        candidate_combos=[tuple(c) for c in candidate_combos],
        verbose=False)

    assumptions = [
        "Knockouts are built by standard marker-based homologous "
        "recombination in the project's parent strain.",
        "Fatty-acid profile is read by GC-MS of derivatised total lipid "
        "extract, reported as fraction of total peak area.",
        "The wet lab's actual throughput and profiling method were not "
        "confirmed with them; effort and test columns are the dry lab's "
        "estimate.",
        "Point mutations in the project's mutation table are NOT applied by "
        "this evaluator, so these predictions describe knockout effects only.",
    ] + list(extra_assumptions or [])

    if not report.usable:
        reasons = []
        for r in report.blocking_failures:
            reasons.append(f"{r.check_id} {r.name}: {r.headline}")
            if r.fix:
                reasons.append(f"  -> {r.fix}")
        return ProposalSheet(
            target_profile=dict(target_profile), diagnostics=report,
            proposals=[], gated=True, gate_reason="\n".join(reasons),
            assumptions=assumptions)

    null = float(evaluator.evaluate((), target_profile))
    scored = sorted(
        ((tuple(c), float(evaluator.evaluate(tuple(c), target_profile)))
         for c in candidate_combos), key=lambda kv: -kv[1])

    proposals = []
    scale = max(abs(null), 1e-12)
    for i, (genes, score) in enumerate(scored[:top_n], start=1):
        margin = score - null
        rel = margin / scale
        # Tiers are set by whether an experiment could resolve the margin, not
        # by whether the solver could. See diagnostics.py on why.
        if rel < DETECTABLE_RELATIVE_MARGIN:
            tier = "C"
        elif rel < MEASURABLE_RELATIVE_MARGIN:
            tier = "B"
        elif mechanism_expected and set(genes) <= set(mechanism_expected):
            tier = "A"
        else:
            tier = "B"
        proposals.append(Proposal(
            rank=i, genes=genes, score=score, margin_over_null=margin,
            tier=tier, rationale=_rationale(genes),
            expected_signature=_signature(genes, target_profile),
            falsified_if=_falsified(genes, target_profile),
            effort=_effort(genes)))

    return ProposalSheet(target_profile=dict(target_profile),
                          diagnostics=report, proposals=proposals,
                          gated=False, assumptions=assumptions)
