# preflight

**Checks whether a metabolic model's scoring function is measuring anything,
before anyone acts on what it says.**

We built this for an iGEM 2026 dry-lab project engineering *S. cerevisiae*
toward a C14:0-rich fatty-acid profile, then generalised it so any team with a
COBRApy-compatible model can use it by writing a config file.

---

## Why we wrote it

Our pipeline scored gene knockouts by how closely the predicted fatty-acid
profile matched a target. Brute-force search over all 175 combinations gave us
a clear winner: `FAA1 + FAA4`, score 0.202529. The unedited strain scored
0.202509.

That is a margin of **0.010%**. GC-MS fatty-acid composition varies by roughly
1–5% between biological replicates, so no experiment we could run would tell
the winning strain apart from doing nothing. The search was exact. The
arithmetic was right. The recommendation was meaningless and nothing about the
output said so.

There are good test suites for the *model* already —
[memote](https://doi.org/10.1038/s41587-020-0446-y) checks annotations, mass
balance, stoichiometric consistency, energy-generating cycles. We ran it and
we recommend it. What none of them check is the scoring function you write on
top of the model, which is exactly where our problem was.

## The idea we would most like other people to take

**The threshold that matters is the assay's resolution, not the solver's.**

Our first version of these checks compared score differences against solver
tolerance, about one part in a million. It passed the broken pipeline above,
because 2×10⁻⁵ is comfortably bigger than 10⁻⁶. The check was useless until we
tied its thresholds to what the wet lab could actually measure.

So the thresholds here are relative and set by measurement resolution
(`DETECTABLE_RELATIVE_MARGIN` = 1%, `MEASURABLE_RELATIVE_MARGIN` = 5%). If
your assay is tighter than ours, tighten them.

## What it does

Six checks, each returning PASS, WARN or FAIL with an explanation and a
suggested fix.

| | Check | Is it new? |
| --- | --- | --- |
| D1 | Does doing nothing score as well as the best intervention? | We could not find this automated anywhere |
| D2 | Can the score separate candidates at measurable resolution? | Could not find this anywhere either |
| D3 | Does a weighted sum over several products control the blend? | No. FVA and alternate-optima analysis, plus [Das & Dennis 1997](https://doi.org/10.1007/BF01197559) |
| D4 | Are product exchange reactions realistically bounded? | No — **memote already ships this**. Ours just looks at your specific targets. Run memote too |
| D5 | Is the score tracking growth rather than product? | No ([Schuetz 2007](https://doi.org/10.1038/msb4100162)). Growth coupling is usually intended, so this never FAILs — it only asks |
| D6 | Do the scoring parameters change anything? | No. Ordinary sensitivity analysis |

**And the proposal step refuses to run when D1 or D2 fail.** It produces an
empty sheet explaining what is wrong rather than a ranking that looks
authoritative and is not. We went through OptKnock, OptGene, OptForce,
RobustKnock, cameo, MEWpy, OptFlux and StrainDesign: all of them will hand you
designs whether or not their own scoring can support one, and none checks
first.

Four of six checks are not new and we would rather say so than have someone
discover it. What we think is ours: D1, D2, and the refusal.

## It is not a search method

Brute force, greedy and reinforcement learning all answer *which candidate
scores best*. These checks answer *whether "best" means anything*. Those are
different questions and neither replaces the other.

Brute force at iGEM scale is exact and cheap — 175 combinations is under a
minute of FBA. We do not beat it and do not try to. But because brute force is
exact, its output carries more authority than any approximate method, which
makes it the one most likely to hide a broken score behind a confident winner.
That is what it did to us.

## Result on our own project

The checks fail our own shipped config. That is why we trust them.

```
D1 null-intervention parity   FAIL   best intervention beats doing nothing by 0.010%
D2 landscape discrimination   FAIL   20 candidates, all inside one measurement band
D3 scalarization degeneracy   FAIL   weighted sum gave one product at every weight
D4 free production            WARN   all four target exchanges bounded at 1000
D5 growth confound            WARN   growth identical across candidates
D6 parameter sensitivity      WARN   ratio_match_alpha moves no score

VERDICT: NOT USABLE
*** NO PROPOSALS ISSUED ***
```

The cause: Yeast9 places no realistic limit on fatty-acid export, so when
growth is allowed to drop the model already produces our exact target ratio
with no edits at all. Raising the growth floor to 20% and 40% of wild type
changed nothing. Fixing it means redefining the objective, which is a
modelling decision, so the tool reports the problem instead of hiding it.

## What the feasibility work found

Running `examples/run_feasibility.py` answered the question we should have
asked first.

| Growth held at | C14:0 export (mmol/gDW/h) | C14:0 share |
| --- | --- | --- |
| 10% of wild type | 0.211 | 100% |
| 25% | 0.178 | 100% |
| 50% | 0.119 | 100% |
| 75% | 0.060 | 100% |
| 90% | 0.024 | 100% |

Unmodified wild-type yeast already exports 100% pure C14:0 at any growth rate
up to 90% of normal. No candidate knockout improves on it. Two of them —
anything containing POX1 — take production to exactly zero.

The reason is that the model's only C14:0 route is beta-oxidative shortening
of C16:0, with POX1 catalysing the committed step. It does not make C14:0 by
terminating synthesis early. We checked whether that was an artefact (one step
runs a ligase in reverse, which is thermodynamically dubious): forbidding the
reversal dropped production 2%, and loopless FBA returns identical flux, so
the route is real within the model.

Two things follow. A plan to knock down degradation in order to preserve C14:0
is inverted relative to what this model predicts. And the model has no
reaction for the early chain release that an engineered FAS is supposed to
produce, so it cannot represent that intervention at all.

## Does a learned search pay off? No

The package contains a linear Q-learning agent because that is what we built
first. **We do not recommend using it.** Over 20 seeds on a synthetic
landscape with main effects, sparse epistasis and a growth penalty:

| At 10 genes, order ≤3 | Finds exact optimum | Mean shortfall | Evaluations |
| --- | --- | --- | --- |
| brute force | 100% (by definition) | 0% | 175 |
| greedy | 65% | 1.14% | 23 |
| reinforcement learning | 10% | 10.13% | 1,204 |

Random search on the RL agent's budget also matched greedy, which tells us the
landscape has many near-optimal answers — and where that holds, search
sophistication buys very little. Anyone reaching for machine learning on a
combinatorial design problem should test a trivial baseline first.

The landscape is synthetic because our real one is flat, and comparing search
strategies on a flat landscape measures nothing. Our agent is also a simple
one; a better-tuned agent might beat greedy and we can only say ours did not.

## Install and run

Full setup, Windows commands and troubleshooting are in
**[INSTALL.md](INSTALL.md)**. Short version:

```bash
pip install -r requirements.txt
wget https://github.com/SysBioChalmers/yeast-GEM/raw/main/model/yeast-GEM.xml
python3 tests/test_diagnostics.py      # expect: 10 passed, 0 failed
python3 examples/run_feasibility.py
```

Python 3.10+. A solver (GLPK) comes with cobra, so there is nothing else to
install. The model file is third-party and not included.

## Using it on your own project

Nothing in `preflight/` needs editing.

1. Copy `configs/yeast9_methyl_myristate.yaml`.
2. Set `model_path` and `biomass_reaction_id`.
3. List your `candidate_genes` and the `gene_to_model_id` mapping.
4. List `exchange_reactions` for your target products.
5. Write any known-invalid combinations as `pre_fba_filters` rules — three
   types, documented in `preflight/filters.py`.
6. Set the thresholds in `preflight/diagnostics.py` to your own assay's
   replicate error.

The checks only need an object with `evaluate(combo, target) -> float`. No
COBRApy requirement, no need to adopt our architecture. Checks that need more
skip rather than crash.

```python
from preflight import run_diagnostics

class MyEvaluator:
    def evaluate(self, combo, target):
        return my_existing_scoring_function(combo, target)

report = run_diagnostics(MyEvaluator(), {"product_A": 0.7, "product_B": 0.3},
                          candidate_combos=[("geneX",), ("geneY",)])
print(report.to_text())
print(report.usable)
```

## Limitations

- **Nothing here has been checked against wet-lab data.** Every number is a
  statement about the model.
- The search benchmarks are synthetic, for the reason given above.
- Our RL agent is a simple linear Q-learner.
- The thresholds assume GC-MS composition analysis, with assumed rather than
  measured replicate error.
- D5 can flag growth-coupled designs, which are often exactly what you want.
  It asks a question rather than giving a verdict.
- The wet-lab effort and testing columns in a proposal sheet are dry-lab
  estimates that had not been confirmed with anyone when we wrote them.
- We use `r_2111` as the growth reaction. Some of our earlier scripts used
  `r_4041` instead and we were unsure whether that mattered. It does not: in
  Yeast9 the biomass metabolite has exactly one producer (`r_4041`) and one
  consumer (`r_2111`), so mass balance forces identical flux through both at
  every solution. We checked numerically as well — constraining either to 50%
  of wild-type growth and maximising C14:0 gives 0.119191 either way.

## Tests

```bash
python3 -m pytest tests/
```

26 tests, none needing the model file. `test_diagnostics.py` includes the case
that defeated our first attempt: an evaluator whose scores differ by 1e-5,
which is real arithmetic and useless biology.

## License

MIT — see `LICENSE`.
