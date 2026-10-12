# Scoring-setup diagnostic report

**Not usable.** A blocking check failed, so any ranking this setup produces is meaningless.

Target profile: `{'C14:0': 0.6, 'C16:0': 0.2, 'C18:0': 0.1, 'C18:1': 0.1}`  
Candidates tested: 20

| Check | Result | Finding |
| --- | --- | --- |
| D1 null-intervention parity | FAIL | The best intervention beats doing nothing by only 0.010%, which is below what the assay could resolve, so the ranking is untestable even though the arithmetic is sound. |
| D2 landscape discrimination | FAIL | The solver separates these candidates into 2 levels, but all of them fall inside one measurement-resolution band, so no experiment could confirm the ordering. |
| D3 scalarization degeneracy | FAIL | A weighted sum over these products cannot control the blend: at every weight tested the solver produced just one product. (This tests the weighted-sum formulation itself, not your current objective.) |
| D4 free production | WARN | 4 target exchange reaction(s) are effectively unbounded, so the model may produce the target more freely than a cell could. |
| D5 growth confound | WARN | Neither the product-match nor the growth term varies meaningfully across candidates. |
| D6 parameter sensitivity | WARN | Changing ratio_match_alpha does not move any score. |

## Detail

### D1 null-intervention parity — FAIL
The best intervention beats doing nothing by only 0.010%, which is below what the assay could resolve, so the ranking is untestable even though the arithmetic is sound.

```
no-intervention score : 0.202509
best candidate        : ['FAA1', 'FAA4'] -> 0.202529
margin                : 1.973e-05  (0.010% of baseline)
```

**Fix:** The objective is essentially satisfiable without any edit. Constrain the route the product must take (e.g. require flux through the synthesis pathway), bound the product exchange realistically, or score product yield at fixed growth instead of profile shape.

*References: StrainDesign docs note the wild type can be returned as a solution (doi:10.1093/bioinformatics/btac632)*

### D2 landscape discrimination — FAIL
The solver separates these candidates into 2 levels, but all of them fall inside one measurement-resolution band, so no experiment could confirm the ordering.

```
candidates              : 20
score range             : 0.202509 .. 0.202529
spread                  : 1.973e-05 (0.010% of typical score)
std deviation           : 4.300e-06
tiers the solver can see: 2 of 20
tiers an assay can see  : 1 of 20
```

**Fix:** The ranking is numerically real and experimentally meaningless. Redefine the score so differences exceed assay error, or report a single undifferentiated set.

*References: Bernstein et al. 2021 (doi:10.1186/s13059-021-02289-z); Kaste & Shachar-Hill 2024 (arXiv:2303.12651)*

### D3 scalarization degeneracy — FAIL
A weighted sum over these products cannot control the blend: at every weight tested the solver produced just one product. (This tests the weighted-sum formulation itself, not your current objective.)

```
weights [0.25, 0.25, 0.25, 0.25] -> 1 product(s) active: ['C14:0']
weights [0.39, 0.59, 0.01, 0.0] -> 1 product(s) active: ['C16:0']
weights [0.15, 0.45, 0.19, 0.21] -> 1 product(s) active: ['C16:0']
weights [0.23, 0.5, 0.27, 0.0] -> 1 product(s) active: ['C16:0']
weights [0.53, 0.02, 0.25, 0.2] -> 1 product(s) active: ['C14:0']
```

**Fix:** Do not use a weighted sum for a ratio target. Use lexicographic / goal programming: maximize total output first, then minimize deviation from the target ratio. Or inspect the trade-off surface with COBRApy's production_envelope.

*References: Das & Dennis 1997 (doi:10.1007/BF01197559); Mahadevan & Schilling 2003 (doi:10.1016/j.ymben.2003.09.002); Reed & Palsson 2004 (doi:10.1101/gr.2546004)*

### D4 free production — WARN
4 target exchange reaction(s) are effectively unbounded, so the model may produce the target more freely than a cell could.

```
C14:0    r_2193     bounds [0, 1000]   <- effectively unbounded
C16:0    r_1993     bounds [0, 1000]   <- effectively unbounded
C18:0    r_2055     bounds [0, 1000]   <- effectively unbounded
C18:1    r_2189     bounds [0, 1000]   <- effectively unbounded
```

**Fix:** Set product exchange upper bounds from measured or literature secretion rates, and run memote's unbounded-flux and energy-cycle checks on the model.

*References: memote (doi:10.1038/s41587-020-0446-y) -- run it as well; Fritzemeier et al. 2017 (doi:10.1371/journal.pcbi.1005494); Thiele & Palsson 2010 (doi:10.1038/nprot.2009.203)*

### D5 growth confound — WARN
Neither the product-match nor the growth term varies meaningfully across candidates.

```
product-match (cosine) spread : 2.116e-05 (0.002% of its scale, inert)
growth-fraction spread        : 8.231e-15 (0.000% of its scale, inert)
corr(score, growth)           : nan
corr(score, product match)    : 1.000
```

**Fix:** Nothing in the score responds to the intervention. See D1 and D2.

*References: Schuetz et al. 2007 (doi:10.1038/msb4100162); Feist & Palsson 2010 (doi:10.1016/j.mib.2010.03.003); Klamt & Mahadevan 2015 (doi:10.1016/j.ymben.2015.05.006) -- growth coupling is often intended*

### D6 parameter sensitivity — WARN
Changing ratio_match_alpha does not move any score.

```
ratio_match_alpha=0.2   max score change 0.000e+00
ratio_match_alpha=0.8   max score change 9.910e-12
hard_min_growth_fraction=0.05  max score change 1.364e-02
hard_min_growth_fraction=0.3   max score change 1.341e-01
```

**Fix:** That constraint is not binding. Either it is redundant, or the solution is pinned somewhere else and the parameter can never take effect.

*References: standard sensitivity analysis; Bernstein et al. 2021 (doi:10.1186/s13059-021-02289-z)*
