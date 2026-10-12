"""
run_feasibility.py - what can the model actually reach?

This is the analysis we should have run first. Instead of asking which
knockout is best, it asks what C14:0 yield and purity the model can reach at
a growth rate the cell could survive, and whether any knockout moves that.

    python examples/run_feasibility.py

You need yeast-GEM.xml in the working directory:
    wget https://github.com/SysBioChalmers/yeast-GEM/raw/main/model/yeast-GEM.xml

Three parts:
  1. the growth/yield frontier for each candidate knockout set
  2. which route the model uses to make C14:0
  3. whether that route survives a thermodynamic constraint

A reminder that applies to everything this prints: FBA reports what is
POSSIBLE given the stoichiometry and the bounds you set. It is an upper
limit. It is not a prediction, it is not a probability, and it has no enzyme
kinetics or regulation in it.
"""
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from cobra.io import read_sbml_model

from preflight import load_config

CONFIG = os.path.join(os.path.dirname(__file__), "..", "configs",
                       "yeast9_methyl_myristate.yaml")

GROWTH_FRACTIONS = [0.01, 0.10, 0.25, 0.50, 0.75, 0.90]
COMBOS = [(), ("FAA1",), ("FAA1", "FAA4"), ("DGA1", "FAA1", "LRO1"),
          ("POX1", "FAA1"), ("POX1", "PXA1", "FAA1")]


def load():
    cfg = load_config(CONFIG)
    model = read_sbml_model(cfg.model_path)
    wt = model.optimize().objective_value
    print(f"Wild-type growth: {wt:.6f} /h\n")
    return cfg, model, wt


def frontier(cfg, model, wt, combo, growth_fractions=GROWTH_FRACTIONS):
    """Max C14:0 export at each growth floor, then other fatty acids
    minimised so the reported share is the purest profile achievable."""
    c14_id = cfg.exchange_reactions["C14:0"]
    out = []
    with model as m:
        for g in combo:
            mid = cfg.gene_to_model_id.get(g)
            if mid:
                try:
                    m.genes.get_by_id(mid).knock_out()
                except KeyError:
                    pass
        for gf in growth_fractions:
            with m as mm:
                mm.reactions.get_by_id(
                    cfg.biomass_reaction_id).lower_bound = gf * wt
                mm.objective = mm.reactions.get_by_id(c14_id)
                mm.objective.direction = "max"
                c14 = mm.slim_optimize()
                if c14 is None or c14 != c14:
                    out.append({"growth_frac": gf, "feasible": False})
                    continue
                if c14 < 1e-9:
                    out.append({"growth_frac": gf, "feasible": True,
                                "c14": 0.0, "share": 0.0})
                    continue
                mm.reactions.get_by_id(c14_id).lower_bound = 0.99 * c14
                others = [mm.reactions.get_by_id(r)
                          for n, r in cfg.exchange_reactions.items()
                          if n != "C14:0"]
                mm.objective = {r: 1.0 for r in others}
                mm.objective.direction = "min"
                sol = mm.optimize()
                if sol.status != "optimal":
                    out.append({"growth_frac": gf, "feasible": True,
                                "c14": float(c14), "share": None})
                    continue
                flux = {n: max(0.0, sol.fluxes.get(r, 0.0))
                        for n, r in cfg.exchange_reactions.items()}
                tot = sum(flux.values())
                out.append({
                    "growth_frac": gf, "feasible": True,
                    "c14": float(flux["C14:0"]), "total": float(tot),
                    "share": float(flux["C14:0"] / tot) if tot > 1e-12 else 0.0,
                })
    return out


def part1_frontiers(cfg, model, wt):
    print("=" * 76)
    print("1. GROWTH / YIELD FRONTIER")
    print("=" * 76)
    print(f"{'knockouts':<24} {'growth':>7} {'C14:0 mmol/gDW/h':>18} "
          f"{'C14:0 share':>12}")
    print("-" * 76)
    results = {}
    for combo in COMBOS:
        rows = frontier(cfg, model, wt, combo)
        results[" ".join(combo) or "(none)"] = rows
        label = " ".join(combo) or "(no edit)"
        for i, r in enumerate(rows):
            if not r.get("feasible"):
                print(f"{label if i == 0 else '':<24} "
                      f"{r['growth_frac']*100:>6.0f}%   INFEASIBLE")
                continue
            share = r.get("share")
            print(f"{label if i == 0 else '':<24} "
                  f"{r['growth_frac']*100:>6.0f}% {r.get('c14', 0):>18.4f} "
                  f"{(share if share is not None else 0)*100:>11.1f}%")
        print()
    return results


def part2_route(cfg, model, wt, growth_frac=0.50):
    print("=" * 76)
    print(f"2. WHICH ROUTE MAKES THE C14:0 (at {growth_frac*100:.0f}% WT growth)")
    print("=" * 76)
    c14_id = cfg.exchange_reactions["C14:0"]
    with model as m:
        m.reactions.get_by_id(
            cfg.biomass_reaction_id).lower_bound = growth_frac * wt
        m.objective = m.reactions.get_by_id(c14_id)
        m.objective.direction = "max"
        sol = m.optimize()
        print(f"C14:0 export: {sol.fluxes[c14_id]:.4f} mmol/gDW/h\n")
        print("Peroxisomal reactions carrying flux (|v| > 0.01):")
        shown = 0
        for r in m.reactions:
            v = sol.fluxes.get(r.id, 0.0)
            if abs(v) > 0.01 and any(c == "p" for c in r.compartments):
                print(f"  {r.id:10} {v:>9.4f}  {r.name[:52]}")
                shown += 1
                if shown >= 12:
                    print("  ...")
                    break
    print()


def part3_robustness(cfg, model, wt, growth_frac=0.50):
    """One step in the route runs a fatty-acid--CoA ligase in reverse, which
    is thermodynamically questionable. Does production depend on it?"""
    print("=" * 76)
    print("3. IS THE ROUTE A MODELLING ARTEFACT?")
    print("=" * 76)
    c14_id = cfg.exchange_reactions["C14:0"]
    REVERSIBLE = ("r_0412", "r_2196", "r_2201", "r_0402")

    def run(label, forward_only=False, ko=()):
        with model as m:
            for g in ko:
                mid = cfg.gene_to_model_id.get(g)
                if mid:
                    try:
                        m.genes.get_by_id(mid).knock_out()
                    except KeyError:
                        pass
            if forward_only:
                for rid in REVERSIBLE:
                    try:
                        rx = m.reactions.get_by_id(rid)
                        if rx.lower_bound < 0:
                            rx.lower_bound = 0
                    except KeyError:
                        pass
            m.reactions.get_by_id(
                cfg.biomass_reaction_id).lower_bound = growth_frac * wt
            m.objective = m.reactions.get_by_id(c14_id)
            m.objective.direction = "max"
            v = m.slim_optimize()
        v = 0.0 if (v is None or v != v) else float(v)
        print(f"  {label:<50} {v:>9.4f}")
        return v

    print(f"Max C14:0 export at >={growth_frac*100:.0f}% WT growth:")
    base = run("as shipped (ligases reversible)")
    fwd = run("ligases constrained forward-only", forward_only=True)
    run("POX1 deleted", ko=("POX1",))
    run("POX1 deleted AND ligases forward-only", forward_only=True, ko=("POX1",))

    try:
        from cobra.flux_analysis import loopless_solution
        with model as m:
            m.reactions.get_by_id(
                cfg.biomass_reaction_id).lower_bound = growth_frac * wt
            m.objective = m.reactions.get_by_id(c14_id)
            m.objective.direction = "max"
            s = m.optimize()
            ll = loopless_solution(m)
        print(f"\n  standard FBA {s.fluxes[c14_id]:.4f} vs "
              f"loopless {ll.fluxes[c14_id]:.4f}")
    except Exception as exc:  # noqa: BLE001
        print(f"\n  loopless check unavailable: {exc}")

    verdict = ("COLLAPSES -- the route depends on the reversal"
               if fwd < 0.01 * max(base, 1e-9)
               else "SURVIVES -- the route is not a reversal artefact")
    print(f"\n  Verdict: {verdict}\n")


def main():
    cfg, model, wt = load()
    frontiers = part1_frontiers(cfg, model, wt)
    part2_route(cfg, model, wt)
    part3_robustness(cfg, model, wt)

    out = os.path.join(os.path.dirname(__file__), "feasibility_frontier.json")
    with open(out, "w") as f:
        json.dump({"wt_growth": wt, "frontiers": frontiers}, f, indent=1)
    print(f"Raw numbers written to {out}")
    print("\nReminder: FBA gives an upper limit on what is possible, not a")
    print("prediction and not a probability.")


if __name__ == "__main__":
    main()
