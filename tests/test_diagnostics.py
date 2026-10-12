"""
Tests for the checks. Each one pins down a failure mode we actually hit.

Run:  python -m pytest tests/      or:  python tests/test_diagnostics.py

These use tiny fake evaluators with deliberately broken scoring, so they run
in milliseconds and need no model file. The one worth knowing about is
_Imperceptible: its scores differ by 1e-5, which is real arithmetic and
useless biology. Our first version of the gate let that through, which is
exactly how we shipped a broken pipeline, so it is a test now.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from preflight.diagnostics import (
    run_diagnostics, check_null_intervention_parity,
    check_landscape_discrimination, PASS, WARN, FAIL,
)


class _Flat:
    """Every candidate scores identically -- the worst case."""
    def evaluate(self, combo, target):
        return 0.5


class _Imperceptible:
    """Scores differ, but by far less than any assay could measure.
    This is the case that slipped through our first, solver-tolerance-based
    version of the gate."""
    def evaluate(self, combo, target):
        return 0.5 + 1e-5 * len(combo)


class _Good:
    """Clear, measurable separation, with per-gene effects so the landscape
    has more than two levels. An earlier version of this fixture scored only
    by combo size, which produced two levels -- and D2 rightly called that a
    WARN, since a two-level landscape cannot order candidates within a tier."""
    _W = {"A": 0.20, "B": 0.13, "C": 0.07}

    def evaluate(self, combo, target):
        return 0.30 + sum(self._W.get(g, 0.0) for g in combo)


COMBOS = [(), ("A",), ("B",), ("A", "B"), ("C",), ("A", "C")]


def _scores(ev):
    return {c: ev.evaluate(c, {}) for c in COMBOS if c != ()}


def test_flat_landscape_fails_d1():
    r = check_null_intervention_parity(_Flat(), {}, _scores(_Flat()))
    assert r.status == FAIL
    assert r.blocking


def test_flat_landscape_fails_d2():
    r = check_landscape_discrimination(_Flat(), {}, _scores(_Flat()))
    assert r.status == FAIL
    assert r.blocking


def test_imperceptible_margin_fails_d1():
    """A 0.002% margin is real arithmetic and useless biology."""
    r = check_null_intervention_parity(
        _Imperceptible(), {}, _scores(_Imperceptible()))
    assert r.status == FAIL, f"expected FAIL, got {r.status}: {r.headline}"


def test_imperceptible_margin_fails_d2():
    r = check_landscape_discrimination(
        _Imperceptible(), {}, _scores(_Imperceptible()))
    assert r.status == FAIL


def test_good_landscape_passes_both():
    assert check_null_intervention_parity(
        _Good(), {}, _scores(_Good())).status == PASS
    assert check_landscape_discrimination(
        _Good(), {}, _scores(_Good())).status == PASS


def test_report_marks_flat_setup_unusable():
    rep = run_diagnostics(_Flat(), {}, candidate_combos=COMBOS[1:],
                           verbose=False)
    assert not rep.usable
    assert len(rep.blocking_failures) >= 1


def test_report_marks_good_setup_usable():
    rep = run_diagnostics(_Good(), {}, candidate_combos=COMBOS[1:],
                           verbose=False)
    assert rep.usable


def test_missing_optional_hooks_skip_not_crash():
    """An evaluator with only evaluate() must still produce a report."""
    rep = run_diagnostics(_Good(), {}, candidate_combos=COMBOS[1:],
                           verbose=False)
    assert len(rep.results) == 6
    assert rep.to_text()
    assert rep.to_markdown()


def test_proposal_sheet_refuses_on_flat_landscape():
    from preflight.proposals import build_proposal_sheet
    sheet = build_proposal_sheet(_Flat(), {}, COMBOS[1:], top_n=3)
    assert sheet.gated
    assert sheet.proposals == []
    assert "NO PROPOSALS ISSUED" in sheet.to_text()


def test_proposal_sheet_emits_on_good_landscape():
    from preflight.proposals import build_proposal_sheet
    sheet = build_proposal_sheet(_Good(), {}, COMBOS[1:], top_n=3)
    assert not sheet.gated
    assert len(sheet.proposals) == 3
    assert sheet.proposals[0].score >= sheet.proposals[-1].score


if __name__ == "__main__":
    tests = [v for k, v in list(globals().items()) if k.startswith("test_")]
    p = f = 0
    for t in tests:
        try:
            t(); p += 1
        except Exception as e:  # noqa: BLE001
            f += 1; print(f"FAILED {t.__name__}: {e}")
    print(f"\n{p} passed, {f} failed")
