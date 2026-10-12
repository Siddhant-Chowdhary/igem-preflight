"""
Tests for the combination filters.

These lock in the behaviour that used to be hardcoded if-statements in our
earlier scripts, so that reorganising filters.py later cannot quietly change
which combinations count as valid.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from preflight.filters import load_rules, is_valid_combo


def make_rules():
    raw = [
        {"type": "mutual_exclusion_at_order", "genes": ["FAA1", "FAA4"], "min_order": 3},
        {"type": "conditional_exclusion", "gene": "OLE1",
         "condition": {"target_profile_key": "C18:1", "operator": ">", "value": 0}},
        {"type": "gated_subset", "genes": ["POX1", "PXA1", "PXA2"],
         "gate_name": "allow_degradation_knockouts", "default": False},
    ]
    return load_rules(raw)


def test_mutual_exclusion_blocks_triple():
    rules = make_rules()
    combo = frozenset({"FAA1", "FAA4", "DGA1"})
    assert is_valid_combo(combo, {}, rules) is False


def test_mutual_exclusion_allows_pair():
    rules = make_rules()
    combo = frozenset({"FAA1", "FAA4"})
    assert is_valid_combo(combo, {}, rules) is True


def test_conditional_exclusion_blocks_ole1_when_target_wants_c18_1():
    rules = make_rules()
    combo = frozenset({"OLE1"})
    assert is_valid_combo(combo, {"C18:1": 0.3}, rules) is False


def test_conditional_exclusion_allows_ole1_when_target_has_no_c18_1():
    rules = make_rules()
    combo = frozenset({"OLE1"})
    assert is_valid_combo(combo, {"C18:1": 0.0}, rules) is True


def test_gated_subset_blocks_by_default():
    rules = make_rules()
    combo = frozenset({"POX1"})
    assert is_valid_combo(combo, {}, rules) is False


def test_gated_subset_allows_when_gate_open():
    rules = make_rules()
    combo = frozenset({"POX1"})
    assert is_valid_combo(combo, {}, rules, open_gates={"allow_degradation_knockouts"}) is True


def test_empty_combo_always_valid():
    rules = make_rules()
    assert is_valid_combo(frozenset(), {"C18:1": 0.5}, rules) is True


def test_unrelated_combo_valid():
    rules = make_rules()
    combo = frozenset({"DGA1", "LRO1"})
    assert is_valid_combo(combo, {"C18:1": 0.5}, rules) is True


def test_unknown_rule_type_raises():
    try:
        import pytest
        with pytest.raises(ValueError):
            load_rules([{"type": "not_a_real_rule"}])
    except ImportError:
        try:
            load_rules([{"type": "not_a_real_rule"}])
            raise AssertionError("expected ValueError")
        except ValueError:
            pass


if __name__ == "__main__":
    # Lightweight fallback so this still runs without pytest installed.
    tests = [v for k, v in list(globals().items()) if k.startswith("test_")]
    passed, failed = 0, 0
    for t in tests:
        try:
            t()
            passed += 1
        except Exception as e:
            failed += 1
            print(f"FAILED: {t.__name__}: {e}")
    print(f"\n{passed} passed, {failed} failed")
