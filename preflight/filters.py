"""
filters.py - rules for combinations that are not worth running FBA on

Anything we know in advance about which knockout combinations are lethal or
contradictory lives here as DATA, read from a config file, rather than as
if-statements buried in the search code. That is what lets another team reuse
this: they edit their config, not this file.

Three rule types covered everything our project needed.

  mutual_exclusion_at_order
    If all of these genes are knocked out together and the combination has
    reached at least this size, skip it. We use it for a gene pair that is
    lethal once a third deletion is added.

  conditional_exclusion
    Skip a gene whenever something about the target profile makes deleting it
    self-contradictory. We use it so the only desaturase is not deleted when
    the target still asks for the unsaturated product it makes.

  gated_subset
    These genes are off limits unless a named gate is opened. We use it to
    keep degradation-pathway genes out of first-iteration proposals, which
    was a decision our team made about what should reach the wet lab and
    when. The gate defaults closed so that behaviour is the default for
    anyone reusing the config.

No eval() and no string-executed conditions anywhere. Conditions are small
structured comparisons so a config file cannot turn into arbitrary code.
"""

from dataclasses import dataclass
from typing import Any


@dataclass
class Rule:
    kind: str
    params: dict


_OPERATORS = {
    ">": lambda a, b: a > b,
    ">=": lambda a, b: a >= b,
    "<": lambda a, b: a < b,
    "<=": lambda a, b: a <= b,
    "==": lambda a, b: a == b,
    "!=": lambda a, b: a != b,
}


def load_rules(raw_rules: list[dict]) -> list[Rule]:
    """Parse the `pre_fba_filters` section of a config file into Rule objects."""
    rules = []
    for r in raw_rules:
        kind = r.get("type")
        if kind not in ("mutual_exclusion_at_order", "conditional_exclusion", "gated_subset"):
            raise ValueError(
                f"Unknown filter rule type '{kind}'. Must be one of: "
                "mutual_exclusion_at_order, conditional_exclusion, gated_subset"
            )
        rules.append(Rule(kind=kind, params=r))
    return rules


def _check_mutual_exclusion_at_order(combo: frozenset, params: dict) -> bool:
    """True = combo is VALID (not excluded)."""
    genes = set(params["genes"])
    min_order = params["min_order"]
    if genes.issubset(combo) and len(combo) >= min_order:
        return False
    return True


def _check_conditional_exclusion(combo: frozenset, target_profile: dict, params: dict) -> bool:
    """True = combo is VALID (not excluded)."""
    gene = params["gene"]
    if gene not in combo:
        return True
    cond = params["condition"]
    key = cond["target_profile_key"]
    op = cond["operator"]
    value = cond["value"]
    actual = target_profile.get(key, 0)
    if op not in _OPERATORS:
        raise ValueError(f"Unknown operator '{op}' in conditional_exclusion rule")
    condition_true = _OPERATORS[op](actual, value)
    # If the condition holds, this gene's knockout is excluded.
    return not condition_true


def _check_gated_subset(combo: frozenset, params: dict, open_gates: set) -> bool:
    """True = combo is VALID (not excluded)."""
    genes = set(params["genes"])
    gate_name = params["gate_name"]
    default_open = params.get("default", False)
    is_open = gate_name in open_gates if open_gates else default_open
    if not is_open and genes & combo:
        return False
    return True


def is_valid_combo(combo: frozenset, target_profile: dict, rules: list[Rule],
                    open_gates: set | None = None) -> bool:
    """
    Run every rule against a candidate combo. Returns False (invalid) on the
    first rule that rejects it, matching the original project's short-circuit
    behaviour. `open_gates` is the set of gate names a caller has explicitly
    opened (e.g. {"allow_degradation_knockouts"} for a later search
    iteration) -- omit it to use each gate's own config-specified default.
    """
    open_gates = open_gates or set()
    for rule in rules:
        if rule.kind == "mutual_exclusion_at_order":
            if not _check_mutual_exclusion_at_order(combo, rule.params):
                return False
        elif rule.kind == "conditional_exclusion":
            if not _check_conditional_exclusion(combo, target_profile, rule.params):
                return False
        elif rule.kind == "gated_subset":
            if not _check_gated_subset(combo, rule.params, open_gates):
                return False
    return True
