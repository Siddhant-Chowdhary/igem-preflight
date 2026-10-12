"""
config.py - loads a project's genes, reactions, filters and model path

A SearchConfig holds everything specific to one project, organism or target.
The checks, the evaluator and the search code never hardcode biology; they
read it from here. To use this on something else, copy
configs/yeast9_methyl_myristate.yaml, change the values, and leave this
directory alone.
"""

from dataclasses import dataclass, field
import yaml

from .filters import load_rules, Rule


@dataclass
class SearchConfig:
    # --- model / organism ---
    model_path: str
    biomass_reaction_id: str

    # --- candidate gene space ---
    candidate_genes: list[str]
    gene_to_model_id: dict[str, str]   # e.g. common name -> systematic ORF ID
    max_knockout_order: int

    # --- target product(s) ---
    exchange_reactions: dict[str, str]  # e.g. "C14:0" -> "r_2193"

    # --- search-space pruning (see filters.py) ---
    filter_rules: list[Rule]

    # --- scoring / objective shape ---
    ratio_match_alpha: float = 0.5       # min fraction of max total output required
    hard_min_growth_fraction: float = 0.01
    growth_penalty_K: float = 0.5
    growth_penalty_n: float = 2.0

    # --- metadata (not used by the engine, just for provenance) ---
    name: str = "unnamed-project"
    description: str = ""
    sources: list[str] = field(default_factory=list)


def load_config(path: str) -> SearchConfig:
    with open(path) as f:
        raw = yaml.safe_load(f)

    required = ["model_path", "biomass_reaction_id", "candidate_genes",
                "gene_to_model_id", "exchange_reactions", "max_knockout_order"]
    missing = [k for k in required if k not in raw]
    if missing:
        raise ValueError(f"Config {path} is missing required keys: {missing}")

    return SearchConfig(
        model_path=raw["model_path"],
        biomass_reaction_id=raw["biomass_reaction_id"],
        candidate_genes=raw["candidate_genes"],
        gene_to_model_id=raw["gene_to_model_id"],
        max_knockout_order=raw["max_knockout_order"],
        exchange_reactions=raw["exchange_reactions"],
        filter_rules=load_rules(raw.get("pre_fba_filters", [])),
        ratio_match_alpha=raw.get("ratio_match_alpha", 0.5),
        hard_min_growth_fraction=raw.get("hard_min_growth_fraction", 0.01),
        growth_penalty_K=raw.get("growth_penalty_K", 0.5),
        growth_penalty_n=raw.get("growth_penalty_n", 2.0),
        name=raw.get("name", "unnamed-project"),
        description=raw.get("description", ""),
        sources=raw.get("sources", []),
    )
