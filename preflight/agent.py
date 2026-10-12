"""
agent.py - linear contextual bandit / Q-learning agent

WE DO NOT RECOMMEND USING THIS. It is here so the benchmark that rejects it
can be reproduced. See benchmark.py and the README: greedy hill-climbing beat
it at every epistasis level we tested, random search on the same budget
matched greedy, and at our scale brute force beats all of them because it is
exact and cheap. Use search.greedy_search() or search.brute_force_best().

It is deliberately not deep RL. At the scale this was built for (around ten
candidate genes, three knockouts) a linear Q-function trains fast and every
weight means something, which matters when the output is supposed to justify
bench work. Everything project-specific comes from a SearchConfig.
"""

from dataclasses import dataclass, field
import random
import numpy as np

from .config import SearchConfig
from .filters import is_valid_combo


def feature_dim(config: SearchConfig) -> int:
    n_products = len(config.exchange_reactions)
    n_genes = len(config.candidate_genes)
    # [target vec][action one-hot][combo one-hot][bias][progress]
    return n_products + n_genes * 2 + 2


def build_features(config: SearchConfig, target_profile: dict,
                    current_combo: frozenset, action_gene: str) -> np.ndarray:
    """
    Feature vector for the linear Q-function. Layout (sizes depend on the
    config passed in):
      [0 : P]         target profile weights, in the order config.exchange_reactions
                       was declared
      [P : P+G]       one-hot: which gene is being considered as the action
      [P+G : P+2G]    one-hot: which genes are already in current_combo
      [P+2G]          bias term
      [P+2G+1]        len(current_combo) / max_knockout_order  (episode progress)
    """
    products = list(config.exchange_reactions.keys())
    genes = config.candidate_genes
    target_vec = np.array([target_profile.get(p, 0.0) for p in products])
    action_onehot = np.array([1.0 if g == action_gene else 0.0 for g in genes])
    combo_onehot = np.array([1.0 if g in current_combo else 0.0 for g in genes])
    bias = np.array([1.0])
    progress = np.array([len(current_combo) / config.max_knockout_order])
    return np.concatenate([target_vec, action_onehot, combo_onehot, bias, progress])


@dataclass
class LinearQAgent:
    """
    Q(s, a) = w . phi(s, a), trained with 1-step TD updates (SARSA-max,
    i.e. Q-learning). `config` supplies the gene list, max combo order,
    and feature dimension; `open_gates` controls which gated_subset filter
    rules (see filters.py) are relaxed for this agent instance -- e.g. a
    "first iteration" agent leaves a gate closed, a "later iteration" agent
    can be constructed with it open.
    """
    config: SearchConfig
    lr: float = 0.05
    gamma: float = 0.9
    epsilon: float = 0.2
    weights: np.ndarray = None
    open_gates: frozenset = field(default_factory=frozenset)

    def __post_init__(self):
        if self.weights is None:
            self.weights = np.zeros(feature_dim(self.config))

    def q_value(self, target_profile, current_combo, action_gene):
        phi = build_features(self.config, target_profile, current_combo, action_gene)
        return float(self.weights @ phi), phi

    def available_actions(self, target_profile, current_combo):
        remaining = [g for g in self.config.candidate_genes if g not in current_combo]
        valid = [
            g for g in remaining
            if is_valid_combo(frozenset(current_combo | {g}), target_profile,
                               self.config.filter_rules, self.open_gates)
        ]
        return valid

    def select_action(self, target_profile, current_combo, greedy=False):
        actions = self.available_actions(target_profile, current_combo)
        if not actions:
            return None
        if (not greedy) and random.random() < self.epsilon:
            return random.choice(actions)
        q_vals = [self.q_value(target_profile, current_combo, a)[0] for a in actions]
        return actions[int(np.argmax(q_vals))]

    def update(self, target_profile, current_combo, action_gene, reward, next_combo, done):
        q_sa, phi = self.q_value(target_profile, current_combo, action_gene)
        if done:
            target = reward
        else:
            next_actions = self.available_actions(target_profile, next_combo)
            if not next_actions:
                target = reward
            else:
                next_q = max(self.q_value(target_profile, next_combo, a)[0] for a in next_actions)
                target = reward + self.gamma * next_q
        td_error = target - q_sa
        self.weights += self.lr * td_error * phi
        return td_error
