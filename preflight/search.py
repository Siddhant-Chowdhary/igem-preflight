"""
search.py - training loop, episode runner, brute force and greedy

Works with any evaluator object that has:
    evaluate(knockout_combo, target_profile) -> float

CobraEvaluator in cobra_evaluator.py is the real one. Anything with that
method works, including a few-line mock for testing the search logic on its
own without loading a model.

WHICH ONE TO USE

brute_force_best() at iGEM scale. It is exact and 175 combinations takes
under a minute on our evaluator.

greedy_search() once enumeration stops being affordable - tens of genes at
order 4 and up. It finds the exact optimum roughly two thirds of the time and
lands within about 1% when it misses.

train_agent() / run_episode() only to reproduce the benchmark. The RL agent
lost to both of the above.
"""

import itertools
import random
import numpy as np

from .config import SearchConfig
from .agent import LinearQAgent
from .filters import is_valid_combo


def sample_random_target_profile(config: SearchConfig) -> dict:
    """Random plausible target profile over the config's tracked products,
    for training diversity. Drawn from a Dirichlet so entries sum to 1."""
    products = list(config.exchange_reactions.keys())
    weights = np.random.dirichlet(np.ones(len(products)))
    return {p: round(float(w), 3) for p, w in zip(products, weights)}


def run_episode(agent: LinearQAgent, evaluator, target_profile: dict, train: bool = True):
    """
    One episode: greedily/eps-greedily build up a knockout combo one gene at
    a time, up to config.max_knockout_order genes, rewarding the agent for
    each incremental improvement in evaluator score (reward shaping -- the
    episode isn't really a long-horizon MDP, it's "sequentially building a
    set," so shaping on the delta keeps the signal dense).
    """
    config = agent.config
    combo = frozenset()
    prev_score = evaluator.evaluate(tuple(sorted(combo)), target_profile)

    trace = []
    for step in range(config.max_knockout_order):
        action = agent.select_action(target_profile, combo, greedy=not train)
        if action is None:
            break
        next_combo = frozenset(combo | {action})
        new_score = evaluator.evaluate(tuple(sorted(next_combo)), target_profile)
        reward = new_score - prev_score
        done = (step == config.max_knockout_order - 1)

        if train:
            agent.update(target_profile, combo, action, reward, next_combo, done)

        trace.append((action, new_score))
        combo, prev_score = next_combo, new_score

    return combo, prev_score, trace


def train_agent(config: SearchConfig, evaluator, n_episodes: int = 2000, seed: int = 0,
                 open_gates=frozenset(), epsilon_decay: float = 0.999, epsilon_min: float = 0.02):
    random.seed(seed)
    np.random.seed(seed)
    agent = LinearQAgent(config=config, open_gates=open_gates)
    history = []
    for _ in range(n_episodes):
        target_profile = sample_random_target_profile(config)
        _, score, _ = run_episode(agent, evaluator, target_profile, train=True)
        history.append(score)
        agent.epsilon = max(epsilon_min, agent.epsilon * epsilon_decay)
    return agent, history


def brute_force_best(config: SearchConfig, evaluator, target_profile: dict,
                      max_order: int = None, open_gates=frozenset()):
    """Exhaustively try every combo up to max_order (defaults to
    config.max_knockout_order), applying the same filter rules the agent
    uses, and return the best-scoring one. The ground-truth baseline the RL
    agent is checked against."""
    max_order = max_order or config.max_knockout_order
    best_combo, best_score = None, -1.0
    for order in range(1, max_order + 1):
        for combo in itertools.combinations(config.candidate_genes, order):
            if not is_valid_combo(frozenset(combo), target_profile,
                                   config.filter_rules, open_gates):
                continue
            score = evaluator.evaluate(combo, target_profile)
            if score > best_score:
                best_combo, best_score = combo, score
    return best_combo, best_score


def evaluate_agent_vs_bruteforce(agent: LinearQAgent, evaluator, n_test_profiles: int = 25,
                                  seed: int = 1, open_gates=frozenset(), verbose: bool = True):
    """Sanity-checks the trained agent against the brute-force optimum on
    fresh random target profiles. Returns the list of (target, rl_combo,
    rl_score, bf_combo, bf_score, gap) tuples; gap = bf_score - rl_score,
    0 meaning the agent matched brute force exactly."""
    config = agent.config
    random.seed(seed)
    np.random.seed(seed)
    rows = []
    for _ in range(n_test_profiles):
        target_profile = sample_random_target_profile(config)
        rl_combo, rl_score, _ = run_episode(agent, evaluator, target_profile, train=False)
        bf_combo, bf_score = brute_force_best(config, evaluator, target_profile,
                                               open_gates=open_gates)
        gap = bf_score - rl_score
        rows.append((target_profile, rl_combo, rl_score, bf_combo, bf_score, gap))

    if verbose:
        print(f"\n{'Target Profile':<45} {'RL combo':<25} {'RL score':<10} "
              f"{'Brute-force combo':<25} {'BF score':<10} {'Gap'}")
        print("-" * 130)
        for target_profile, rl_combo, rl_score, bf_combo, bf_score, gap in rows:
            tp_str = ", ".join(f"{k}:{v:.2f}" for k, v in target_profile.items() if v > 0.01)
            print(f"{tp_str:<45} {str(sorted(rl_combo)):<25} {rl_score:<10.3f} "
                  f"{str(sorted(bf_combo)):<25} {bf_score:<10.3f} {gap:.3f}")
        print("-" * 130)
        gaps = [r[-1] for r in rows]
        print(f"Mean gap (BF score - RL score): {np.mean(gaps):.4f}  "
              f"(closer to 0 = RL agent nearly matches brute force)")

    return rows


def greedy_search(config: SearchConfig, evaluator, target_profile: dict,
                   max_order: int = None, open_gates=frozenset()):
    """Add the single best next gene at each step, stopping when none helps.

    RECOMMENDED DEFAULT. The benchmark in benchmark.py found that on
    synthetic landscapes with main effects, sparse pairwise epistasis and a
    growth penalty, this strategy reached 99.4-99.9% of the exhaustive
    optimum across five epistasis levels and five seeds, using roughly
    2 x n_genes evaluations -- while the linear Q-learning agent in agent.py
    reached only 80-85% with much higher variance and far more evaluations.

    Use the RL agent only if you can show on your own landscape that it beats
    this. On present evidence it does not.
    """
    max_order = max_order or config.max_knockout_order
    combo = ()
    best = float(evaluator.evaluate(combo, target_profile))
    for _ in range(max_order):
        cand, cand_score = None, best
        for g in config.candidate_genes:
            if g in combo:
                continue
            trial = tuple(sorted(combo + (g,)))
            if not is_valid_combo(frozenset(trial), target_profile,
                                   config.filter_rules, open_gates):
                continue
            s = float(evaluator.evaluate(trial, target_profile))
            if s > cand_score:
                cand, cand_score = trial, s
        if cand is None:
            break
        combo, best = cand, cand_score
    return combo, best
