"""
benchmark.py - we checked whether our own RL agent was worth it

The obvious objection to using reinforcement learning for knockout search at
our scale is that it is unnecessary. With 10 candidate genes and up to 3
knockouts there are only 175 combinations, so you can just try them all. We
thought the objection deserved a real test rather than an argument, so we
built one.

WHY THE LANDSCAPE HERE IS SYNTHETIC

We could not run this on our own FBA landscape, because that landscape turned
out to be flat (see diagnostics.py - doing nothing scored as well as the best
knockout). Comparing search strategies on a flat landscape measures nothing.
So this builds an explicit synthetic landscape with main effects, sparse
pairwise epistasis, a growth penalty, and a known best answer.

What that does show: how each strategy's cost and solution quality change as
the candidate gene set grows.

What it does not show: that any of this transfers to a real genome-scale
model. That depends on the real landscape having similar structure, and for
our project we have not established that. Anyone quoting these numbers needs
to carry that caveat with them.

WHAT WE FOUND

Greedy hill-climbing beat our RL agent at every epistasis level we tried, and
random search on the same evaluation budget matched greedy. At our actual
scale brute force beats both, because it is exact and takes under a minute.
The numbers are in the README. We left the RL agent in the repository so this
benchmark can be reproduced, not because we recommend using it.
"""

from dataclasses import dataclass, field
import itertools
import random
import time

import numpy as np


# ---------------------------------------------------------------------------
# synthetic landscape
# ---------------------------------------------------------------------------

@dataclass
class SyntheticLandscape:
    """A deterministic fitness landscape over gene knockouts.

    Structure, chosen to mimic the qualitative features of a product-ratio
    FBA objective rather than any specific organism:

      * each gene has a main effect on the product vector (some genes push
        toward short chains, some toward long, some do nothing much)
      * each gene carries a growth cost
      * selected gene PAIRS interact (epistasis): their combined effect is
        not the sum of their parts, which is what makes the landscape worth
        searching rather than sorting
      * score = cosine(product vector, target) * growth penalty

    Every evaluation is counted, so search strategies can be compared on
    evaluation budget rather than wall-clock time.
    """
    n_genes: int
    n_products: int = 4
    seed: int = 0
    epistasis_pairs: int = 0          # 0 -> n_genes (one per gene, on average)
    epistasis_strength: float = 0.6
    growth_cost_scale: float = 0.08

    genes: list = field(default_factory=list, init=False)
    n_evals: int = field(default=0, init=False)

    def __post_init__(self):
        rng = np.random.default_rng(self.seed)
        self.genes = [f"G{i:02d}" for i in range(self.n_genes)]
        self._rng = rng
        # main effect of deleting each gene on the product vector
        self._main = rng.normal(0, 1, size=(self.n_genes, self.n_products))
        # growth cost of deleting each gene (always non-negative)
        self._cost = np.abs(rng.normal(0, 1, size=self.n_genes)) * \
            self.growth_cost_scale
        # baseline product vector of the unedited strain
        self._base = np.abs(rng.normal(1.0, 0.25, size=self.n_products))
        # pairwise interactions on a sparse set of pairs
        n_pairs = self.epistasis_pairs or self.n_genes
        all_pairs = list(itertools.combinations(range(self.n_genes), 2))
        rng.shuffle(all_pairs)
        self._pairs = {}
        for (i, j) in all_pairs[:min(n_pairs, len(all_pairs))]:
            self._pairs[(i, j)] = rng.normal(
                0, self.epistasis_strength, size=self.n_products)
        self._cache = {}

    # -- the evaluator interface the search code expects --------------------

    def evaluate(self, knockout_combo, target_profile: dict) -> float:
        key = (tuple(sorted(knockout_combo)),
               tuple(sorted(target_profile.items())))
        if key in self._cache:
            return self._cache[key]
        self.n_evals += 1

        idx = [self.genes.index(g) for g in knockout_combo if g in self.genes]
        vec = self._base.copy()
        for i in idx:
            vec = vec + self._main[i]
        for a, b in itertools.combinations(sorted(idx), 2):
            if (a, b) in self._pairs:
                vec = vec + self._pairs[(a, b)]
        vec = np.clip(vec, 0.0, None)

        growth = float(np.clip(1.0 - sum(self._cost[i] for i in idx), 0.0, 1.0))

        t = np.array([target_profile.get(p, 0.0)
                      for p in self.product_names()])
        if np.linalg.norm(vec) < 1e-12 or np.linalg.norm(t) < 1e-12:
            score = 0.0
        else:
            cos = float(np.dot(vec, t) /
                        (np.linalg.norm(vec) * np.linalg.norm(t)))
            # Hill-shaped growth penalty, same form as the FBA evaluator
            d = 1.0 - growth
            penalty = (d ** 2) / (0.5 ** 2 + d ** 2)
            score = cos * (1.0 - penalty)

        self._cache[key] = score
        return score

    def product_names(self):
        return [f"P{i}" for i in range(self.n_products)]

    def random_target(self, rng=None):
        rng = rng or self._rng
        w = rng.dirichlet(np.ones(self.n_products))
        return {p: round(float(x), 3)
                for p, x in zip(self.product_names(), w)}

    def reset_counter(self):
        self.n_evals = 0
        self._cache = {}


# ---------------------------------------------------------------------------
# search strategies, all on the same evaluation-budget footing
# ---------------------------------------------------------------------------

def search_exhaustive(landscape, target, max_order):
    """Every combination up to max_order. Exact, and the cost baseline."""
    best, best_score = None, -np.inf
    for order in range(1, max_order + 1):
        for combo in itertools.combinations(landscape.genes, order):
            s = landscape.evaluate(combo, target)
            if s > best_score:
                best, best_score = combo, s
    return best, best_score


def search_random(landscape, target, max_order, budget, seed=0):
    """Random valid combos, same evaluation budget as the learned search."""
    rng = random.Random(seed)
    best, best_score = None, -np.inf
    for _ in range(budget):
        order = rng.randint(1, max_order)
        combo = tuple(rng.sample(landscape.genes, order))
        s = landscape.evaluate(combo, target)
        if s > best_score:
            best, best_score = combo, s
    return best, best_score


def search_greedy(landscape, target, max_order):
    """Add the single best next gene each step. Cheap, and the strategy
    epistasis is specifically expected to defeat."""
    combo, best_score = (), landscape.evaluate((), target)
    for _ in range(max_order):
        cand, cand_score = None, best_score
        for g in landscape.genes:
            if g in combo:
                continue
            s = landscape.evaluate(combo + (g,), target)
            if s > cand_score:
                cand, cand_score = combo + (g,), s
        if cand is None:
            break
        combo, best_score = cand, cand_score
    return combo, best_score


def search_rl(landscape, target_sampler, eval_target, max_order,
               n_episodes=300, seed=0, lr=0.08, epsilon=0.25, gamma=0.9):
    """Linear Q-learning over (target profile, partial combo, next gene).

    Deliberately the same algorithm class as the project's own agent, written
    standalone here so the benchmark does not depend on a SearchConfig or a
    COBRA model.
    """
    rng = random.Random(seed)
    nprng = np.random.default_rng(seed)
    genes = landscape.genes
    n_p = landscape.n_products
    dim = n_p + 2 * len(genes) + 2
    w = np.zeros(dim)

    def phi(target, combo, action):
        tv = np.array([target.get(p, 0.0) for p in landscape.product_names()])
        a = np.array([1.0 if g == action else 0.0 for g in genes])
        c = np.array([1.0 if g in combo else 0.0 for g in genes])
        return np.concatenate([tv, a, c, [1.0], [len(combo) / max_order]])

    def q(target, combo, action):
        p = phi(target, combo, action)
        return float(w @ p), p

    eps = epsilon
    for _ in range(n_episodes):
        target = target_sampler(nprng)
        combo = ()
        prev = landscape.evaluate(combo, target)
        for step in range(max_order):
            avail = [g for g in genes if g not in combo]
            if not avail:
                break
            if rng.random() < eps:
                action = rng.choice(avail)
            else:
                action = max(avail, key=lambda g: q(target, combo, g)[0])
            nxt = combo + (action,)
            s = landscape.evaluate(nxt, target)
            reward = s - prev
            done = step == max_order - 1
            q_sa, p = q(target, combo, action)
            if done:
                tgt = reward
            else:
                rest = [g for g in genes if g not in nxt]
                tgt = reward + (gamma * max(q(target, nxt, g)[0] for g in rest)
                                if rest else 0.0)
            w += lr * (tgt - q_sa) * p
            combo, prev = nxt, s
        eps = max(0.02, eps * 0.995)

    # greedy rollout on the target we actually care about
    combo = ()
    best_score = landscape.evaluate(combo, eval_target)
    for _ in range(max_order):
        avail = [g for g in genes if g not in combo]
        if not avail:
            break
        action = max(avail, key=lambda g: q(eval_target, combo, g)[0])
        nxt = combo + (action,)
        s = landscape.evaluate(nxt, eval_target)
        if s > best_score:
            combo, best_score = nxt, s
        else:
            combo = nxt
            best_score = max(best_score, s)
    return combo, best_score


# ---------------------------------------------------------------------------
# the scaling experiment
# ---------------------------------------------------------------------------

@dataclass
class BenchRow:
    n_genes: int
    max_order: int
    strategy: str
    best_score: float
    frac_of_exhaustive: float
    n_evals: int
    seconds: float
    combo: tuple


def combination_count(n_genes, max_order):
    from math import comb
    return sum(comb(n_genes, k) for k in range(1, max_order + 1))


def run_scaling_benchmark(gene_counts=(10, 20, 30, 40, 50), max_order=3,
                           n_episodes=300, seed=0, fba_seconds_per_eval=0.30,
                           verbose=True):
    """Compare strategies as the candidate gene set grows.

    `fba_seconds_per_eval` converts evaluation counts into the wall-clock
    cost they would carry on a real genome-scale model. 0.30 s is measured
    from this project's own 3-stage lexicographic Yeast9 evaluator; set it to
    your own measured figure before quoting any projected time.
    """
    rows = []
    for n in gene_counts:
        land = SyntheticLandscape(n_genes=n, seed=seed)
        target = land.random_target(np.random.default_rng(seed + 1))
        total_combos = combination_count(n, max_order)

        # exhaustive, the exact reference
        land.reset_counter()
        t0 = time.time()
        ex_combo, ex_score = search_exhaustive(land, target, max_order)
        ex_evals, ex_sec = land.n_evals, time.time() - t0
        rows.append(BenchRow(n, max_order, "exhaustive", ex_score, 1.0,
                             ex_evals, ex_sec, ex_combo))

        # greedy
        land.reset_counter()
        t0 = time.time()
        g_combo, g_score = search_greedy(land, target, max_order)
        rows.append(BenchRow(n, max_order, "greedy", g_score,
                             g_score / ex_score if ex_score else 0.0,
                             land.n_evals, time.time() - t0, g_combo))

        # RL
        land.reset_counter()
        t0 = time.time()
        rl_combo, rl_score = search_rl(
            land, land.random_target, target, max_order,
            n_episodes=n_episodes, seed=seed)
        rl_evals, rl_sec = land.n_evals, time.time() - t0
        rows.append(BenchRow(n, max_order, "rl", rl_score,
                             rl_score / ex_score if ex_score else 0.0,
                             rl_evals, rl_sec, rl_combo))

        # random, matched to the RL evaluation budget -- the control that
        # separates "learning helped" from "sampling that many points helped"
        land.reset_counter()
        t0 = time.time()
        r_combo, r_score = search_random(land, target, max_order,
                                          budget=rl_evals, seed=seed)
        rows.append(BenchRow(n, max_order, "random_matched", r_score,
                             r_score / ex_score if ex_score else 0.0,
                             land.n_evals, time.time() - t0, r_combo))

        if verbose:
            print(f"n={n:3d}  combos={total_combos:>8,}  "
                  f"exhaustive={ex_score:.4f} ({ex_evals:,} evals)  "
                  f"rl={rl_score:.4f} ({rl_evals:,})  "
                  f"greedy={g_score:.4f}  random={r_score:.4f}", flush=True)
    return rows


def benchmark_table(rows, fba_seconds_per_eval=0.30) -> str:
    """Plain-text table, including what each evaluation count would cost on a
    real GEM at the given seconds-per-evaluation."""
    def hours(n):
        h = n * fba_seconds_per_eval / 3600.0
        if h < 1 / 60:
            return f"{n * fba_seconds_per_eval:.0f}s"
        if h < 1:
            return f"{h*60:.0f}min"
        return f"{h:.1f}h"

    L = [f"{'genes':>5} {'strategy':<16} {'score':>8} {'% of exact':>11} "
         f"{'evals':>9} {'projected FBA time':>19}",
         "-" * 74]
    for r in rows:
        L.append(f"{r.n_genes:>5} {r.strategy:<16} {r.best_score:>8.4f} "
                 f"{r.frac_of_exhaustive*100:>10.1f}% {r.n_evals:>9,} "
                 f"{hours(r.n_evals):>19}")
    L.append("-" * 74)
    L.append(f"projected time = evals x {fba_seconds_per_eval}s, this "
             "project's measured Yeast9 evaluator cost")
    return "\n".join(L)


def crossover_summary(rows, fba_seconds_per_eval=0.30,
                       practical_hours=8.0) -> str:
    """Where does enumeration stop being the right answer?

    `practical_hours` is the wall-clock budget past which exhaustive search
    stops being something a student can run between lab sessions.
    """
    by_n = {}
    for r in rows:
        by_n.setdefault(r.n_genes, {})[r.strategy] = r
    L = []
    crossover = None
    for n in sorted(by_n):
        ex = by_n[n].get("exhaustive")
        rl = by_n[n].get("rl")
        if not ex or not rl:
            continue
        ex_h = ex.n_evals * fba_seconds_per_eval / 3600.0
        if crossover is None and ex_h > practical_hours:
            crossover = n
        L.append(f"n={n:3d}: exhaustive needs {ex.n_evals:,} evaluations "
                 f"({ex_h:.1f} h on a real GEM); RL used {rl.n_evals:,} "
                 f"({rl.n_evals*fba_seconds_per_eval/3600.0:.1f} h) and "
                 f"reached {rl.frac_of_exhaustive*100:.1f}% of the exact "
                 "optimum")
    head = (f"Enumeration passes the {practical_hours:.0f}-hour practical "
            f"budget at about {crossover} candidate genes."
            if crossover else
            f"Enumeration stayed within the {practical_hours:.0f}-hour budget "
            "at every size tested, so a learned search is not needed in this "
            "range.")
    return head + "\n" + "\n".join(L)
