"""
simulate_demo.py - synthetic-data demonstration of the three-stage framework and OAENet.

Data-generating process (Section 6.1 of the three-stage paper; Shortreed & Ertefaie 2017):
    X ~ N(0, Sigma), Sigma_jk = rho (j != k), p = 20
    T ~ Bernoulli(expit(X beta))
    continuous outcome: Y = 0.5 T + X theta + eps,            eps ~ N(0, 1)
    binary outcome:     Y ~ Bernoulli(expit(-1 + 0.5 T + X theta))
Confounders X1-X2, pure outcome predictors X3-X4, pure treatment predictors X5-X6, noise X7-X20.
The desired selection is {X1, X2, X3, X4}.

Usage:  python simulate_demo.py --reps 10 --n 500
"""
import argparse
import numpy as np
import pandas as pd
from scipy.special import expit
from causal_fs import three_stage, oaenet, att_matching

SCENARIOS = {
    1: (np.r_[.6, .6, .6, .6, np.zeros(16)], np.r_[1, 1, 0, 0, 1, 1, np.zeros(14)]),
    2: (np.r_[.6, .6, .6, .6, np.zeros(16)], np.r_[.4, .4, 0, 0, 1, 1, np.zeros(14)]),
    3: (np.r_[.2, .2, .6, .6, np.zeros(16)], np.r_[1, 1, 0, 0, 1, 1, np.zeros(14)]),
    4: (np.r_[.6, .6, .6, .6, np.zeros(16)], np.r_[1, 1, 0, 0, 1.8, 1.8, np.zeros(14)]),
}


def simulate(n, rho, scenario, seed, outcome="continuous", p=20):
    rng = np.random.default_rng(seed)
    theta, beta = SCENARIOS[scenario]
    S = np.full((p, p), rho); np.fill_diagonal(S, 1.0)
    X = rng.multivariate_normal(np.zeros(p), S, size=n)
    T = rng.binomial(1, expit(X @ beta))
    if outcome == "continuous":
        Y = 0.5 * T + X @ theta + rng.normal(size=n)
        true_att = 0.5
    else:
        eta0 = -1 + X @ theta
        Y = rng.binomial(1, expit(eta0 + 0.5 * T)).astype(float)
        true_att = float(np.mean(expit(eta0[T == 1] + 0.5) - expit(eta0[T == 1])))
    return X, Y, T, true_att


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--rho", type=float, nargs="+", default=[0.25, 0.75])
    ap.add_argument("--scenarios", type=int, nargs="+", default=[1, 2, 3, 4])
    ap.add_argument("--outcomes", nargs="+", default=["continuous", "binary"])
    a = ap.parse_args()

    methods = {
        "Enh-ESVMS": lambda X, Y, T, s: three_stage(X, Y, T, "SVM", "sigmoid", seed=s, verbose=False)["selected"],
        "Enh-ELRT":  lambda X, Y, T, s: three_stage(X, Y, T, "LR", "tanh", seed=s, verbose=False)["selected"],
        "OAENet":    lambda X, Y, T, s: oaenet(X, Y, T, verbose=False)["selected"],
    }
    rows = []
    for outcome in a.outcomes:
        for sc in a.scenarios:
            for rho in a.rho:
                for seed in range(1, a.reps + 1):
                    X, Y, T, true_att = simulate(a.n, rho, sc, seed, outcome)
                    target = att_matching(X[:, :4], Y, T)                       # target model: X1-X4
                    for name, f in methods.items():
                        sel = set(f(X, Y, T, seed).tolist())
                        att = att_matching(X[:, sorted(sel)], Y, T) if sel else np.nan
                        rows.append(dict(outcome=outcome, scenario=sc, rho=rho, seed=seed, method=name,
                                         X1_X4_all=int({0, 1, 2, 3} <= sel), X5_X6_any=int(bool(sel & {4, 5})),
                                         n_noise=len(sel - set(range(6))), att=att, target_att=target, true_att=true_att,
                                         abs_bias_vs_target=abs(att - target), abs_bias_vs_truth=abs(att - true_att)))
                print(f"{outcome:10s} scenario {sc} done", flush=True)
    res = pd.DataFrame(rows)
    res.to_csv("simulation_results.csv", index=False)
    summary = res.groupby(["outcome", "method", "scenario"]).agg(
        X1_X4_all=("X1_X4_all", "mean"), X5_X6_any=("X5_X6_any", "mean"), n_noise=("n_noise", "mean"),
        abs_bias_vs_target=("abs_bias_vs_target", "mean"), abs_bias_vs_truth=("abs_bias_vs_truth", "mean")).round(3)
    pd.set_option("display.width", 200)
    print(summary)
    summary.to_csv("simulation_summary.csv")


if __name__ == "__main__":
    main()
