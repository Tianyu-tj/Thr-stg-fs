"""
simulate_demo.py - synthetic-data experiments for the three-stage framework and OAENet.

Two modes
---------
--mode paper   (default) replicates FS_CI_Final.R / Section 6.1 of the three-stage paper:
    X ~ N(0, Sigma), Sigma_jk = rho, p = 20, X standardised after generation
    T = 1{U < expit(X beta)},  Y = X theta + N(0, 1) + 0.5 T
    coefficients multiplied by --coef-mult (2 in FS_CI_Final.R; 1 = values printed in Section 6.1)
    methods: Enh-ELRT, Enh-ESVMS, ONet (benchmark of FS_CI_Final.R) in paper mode
    metric : percentage selection bias |ATT - ATT_target| / 0.5 * 100 and per-variable selection rates
--mode general  the earlier demonstration (simulate_demo_general.py; continuous and binary outcomes).

Confounders X1-X2, pure outcome predictors X3-X4, pure treatment predictors X5-X6, noise X7-X20.

--from-r DIR  runs the paper-mode methods on datasets exported by simulate_demo.R (export_dir) and
              compares the selected variables with R's, data set by data set.

Usage:
  python simulate_demo.py --mode paper --coef-mult 2 --n 200 --reps 30
  python simulate_demo.py --mode paper --from-r r_export
"""
import argparse, glob, os
import numpy as np
import pandas as pd
from scipy.special import expit
from causal_fs import three_stage, oaenet, att_matching

SCENARIOS = {   # (theta: outcome, beta: treatment), as printed in Section 6.1
    1: (np.r_[.6, .6, .6, .6, np.zeros(16)], np.r_[1, 1, 0, 0, 1, 1, np.zeros(14)]),
    2: (np.r_[.6, .6, .6, .6, np.zeros(16)], np.r_[.4, .4, 0, 0, 1, 1, np.zeros(14)]),
    3: (np.r_[.2, .2, .6, .6, np.zeros(16)], np.r_[1, 1, 0, 0, 1, 1, np.zeros(14)]),
    4: (np.r_[.6, .6, .6, .6, np.zeros(16)], np.r_[1, 1, 0, 0, 1.8, 1.8, np.zeros(14)]),
}
PAPER_METHODS = {
    "Enh-ELRT":  lambda X, Y, T, s: three_stage(X, Y, T, "LR", "tanh", seed=s, mode="paper")["selected"],
    "Enh-ESVMS": lambda X, Y, T, s: three_stage(X, Y, T, "SVM", "sigmoid", seed=s, mode="paper")["selected"],
    "ONet":      lambda X, Y, T, s: oaenet(X, Y, T, mode="paper")["selected"],
}


def simulate_paper(n, rho, scenario, seed, coef_mult=2.0, p=20):
    rng = np.random.default_rng(seed)
    theta, beta = SCENARIOS[scenario]
    theta, beta = coef_mult * theta[:p], coef_mult * beta[:p]
    S = np.full((p, p), rho); np.fill_diagonal(S, 1.0)
    X = rng.multivariate_normal(np.zeros(p), S, size=n)
    T = (rng.uniform(size=n) < expit(X @ beta)).astype(int)
    Y = X @ theta + rng.normal(size=n) + 0.5 * T
    X = (X - X.mean(0)) / X.std(0, ddof=1)
    return X, Y, T


def record(X, Y, T, sel, **keys):
    att = att_matching(X[:, sel], Y, T, mode="paper")
    tgt = att_matching(X[:, :4], Y, T, mode="paper")
    r = dict(keys, n_selected=len(sel), att=att, att_target=tgt, pct_bias=100 * abs(att - tgt) / 0.5)
    for j in range(X.shape[1]):
        r[f"X{j + 1}"] = int(j in set(np.asarray(sel).tolist()))
    return r


def summarise(res, by):
    res = res.copy()
    res["XC"] = res[["X1", "X2"]].mean(axis=1); res["XP"] = res[["X3", "X4"]].mean(axis=1)
    res["XT"] = res[["X5", "X6"]].mean(axis=1)
    res["XS"] = res[[c for c in res.columns if c.startswith("X") and c[1:].isdigit() and int(c[1:]) > 6]].mean(axis=1)
    g = res.groupby(by)
    out = g[["XC", "XP", "XT", "XS", "pct_bias"]].mean().round(3)
    out["n_failed"] = g.att.apply(lambda s: int(s.isna().sum()))
    return out


def run_paper(a):
    rows = []
    for sc in a.scenarios:
        for rho in a.rho:
            for seed in range(1, a.reps + 1):
                X, Y, T = simulate_paper(a.n, rho, sc, seed, a.coef_mult)
                for name, f in PAPER_METHODS.items():
                    rows.append(record(X, Y, T, f(X, Y, T, seed), scenario=sc, rho=rho, seed=seed, method=name))
            print(f"scenario {sc} rho {rho} done", flush=True)
            pd.DataFrame(rows).to_csv(a.out, index=False)
    res = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print(summarise(res, ["scenario", "rho", "method"]))
    print("\nPooled over rho (Figure 3 / Figure 4 layout):")
    print(summarise(res, ["scenario", "method"]))


def run_from_r(a):
    """Datasets data_*.csv (columns X1..Xp, A, Y) and r_selections.csv exported by simulate_demo.R."""
    rsel = pd.read_csv(os.path.join(a.from_r, "r_selections.csv"))
    rows = []
    for f in sorted(glob.glob(os.path.join(a.from_r, "data_*.csv"))):
        key = os.path.basename(f)[5:-4]
        d = pd.read_csv(f)
        X = d.filter(regex=r"^X\d+$").to_numpy(float); T = d["A"].to_numpy(int); Y = d["Y"].to_numpy(float)
        seed = int(key.split("_itr")[-1])
        for name, fn in PAPER_METHODS.items():
            sel = set((np.asarray(fn(X, Y, T, seed)) + 1).tolist())
            r_row = rsel[(rsel.dataset == key) & (rsel.method == name)]
            rs = set(int(v) for v in str(r_row.selected.iloc[0]).split(";") if v not in ("", "nan")) if len(r_row) else None
            rows.append(dict(dataset=key, method=name, python=";".join(map(str, sorted(sel))),
                             r=None if rs is None else ";".join(map(str, sorted(rs))),
                             identical=None if rs is None else sel == rs))
    res = pd.DataFrame(rows); res.to_csv(a.out, index=False)
    print(res.groupby("method").identical.mean().rename("share of datasets with identical selection"))
    print(res[res.identical == False].head(20).to_string(index=False))


def run_general():
    import sys
    from simulate_demo_general import main as general_main   # earlier demo, kept unchanged
    argv = sys.argv[1:]
    if "--mode" in argv:                                     # pass every other argument through
        i = argv.index("--mode"); del argv[i:i + 2]
    sys.argv = [sys.argv[0]] + argv
    general_main()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["paper", "general"], default="paper")
    ap.add_argument("--coef-mult", type=float, default=2.0)
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--rho", type=float, nargs="+", default=[0, 0.25, 0.5, 0.75])
    ap.add_argument("--scenarios", type=int, nargs="+", default=[1, 2, 3, 4])
    ap.add_argument("--from-r", default=None)
    ap.add_argument("--out", default="simulation_results_paper.csv")
    a, _ = ap.parse_known_args()
    if a.from_r:
        run_from_r(a)
    elif a.mode == "paper":
        run_paper(a)
    else:
        run_general()
