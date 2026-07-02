"""
biosafety_lhs.py — Pillar 4c: Global uncertainty quantification for the
Luria-Delbrück evolutionary biosafety model via Latin Hypercube Sampling.

The single-point P_escape estimate of Pillar 4 (6.0e-17 at the
assumed baseline parameters) does not propagate the joint uncertainty
in the input parameters. Pillar 4b (biosafety_sensitivity.py) addresses
this with a one-factor-at-a-time local sweep, but a local sweep does
not capture joint-tail sensitivity — a parameter can be fine at ±10x
locally and still fail badly in combination with other uncertain
parameters at the joint tails.

This script replaces the local sweep with a global LHS analysis:

  - 8 input parameters are sampled log-uniformly over literature-
    bounded ranges (Drake 1991, Lee 2012 for mutation rate; Licht
    1999, Dahlberg 1998 for HGT; deployment-scale ranges for
    reactor volume, cell density, generation time, deployment
    window, and kill-switch target size).
  - N = 10,000 LHS samples are drawn using scipy.stats.qmc.LatinHypercube
    with a fixed random seed (42) for reproducibility.
  - For each sample, P_escape is computed via the exact Poisson
    formula P = 1 - exp(-N*G*p), evaluated using scipy.special.expm1
    to avoid the IEEE-754 catastrophic cancellation that the
    Taylor-approximation branch in biosafety_mutation_model.py
    works around. expm1(x) computes exp(x) - 1 with full precision
    even for |x| ~ 1e-17, so P = -expm1(-N*G*p) is the exact,
    IEEE-754-safe expression for the Poisson escape probability.
  - Outputs: a CSV of all 10,000 samples; a 3-panel figure (CDF,
    histogram, Spearman correlation heatmap); and a stdout summary
    reporting median, 95% CI, fraction of samples above the
    de-minimis threshold, and the top-3 parameters most strongly
    correlated with P_escape.

All numerical constants and ranges are imported from parameters.py.

Outputs:
  - results/biosafety_lhs.csv       (10,000 rows)
  - results/biosafety_lhs.png       (3-panel figure)
  - stdout summary                  (median, 95% CI, top correlations)
"""

import csv
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import qmc, spearmanr
from scipy.special import expm1

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parameters import (
    # LHS configuration
    LHS_N_SAMPLES,
    LHS_RANDOM_SEED,
    LHS_DISTRIBUTION,
    # LHS parameter bounds (log10 values)
    LHS_LOG10_MUTATION_RATE_PER_BP_PER_GEN,
    LHS_LOG10_HGT_P_REVERT_THYA_PER_DIV,
    LHS_LOG10_HGT_P_REVERT_DAPA_PER_DIV,
    LHS_LOG10_POPULATION_VOLUME_L,
    LHS_LOG10_MAX_CELL_DENSITY_PER_L,
    LHS_LOG10_GENERATIONS_PER_DAY,
    LHS_LOG10_DEPLOYMENT_WINDOW_DAYS,
    LHS_LOG10_TARGET_SIZE_BP,
    # Threshold for reporting
    DE_MINIMIS_RISK_THRESHOLD,
)


# Parameter names in fixed order, matching the LHS column layout
PARAM_NAMES = [
    "mu_bp_per_gen",            # 1: per-bp per-gen mutation rate
    "p_HGT_thyA_per_div",       # 2: HGT reversion, thyA
    "p_HGT_dapA_per_div",       # 3: HGT reversion, dapA
    "V_reactor_L",              # 4: reactor volume
    "rho_max_per_L",            # 5: max cell density
    "g_per_day",                # 6: generations per day
    "T_days",                   # 7: deployment window
    "L_target_bp",              # 8: kill-switch target size
]

# Log10 bounds in the same order
LOG10_BOUNDS = np.array([
    LHS_LOG10_MUTATION_RATE_PER_BP_PER_GEN,
    LHS_LOG10_HGT_P_REVERT_THYA_PER_DIV,
    LHS_LOG10_HGT_P_REVERT_DAPA_PER_DIV,
    LHS_LOG10_POPULATION_VOLUME_L,
    LHS_LOG10_MAX_CELL_DENSITY_PER_L,
    LHS_LOG10_GENERATIONS_PER_DAY,
    LHS_LOG10_DEPLOYMENT_WINDOW_DAYS,
    LHS_LOG10_TARGET_SIZE_BP,
])

# Pre-flight check on bounds
assert LOG10_BOUNDS.shape == (8, 2), f"Expected 8 (lo,hi) bounds, got {LOG10_BOUNDS.shape}"
assert np.all(LOG10_BOUNDS[:, 0] < LOG10_BOUNDS[:, 1]), "Each lo bound must be < hi bound"
assert LHS_DISTRIBUTION == "log-uniform", f"Only log-uniform supported; got {LHS_DISTRIBUTION!r}"


def _p_escape_exact(N, G, p):
    """Exact Poisson escape probability P = 1 - exp(-N*G*p).

    Uses scipy.special.expm1 for IEEE-754-safe evaluation. expm1(x)
    computes exp(x) - 1 with full precision even for |x| ~ 1e-17,
    so the formula P = -expm1(-N*G*p) is exact for the entire
    range of arguments arising in this study (N*G*p in [~1e-30, ~1e-5]).

    This is the mathematically exact replacement for the manual
    Taylor-approximation branch in biosafety_mutation_model.py's
    _prob_escape(), which uses np.where(x < 1e-10, x, 1-exp(-x))
    to work around the same cancellation issue.

    Accepts numpy arrays; broadcasts N and p with G.
    """
    x = N * G * p
    # expm1 handles x ~ 0 exactly; for large x, -expm1(-x) -> 1 - exp(-x) -> 1
    return -expm1(-x)


def run_lhs():
    """Draw LHS samples, compute P_escape for each, summarize, and
    write CSV + 3-panel figure.

    Returns
    -------
    summary : dict
        Keys: median, ci_lo, ci_hi, frac_above_threshold,
        top_corr_param_1, top_corr_param_2, top_corr_param_3,
        rho_1, rho_2, rho_3.
    """
    # -----------------------------------------------------------------
    # 1. Draw LHS samples in [0, 1]^8 and invert through the
    #    log10-bounds to get the parameter values.
    # -----------------------------------------------------------------
    sampler = qmc.LatinHypercube(d=8, seed=LHS_RANDOM_SEED)
    unit_samples = sampler.random(n=LHS_N_SAMPLES)  # shape (N, 8), values in [0, 1)

    # Convert unit-space to log10 space, then to linear space
    lo = LOG10_BOUNDS[:, 0]
    hi = LOG10_BOUNDS[:, 1]
    log10_samples = qmc.scale(unit_samples, lo, hi)  # shape (N, 8)
    samples = 10.0 ** log10_samples                  # shape (N, 8)

    # Unpack for clarity
    mu_bp      = samples[:, 0]   # per-bp per-gen
    p_thyA     = samples[:, 1]   # per-division
    p_dapA     = samples[:, 2]
    V          = samples[:, 3]   # L
    rho_max    = samples[:, 4]   # cells/L
    g_per_day  = samples[:, 5]
    T_days     = samples[:, 6]
    L_target   = samples[:, 7]   # bp

    # -----------------------------------------------------------------
    # 2. Compute P_escape for each sample using the exact Poisson
    #    formula with IEEE-754-safe expm1.
    # -----------------------------------------------------------------
    N = V * rho_max                 # total cell count
    G = g_per_day * T_days          # total generations
    p_ks_fail = mu_bp * L_target    # kill-switch failure probability
    p_combined = p_ks_fail * p_thyA * p_dapA   # AND of all three failure modes

    P_escape = _p_escape_exact(N, G, p_combined)

    # Clamp at 1.0 to handle floating-point overshoot (P_escape <= 1
    # by definition; for very large N*G*p the value can round slightly
    # above 1.0 in IEEE-754, which is meaningless).
    P_escape = np.minimum(P_escape, 1.0)

    # -----------------------------------------------------------------
    # 3. Summary statistics
    # -----------------------------------------------------------------
    median = float(np.median(P_escape))
    ci_lo = float(np.quantile(P_escape, 0.025))
    ci_hi = float(np.quantile(P_escape, 0.975))
    frac_above = float(np.mean(P_escape > DE_MINIMIS_RISK_THRESHOLD))

    # Spearman rank correlation of each input vs P_escape
    # log-transform inputs and P_escape for the rank correlation
    log_samples = np.log10(samples)         # (N, 8)
    log_P = np.log10(np.maximum(P_escape, 1e-50))  # clamp at 1e-50 to avoid -inf
    rhos = []
    for j in range(8):
        rho, _ = spearmanr(log_samples[:, j], log_P)
        rhos.append((PARAM_NAMES[j], float(rho)))
    rhos.sort(key=lambda x: abs(x[1]), reverse=True)

    summary = {
        "median": median,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "frac_above_threshold": frac_above,
        "top_corr_param_1": rhos[0][0],
        "top_corr_param_2": rhos[1][0],
        "top_corr_param_3": rhos[2][0],
        "rho_1": rhos[0][1],
        "rho_2": rhos[1][1],
        "rho_3": rhos[2][1],
    }

    # -----------------------------------------------------------------
    # 4. Write CSV of all 10,000 samples
    # -----------------------------------------------------------------
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, "biosafety_lhs.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "mu_bp_per_gen",
            "p_HGT_thyA_per_div",
            "p_HGT_dapA_per_div",
            "V_reactor_L",
            "rho_max_per_L",
            "g_per_day",
            "T_days",
            "L_target_bp",
            "N_total_cells",
            "G_total_generations",
            "p_combined_per_div",
            "P_escape_30d",
        ])
        for i in range(LHS_N_SAMPLES):
            writer.writerow([
                f"{samples[i, 0]:.6e}",
                f"{samples[i, 1]:.6e}",
                f"{samples[i, 2]:.6e}",
                f"{samples[i, 3]:.6e}",
                f"{samples[i, 4]:.6e}",
                f"{samples[i, 5]:.6e}",
                f"{samples[i, 6]:.6e}",
                f"{samples[i, 7]:.6e}",
                f"{N[i]:.6e}",
                f"{G[i]:.6e}",
                f"{p_combined[i]:.6e}",
                f"{P_escape[i]:.6e}",
            ])

    # -----------------------------------------------------------------
    # 5. 3-panel figure
    # -----------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # Panel 1: CDF of P_escape (log-x axis)
    sorted_P = np.sort(P_escape)
    cdf = np.arange(1, LHS_N_SAMPLES + 1) / LHS_N_SAMPLES
    axes[0].plot(sorted_P, cdf, lw=1.5, color="C0")
    axes[0].axvline(DE_MINIMIS_RISK_THRESHOLD, color="orange", linestyle=":", lw=2,
                    label=f"10^-15 threshold")
    axes[0].axvline(median, color="C2", linestyle="--", lw=1.5,
                    label=f"median = {median:.2e}")
    axes[0].set_xscale("log")
    axes[0].set_xlabel("P_escape (30-day)")
    axes[0].set_ylabel("CDF")
    axes[0].set_title("(a) CDF of P_escape over 10,000 LHS samples")
    axes[0].legend(loc="lower right", fontsize=9)
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].set_xlim(1e-35, 2.0)

    # Panel 2: histogram of log10(P_escape)
    log_P_plot = np.log10(np.maximum(P_escape, 1e-50))
    axes[1].hist(log_P_plot, bins=60, color="C0", edgecolor="black", linewidth=0.4)
    axes[1].axvline(np.log10(DE_MINIMIS_RISK_THRESHOLD), color="orange",
                    linestyle=":", lw=2, label=f"10^-15 threshold")
    axes[1].axvline(np.log10(median), color="C2", linestyle="--", lw=1.5,
                    label=f"median = {median:.2e}")
    axes[1].set_xlabel("log10(P_escape)")
    axes[1].set_ylabel("count")
    axes[1].set_title("(b) Distribution of log10(P_escape)")
    axes[1].legend(loc="upper left", fontsize=9)
    axes[1].grid(True, alpha=0.3)

    # Panel 3: Spearman rank correlation heatmap
    rho_vec = np.array([r[1] for r in rhos[::-1]])  # restore input order
    labels_short = [
        "mu_bp", "p_thyA", "p_dapA", "V", "rho_max", "g", "T", "L_target"
    ]
    im = axes[2].imshow(
        rho_vec.reshape(1, -1),
        aspect="auto",
        cmap="RdBu_r",
        vmin=-1.0,
        vmax=1.0,
    )
    axes[2].set_xticks(range(8))
    axes[2].set_xticklabels(labels_short, rotation=45, ha="right", fontsize=9)
    axes[2].set_yticks([0])
    axes[2].set_yticklabels(["rho_spearman(log10(input), log10(P_esc))"])
    axes[2].set_title("(c) Rank correlation: input -> P_escape")
    for j in range(8):
        axes[2].text(j, 0, f"{rho_vec[j]:+.2f}", ha="center", va="center",
                     color="white" if abs(rho_vec[j]) > 0.5 else "black", fontsize=9)
    fig.colorbar(im, ax=axes[2], fraction=0.06, pad=0.08)

    fig.suptitle(
        f"Global LHS Uncertainty Quantification (Pillar 4c): N={LHS_N_SAMPLES} samples, "
        f"8 parameters, log-uniform",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()
    png_path = os.path.join(out_dir, "biosafety_lhs.png")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # -----------------------------------------------------------------
    # 6. Stdout summary
    # -----------------------------------------------------------------
    print("=" * 72)
    print("Pillar 4c: LHS Global Uncertainty Quantification")
    print("=" * 72)
    print(f"  N samples:               {LHS_N_SAMPLES}")
    print(f"  Random seed:             {LHS_RANDOM_SEED}")
    print(f"  Distribution:            {LHS_DISTRIBUTION}")
    print(f"  Number of parameters:    {len(PARAM_NAMES)}")
    print("-" * 72)
    print(f"  Median P_escape (30d):   {median:.3e}")
    print(f"  95% CI:                  [{ci_lo:.3e}, {ci_hi:.3e}]")
    print(f"  Fraction above 10^-15:   {frac_above:.4%}")
    print("-" * 72)
    print("  Top-3 Spearman rank correlations |log10(input), log10(P_esc)|:")
    for k, (name, rho) in enumerate(rhos[:3], 1):
        print(f"    {k}. {name:25s}  rho = {rho:+.3f}")
    print("-" * 72)
    print(f"  CSV:    {csv_path}")
    print(f"  Figure: {png_path}")
    print("=" * 72)

    return summary


if __name__ == "__main__":
    run_lhs()
