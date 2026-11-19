"""Render the published OGLE fit and photometric errors as a static figure."""

import json
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
data = json.loads((ROOT / "data/ogle-2023-blg-0100-fit.json").read_text(encoding="utf-8"))
observations = [row for row in data["observations"] if abs(row["day"]) <= 90]
curve = data["curve"]
day = [row["day"] for row in observations]
residual = [row["standardised_residual"] for row in observations]

with plt.rc_context({"svg.fonttype": "none", "svg.hashsalt": "ogle-0100-evidence", "font.size": 12}):
    figure, (flux_ax, residual_ax) = plt.subplots(
        2, 1, figsize=(8.4, 6.2), sharex=True, height_ratios=(2, 1), layout="constrained"
    )
    figure.suptitle("OGLE-2023-BLG-0100 · I-band evidence", x=0.01, ha="left", fontweight="bold")
    flux_ax.set_title(f"{len(observations)} measurements within 90 days of peak", loc="left", fontsize=11)
    flux_ax.errorbar(
        day, [row["relative_flux"] for row in observations],
        yerr=[row["flux_error"] for row in observations],
        fmt="o", ms=2.8, color="#223c50", ecolor="#8395a1", lw=0.75, label="OGLE-IV ±1σ",
    )
    flux_ax.plot(
        [row["day"] for row in curve], [row["relative_flux"] for row in curve],
        color="#b14930", lw=1.7, label="Point lens + constant blend",
    )
    flux_ax.set_ylabel("Relative flux")
    flux_ax.legend(loc="upper right", frameon=False, fontsize=9)
    residual_ax.axhspan(-3, 3, color="#eef1f3")
    residual_ax.axhline(0, color="#68747b", lw=0.8)
    residual_ax.scatter(
        day, residual, s=14,
        color=["#b14930" if abs(value) > 3 else "#223c50" for value in residual],
    )
    residual_ax.set(xlabel="Days from fitted peak", ylabel="Residual / σ", xlim=(-90, 90))
    for axis in (flux_ax, residual_ax):
        axis.spines[["top", "right"]].set_visible(False)
    figure.savefig(ROOT / "assets/ogle-lightcurve-evidence.svg", metadata={"Date": None})
    plt.close(figure)
