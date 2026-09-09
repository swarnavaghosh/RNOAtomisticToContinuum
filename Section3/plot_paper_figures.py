#!/usr/bin/env python3

from pathlib import Path
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "paper_figure_data"
OUTPUT_DIR = HERE / "paper_figures_output"
OUTPUT_DIR.mkdir(exist_ok=True)

# Default: use LaTeX font rendering.
USE_TEX = os.environ.get("PAPER_USE_TEX", "1").strip() not in {"0", "false", "False"}

PAPER_STYLE = {
    "figure_size": (3.55, 2.85),
    "screen_dpi": 120,
    "save_dpi": 600,
    "base_font": 9.0,
    "axis_label_font": 9.0,
    "tick_font": 8.0,
    "legend_font": 8.0,
    "inset_label_font": 6.5,
    "inset_tick_font": 6.0,
    "line_width": 1.45,
    "reference_line_width": 1.10,
    "marker_size": 4.2,
    "axes_line_width": 0.75,
    "tick_width": 0.75,
    "left": 0.16,
    "right": 0.98,
    "top": 0.98,
    "bottom": 0.18,
}

matplotlib.rcParams.update({
    "text.usetex": USE_TEX,
    "text.latex.preamble": r"\usepackage{amsmath}",
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif"],
    "font.size": PAPER_STYLE["base_font"],
    "axes.labelsize": PAPER_STYLE["axis_label_font"],
    "xtick.labelsize": PAPER_STYLE["tick_font"],
    "ytick.labelsize": PAPER_STYLE["tick_font"],
    "legend.fontsize": PAPER_STYLE["legend_font"],
    "figure.dpi": PAPER_STYLE["screen_dpi"],
    "savefig.dpi": PAPER_STYLE["save_dpi"],
    "axes.linewidth": PAPER_STYLE["axes_line_width"],
    "xtick.major.width": PAPER_STYLE["tick_width"],
    "ytick.major.width": PAPER_STYLE["tick_width"],
    "xtick.direction": "out",
    "ytick.direction": "out",
    "axes.grid": False,
    "legend.frameon": False,
    "lines.solid_capstyle": "round",
    "lines.dash_capstyle": "butt",
})

def load_csv(filename):
    path = DATA_DIR / filename
    if not path.exists():
        raise FileNotFoundError(f"Missing required data file: {path}")
    return np.genfromtxt(path, delimiter=",", names=True, encoding="utf-8")

def save_figure(fig, stem):
    fig.savefig(OUTPUT_DIR / f"{stem}.pdf", bbox_inches="tight", dpi=PAPER_STYLE["save_dpi"])
    fig.savefig(OUTPUT_DIR / f"{stem}.png", bbox_inches="tight", dpi=PAPER_STYLE["save_dpi"])

def _finish_single_axes(fig):
    fig.subplots_adjust(
        left=PAPER_STYLE["left"],
        right=PAPER_STYLE["right"],
        top=PAPER_STYLE["top"],
        bottom=PAPER_STYLE["bottom"],
    )

def make_bulk_modulus():
    d = load_csv("bulk_modulus.csv")
    fig, ax = plt.subplots(figsize=PAPER_STYLE["figure_size"])
    ax.plot(
        d["cell_size_A"], d["bulk_modulus_GPa"], "-o",
        color="red", linewidth=PAPER_STYLE["line_width"],
        markersize=PAPER_STYLE["marker_size"],
        markerfacecolor="red", markeredgecolor="black", markeredgewidth=0.45,
        label=r"MD simulation"
    )
    ax.axhline(2.20, color="black", linestyle="--",
               linewidth=PAPER_STYLE["reference_line_width"], label=r"Converged value")
    ax.axvline(120, color="black", linestyle="--",
               linewidth=PAPER_STYLE["reference_line_width"])
    ax.set_xlim(0, 150)
    ax.set_ylim(1.90, 2.22)
    ax.set_xticks(np.arange(0, 141, 20))
    ax.set_yticks(np.arange(1.90, 2.201, 0.05))
    ax.set_xlabel(r"Cell size ($\mathrm{\AA}$)")
    ax.set_ylabel(r"Bulk modulus (GPa)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.58, 0.05),
              handlelength=2.5, borderaxespad=0.0)
    _finish_single_axes(fig)
    save_figure(fig, "bulk_modulus")
    plt.close(fig)

def make_stress_strain():
    d = load_csv("volumetric_stress_strain.csv")
    fig, ax = plt.subplots(figsize=PAPER_STYLE["figure_size"])
    ax.plot(d["volumetric_strain"], d["stress_300K_GPa"],
            color="red", linewidth=PAPER_STYLE["line_width"], label=r"$300\,\mathrm{K}$")
    ax.plot(d["volumetric_strain"], d["stress_400K_GPa"],
            color="blue", linewidth=PAPER_STYLE["line_width"], label=r"$400\,\mathrm{K}$")
    ax.plot(d["volumetric_strain"], d["stress_500K_GPa"],
            color="black", linewidth=PAPER_STYLE["line_width"], label=r"$500\,\mathrm{K}$")
    ax.set_xlim(-0.10, 0.11)
    ax.set_ylim(-0.36, 0.46)
    ax.set_xticks(np.arange(-0.10, 0.1001, 0.025))
    ax.set_yticks(np.arange(-0.30, 0.401, 0.10))
    ax.set_xlabel(r"Volumetric strain")
    ax.set_ylabel(r"Volumetric stress (GPa)")
    ax.legend(loc="upper left", handlelength=2.3, borderaxespad=0.4)
    _finish_single_axes(fig)
    save_figure(fig, "stress_strain_plot")
    plt.close(fig)

def make_internal_variables():
    d = load_csv("internal_variables.csv")
    fig, ax = plt.subplots(figsize=PAPER_STYLE["figure_size"])
    ax.plot(d["internal_variables"], d["relative_error_300K"], "-o",
            color="red", linewidth=PAPER_STYLE["line_width"],
            markersize=PAPER_STYLE["marker_size"], label=r"$300\,\mathrm{K}$")
    ax.plot(d["internal_variables"], d["relative_error_400K"], "-s",
            color="blue", linewidth=PAPER_STYLE["line_width"],
            markersize=PAPER_STYLE["marker_size"], label=r"$400\,\mathrm{K}$")
    ax.plot(d["internal_variables"], d["relative_error_500K"], "-^",
            color="black", linewidth=PAPER_STYLE["line_width"],
            markersize=PAPER_STYLE["marker_size"], label=r"$500\,\mathrm{K}$")
    ax.set_xlim(0, 50)
    ax.set_ylim(0, 0.63)
    ax.set_xticks(np.arange(0, 51, 10))
    ax.set_yticks(np.arange(0, 0.61, 0.10))
    ax.set_xlabel(r"Internal variables")
    ax.set_ylabel(r"Relative error")
    ax.legend(loc="upper right", handlelength=2.3, borderaxespad=0.35)
    _finish_single_axes(fig)
    save_figure(fig, "internal_variables")
    plt.close(fig)

def make_relative_error_epoch():
    d300 = load_csv("figure4b_300K.csv")
    d400 = load_csv("figure4b_400K.csv")
    d500 = load_csv("figure4b_500K.csv")
    fig, ax = plt.subplots(figsize=PAPER_STYLE["figure_size"])
    ax.plot(d300["epoch"], d300["train_error"], color="red", linestyle="-",
            linewidth=PAPER_STYLE["line_width"], label=r"$300\,\mathrm{K}$ Train")
    ax.plot(d300["epoch"], d300["test_error"], color="red", linestyle="--",
            linewidth=PAPER_STYLE["line_width"], label=r"$300\,\mathrm{K}$ Test")
    ax.plot(d400["epoch"], d400["train_error"], color="blue", linestyle="-",
            linewidth=PAPER_STYLE["line_width"], label=r"$400\,\mathrm{K}$ Train")
    ax.plot(d400["epoch"], d400["test_error"], color="blue", linestyle="--",
            linewidth=PAPER_STYLE["line_width"], label=r"$400\,\mathrm{K}$ Test")
    ax.plot(d500["epoch"], d500["train_error"], color="black", linestyle="-",
            linewidth=PAPER_STYLE["line_width"], label=r"$500\,\mathrm{K}$ Train")
    ax.plot(d500["epoch"], d500["test_error"], color="black", linestyle="--",
            linewidth=PAPER_STYLE["line_width"], label=r"$500\,\mathrm{K}$ Test")
    ax.set_xlim(-5, 605)
    ax.set_ylim(-0.0025, 0.0715)
    ax.set_xticks(np.arange(0, 601, 100))
    ax.set_yticks(np.arange(0.00, 0.071, 0.01))
    ax.set_xlabel(r"Epoch")
    ax.set_ylabel(r"Relative error")
    ax.legend(loc="upper right", ncol=2, columnspacing=1.05,
              handlelength=2.7, handletextpad=0.55,
              borderaxespad=0.45, fontsize=7.2)

    inset = ax.inset_axes([0.47, 0.36, 0.49, 0.49])
    for d, color, ls in [
        (d300, "red", "-"), (d300, "red", "--"),
        (d400, "blue", "-"), (d400, "blue", "--"),
        (d500, "black", "-"), (d500, "black", "--"),
    ]:
        yfield = "train_error" if ls == "-" else "test_error"
        mask = (d["epoch"] >= 0) & (d["epoch"] <= 10)
        inset.plot(d["epoch"][mask], d[yfield][mask],
                   color=color, linestyle=ls, linewidth=0.95)

    inset.set_xlim(0, 10)
    inset.set_ylim(0, 0.0715)
    inset.set_xticks(np.arange(0, 11, 2))
    inset.set_yticks(np.arange(0.00, 0.071, 0.01))
    inset.set_xlabel(r"Epoch", fontsize=PAPER_STYLE["inset_label_font"])
    inset.set_ylabel(r"Relative error", fontsize=PAPER_STYLE["inset_label_font"])
    inset.tick_params(axis="both", labelsize=PAPER_STYLE["inset_tick_font"],
                      width=0.55, length=2.3)
    for spine in inset.spines.values():
        spine.set_linewidth(0.60)

    _finish_single_axes(fig)
    save_figure(fig, "relative_error_epoch")
    plt.close(fig)

if __name__ == "__main__":
    make_bulk_modulus()
    make_stress_strain()
    make_internal_variables()
    make_relative_error_epoch()
    print("Done.")
