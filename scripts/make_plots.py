#!/usr/bin/env python3
"""Regenerate every figure in docs/images from the simulator's CSV output.

Usage (from the repository root, after building):
    python3 scripts/make_plots.py --build build

Runs the study and CLI executables, then plots:
    convergence.png   phase-space error vs step size (observed order of accuracy)
    energy.png        relative energy error over 1000 orbits at equal cost
    scaling.png       direct-summation force cost vs number of bodies
    orbits.png        figure-eight and inner solar system trajectories
    figure8.gif       animation of the figure-eight three-body orbit
"""

import argparse
import csv
import math
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
IMAGES = ROOT / "docs" / "images"

# Categorical palette, fixed order (validated for colour-vision deficiency on the
# adjacent pairs used by line charts). Marker shapes and direct labels carry
# identity too, so no series relies on colour alone.
SERIES = {
    "euler": ("#2a78d6", "o", "Explicit Euler"),
    "symplectic-euler": ("#eb6834", "s", "Symplectic Euler"),
    "verlet": ("#1baf7a", "^", "Velocity Verlet"),
    "rk4": ("#eda100", "D", "RK4"),
    "yoshida4": ("#e87ba4", "v", "Yoshida 4"),
}
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.edgecolor": GRID,
    "axes.labelcolor": INK_2,
    "axes.titlecolor": INK,
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "xtick.color": INK_2,
    "ytick.color": INK_2,
    "font.size": 10.5,
    "legend.frameon": False,
    "lines.linewidth": 2.0,
    "lines.markersize": 6,
    "savefig.dpi": 160,
    "savefig.bbox": "tight",
})


def run(cmd):
    subprocess.run([str(c) for c in cmd], check=True, cwd=ROOT, stdout=subprocess.DEVNULL)


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def slope(xs, ys):
    lx = [math.log(x) for x in xs]
    ly = [math.log(y) for y in ys]
    n = len(lx)
    mx, my = sum(lx) / n, sum(ly) / n
    return sum((a - mx) * (b - my) for a, b in zip(lx, ly)) / sum((a - mx) ** 2 for a in lx)


def plot_convergence():
    rows = defaultdict(list)
    for r in read_csv(RESULTS / "convergence.csv"):
        rows[r["integrator"]].append(r)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for name, (color, marker, label) in SERIES.items():
        evals = [int(r["force_evals"]) for r in rows[name]]
        err = [float(r["phase_space_error"]) for r in rows[name]]
        dts = [float(r["dt"]) for r in rows[name]]
        # Fit the asymptotic order on points well above the round-off floor.
        fit = [(d, e) for d, e in zip(dts, err) if e > 1e-10][-4:]
        p = slope([d for d, _ in fit], [e for _, e in fit])
        ax.loglog(evals, err, color=color, marker=marker, markeredgecolor=SURFACE,
                  markeredgewidth=1.2, label=f"{label}  (order {p:.2f})")
        ax.annotate(label, (evals[0], err[0]), xytext=(-6, 0), textcoords="offset points",
                    ha="right", va="center", color=INK_2, fontsize=9)

    ax.axhspan(1e-17, 3e-12, color=GRID, alpha=0.45, lw=0)
    ax.text(15, 3e-13, "floating-point round-off floor", color=INK_2, fontsize=9)
    ax.set_ylim(1e-14, 10)
    ax.set_xlim(left=12)
    ax.set_xlabel("force evaluations (cost)")
    ax.set_ylabel("phase-space error |Δr| + |Δv|")
    ax.set_title("Error vs cost on an eccentric Kepler orbit (e = 0.5)")
    ax.legend(loc="lower left", fontsize=9)
    fig.savefig(IMAGES / "convergence.png")
    plt.close(fig)


def plot_energy():
    rows = defaultdict(list)
    for r in read_csv(RESULTS / "longterm.csv"):
        rows[r["integrator"]].append(r)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for name, (color, marker, label) in SERIES.items():
        periods = [int(r["periods"]) for r in rows[name]]
        err = [float(r["max_rel_energy_error"]) for r in rows[name]]
        ax.loglog(periods, err, color=color, label=label)
        ax.annotate(label, (periods[-1], err[-1]), xytext=(6, 0), textcoords="offset points",
                    va="center", color=INK_2, fontsize=9)
    ax.set_xlim(1, 4000)
    ax.set_xlabel("orbits completed")
    ax.set_ylabel("max relative energy error per orbit |ΔE / E₀|")
    ax.set_title("Energy error over 1000 orbits at equal cost (600 force evals / orbit)")
    fig.savefig(IMAGES / "energy.png")
    plt.close(fig)


def plot_scaling():
    rows = read_csv(RESULTS / "scaling.csv")
    n = [int(r["n"]) for r in rows]
    t = [float(r["seconds_per_eval"]) for r in rows]
    p = slope(n, t)

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.loglog(n, t, color=SERIES["euler"][0], marker="o", markeredgecolor=SURFACE,
              markeredgewidth=1.2, label="measured")
    mid = len(n) // 2  # anchor the reference mid-range, away from cache effects at the ends
    ref = [t[mid] * (x / n[mid]) ** 2 for x in n]
    ax.loglog(n, ref, color=INK_2, lw=1.2, ls="--", label="ideal O(N²)")
    ax.annotate(f"fitted exponent {p:.2f}", (n[-1], t[-1]), xytext=(-10, 6),
                textcoords="offset points", ha="right", color=INK_2, fontsize=9)
    ax.set_xlabel("number of bodies N")
    ax.set_ylabel("seconds per force evaluation")
    ax.set_title("Direct-summation cost (single core)")
    ax.legend(loc="upper left", fontsize=9)
    fig.savefig(IMAGES / "scaling.png")
    plt.close(fig)


def trajectory(path):
    tracks = defaultdict(lambda: ([], []))
    for r in read_csv(path):
        xs, ys = tracks[int(r["body"])]
        xs.append(float(r["x"]))
        ys.append(float(r["y"]))
    return tracks


def plot_orbits(build):
    cli = build / "nbody"
    fig8_csv = RESULTS / "figure8_traj.csv"
    solar_csv = RESULTS / "solar_traj.csv"
    run([cli, "--scenario", "figure8", "--integrator", "yoshida4", "--dt", "0.001",
         "--t-end", "6.3259", "--output-every", "10", "--trajectory", fig8_csv,
         "--diagnostics", "", "--quiet"])
    run([cli, "--scenario", "solar", "--integrator", "yoshida4", "--dt", "0.0005",
         "--t-end", "2", "--output-every", "4", "--trajectory", solar_csv,
         "--diagnostics", "", "--quiet"])

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.6))
    tracks = trajectory(fig8_csv)
    # All three bodies follow the same curve, so draw it once and colour the bodies.
    a1.plot(*tracks[0], color=INK_2, lw=1.4)
    for i, (xs, ys) in tracks.items():
        color = list(SERIES.values())[i][0]
        a1.plot(xs[0], ys[0], "o", color=color, ms=11, markeredgecolor=SURFACE, markeredgewidth=2)
    a1.set_title("Figure-eight: three bodies, one shared path")
    a1.set_aspect("equal")
    a1.set_xlabel("x")
    a1.set_ylabel("y")

    tracks = trajectory(solar_csv)
    labels = ["Sun", "Mercury", "Venus", "Earth", "Mars"]
    for i in range(5):
        xs, ys = tracks[i]
        color = INK if i == 0 else list(SERIES.values())[i - 1][0]
        a2.plot(xs, ys, color=color, lw=1.6)
        a2.plot(xs[-1], ys[-1], "o", color=color, ms=7 if i else 10,
                markeredgecolor=SURFACE, markeredgewidth=1.5)
        if i:
            a2.annotate(labels[i], (xs[-1], ys[-1]), xytext=(6, 4), textcoords="offset points",
                        color=INK_2, fontsize=9)
    a2.set_title("Inner solar system, 2 years")
    a2.set_aspect("equal")
    a2.set_xlim(-1.8, 1.8)
    a2.set_ylim(-1.8, 1.8)
    a2.set_xlabel("x (AU)")
    a2.set_ylabel("y (AU)")
    fig.tight_layout()
    fig.savefig(IMAGES / "orbits.png")
    plt.close(fig)
    return fig8_csv


def animate_figure8(fig8_csv):
    from matplotlib.animation import FuncAnimation, PillowWriter

    tracks = trajectory(fig8_csv)
    colors = [v[0] for v in list(SERIES.values())[:3]]
    frames = len(tracks[0][0])
    step = 3

    fig, ax = plt.subplots(figsize=(5.2, 2.8))
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-0.5, 0.5)
    ax.set_aspect("equal")
    ax.axis("off")
    for i in range(3):
        ax.plot(*tracks[i], color=GRID, lw=1.2, zorder=1)
    trails = [ax.plot([], [], color=c, lw=2, zorder=2)[0] for c in colors]
    dots = [ax.plot([], [], "o", color=c, ms=10, markeredgecolor=SURFACE,
                    markeredgewidth=2, zorder=3)[0] for c in colors]

    def update(f):
        k = f * step
        for i in range(3):
            xs, ys = tracks[i]
            lo = max(0, k - 40)
            trails[i].set_data(xs[lo:k + 1], ys[lo:k + 1])
            dots[i].set_data([xs[k]], [ys[k]])
        return trails + dots

    anim = FuncAnimation(fig, update, frames=frames // step, blit=True)
    anim.save(IMAGES / "figure8.gif", writer=PillowWriter(fps=30), dpi=90)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--build", default="build", help="CMake build directory")
    args = parser.parse_args()
    build = (ROOT / args.build).resolve()

    RESULTS.mkdir(exist_ok=True)
    IMAGES.mkdir(parents=True, exist_ok=True)
    study = build / "nbody_study"
    for mode in ("convergence", "longterm", "scaling"):
        run([study, mode, RESULTS / f"{mode}.csv"])

    plot_convergence()
    plot_energy()
    plot_scaling()
    fig8_csv = plot_orbits(build)
    animate_figure8(fig8_csv)
    print(f"figures written to {IMAGES.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
