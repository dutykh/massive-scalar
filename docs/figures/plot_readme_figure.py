#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: stable spectral roots of the massive scalar field, read from the reports in data/.

From this folder,
    python3 plot_readme_figure.py
writes massive-scalar-roots.png: for the multipoles l = 0, 1, 2 and five masses mu = M m, the complex roots
Omega = M omega of the convergence reports (three resolutions, tolerance 1e-4) that lie off the real axis,
Im Omega < -1e-3. The style is the stamped copy manuscript_style.py beside this script.
"""
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "data"
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

MASSES = (("0_0", "0"), ("0_2", "0.2"), ("0_8", "0.8"), ("1_45", "1.45"), ("2_2", "2.2"))
MULTIPOLES = (0, 1, 2)
OFF_AXIS = 1e-3   # roots closer to the real axis cluster at Omega = mu and are left out


def general_roots(folder):
    """Complex roots listed under "General QNMs" in the report of one case."""
    report = next(folder.glob("*.txt"))
    roots, inside = [], False
    for line in report.read_text().splitlines():
        if re.match(r"^(General|Purely Real|Purely Imaginary) QNMs", line):
            inside = line.startswith("General")
            continue
        parts = line.split()
        if inside and len(parts) == 2:
            try:
                roots.append(complex(float(parts[0]), float(parts[1])))
            except ValueError:
                pass
    return np.array(roots)


hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 265
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig, axes = plt.subplots(1, 3, figsize=(hs.WIDTH, 2.5), sharey=True)
for axis, ell in zip(axes, MULTIPOLES):
    for k, (tag, label) in enumerate(MASSES):
        roots = general_roots(DATA / f"mu_{tag}_L_{ell}")
        roots = roots[roots.imag < -OFF_AXIS]
        style = hs.series(k)
        axis.plot(roots.real, roots.imag, linestyle="none", marker=style["marker"], markersize=2.8,
                  markerfacecolor="white", markeredgecolor=style["color"], markeredgewidth=0.8,
                  label=rf"$\mu = {label}$")
    axis.set_xlim(0.0, 2.5)
    axis.set_ylim(-3.2, 0.1)
    axis.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
    axis.set_xlabel(r"$\mathrm{Re}\,\Omega$")
    hs.panel(axis, rf"({'abc'[ell]}) $\ell = {ell}$")
axes[0].set_ylabel(r"$\mathrm{Im}\,\Omega$")

hs.legend_outside(fig, where="upper", ncol=5)
out = HERE / "massive-scalar-roots.png"
hs.save(fig, out)
print("written", out)
