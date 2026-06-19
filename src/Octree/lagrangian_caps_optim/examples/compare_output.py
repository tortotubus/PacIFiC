#!/usr/bin/env python3
"""
Compare two one_cap output.txt files.

Columns (0-indexed):
  0  i              iteration
  1  t              time
  2  fluid_visc_stress
  3  ppres          particle pressure
  4  pmu            shear stress          <- focus
  5  pN1            first normal stress diff
  6  pN2            second normal stress diff
  7  avg_td         Taylor deformation
  8  avg_inclin     inclination angle
  9  avg_TDmaxmin
 10  avg_TDang
 11  avg_rs.x
 12  avg_rs.y
 13  avg_rs.z
 14  avg_angvel     angular velocity (z)
 15  avg_area
 16  avg_vol

Usage:
  python3 compare_output.py [file1] [file2] [label1] [label2]

Defaults:
  file1  = one_cap_data/output.txt      (optim/multigrid)
  file2  = noncubic_one_cap/output.txt  (reference)
"""

import sys
import os
import numpy as np
import matplotlib
if not os.environ.get("DISPLAY") and not os.environ.get("WAYLAND_DISPLAY"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paths and labels
# ---------------------------------------------------------------------------
file1  = sys.argv[1] if len(sys.argv) > 1 else "one_cap_data/output.txt"
file2  = sys.argv[2] if len(sys.argv) > 2 else "noncubic_one_cap/output.txt"
label1 = sys.argv[3] if len(sys.argv) > 3 else "optim/multigrid"
label2 = sys.argv[4] if len(sys.argv) > 4 else "noncubic reference"

def load(path):
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data

d1 = load(file1)
d2 = load(file2)

# ---------------------------------------------------------------------------
# Column indices
# ---------------------------------------------------------------------------
T        = 1
SHEAR    = 4   # pmu  — shear stress (5th column, 1-indexed)
N1       = 5   # pN1
N2       = 6   # pN2
TD       = 7   # Taylor deformation
INCLIN   = 8   # inclination angle
ANGVEL   = 14  # angular velocity z

# ---------------------------------------------------------------------------
# Figure layout: shear stress main panel + supporting diagnostics
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
fig.suptitle("one_cap comparison", fontsize=13)

def plot_pair(ax, col, ylabel, title, yscale="linear", skip=0):
    ax.plot(d1[skip:, T], d1[skip:, col], label=label1, lw=1.5)
    ax.plot(d2[skip:, T], d2[skip:, col], label=label2, lw=1.5, ls="--")
    ax.set_xlabel("t")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_yscale(yscale)
    ax.legend(fontsize=8)
    ax.grid(True, lw=0.4)

# Main panel: shear stress — larger
ax_main = axes[0, 0]
ax_main.plot(d1[:, T], d1[:, SHEAR], label=label1, lw=2)
ax_main.plot(d2[:, T], d2[:, SHEAR], label=label2, lw=2, ls="--")
ax_main.set_xlabel("t", fontsize=11)
ax_main.set_ylabel(r"$\Sigma_{xy}$ (shear stress)", fontsize=11)
ax_main.set_title("Shear stress vs time", fontsize=11)
ax_main.legend()
ax_main.grid(True, lw=0.4)

plot_pair(axes[0, 1], N1,     r"$N_1$",              "1st normal stress diff")
plot_pair(axes[0, 2], N2,     r"$N_2$",              "2nd normal stress diff")
plot_pair(axes[1, 0], TD,     "D (Taylor factor)",   "Taylor deformation")
plot_pair(axes[1, 1], INCLIN, r"$\theta$ (rad)",     "Inclination angle",  skip=1)
plot_pair(axes[1, 2], ANGVEL, r"$\Omega_z$",         "Angular velocity z", skip=1)

plt.tight_layout()
out = "comparison.png"
plt.savefig(out, dpi=150)
print(f"Saved {out}")
try:
    plt.show()
except Exception:
    pass
