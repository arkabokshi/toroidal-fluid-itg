"""
Produces heat map style figures to visualise the GPO output data.
Uses interpolation to fill the gaps between the data.
"""
import csv
from typing import List
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

eta_g: List[float] = []
epsilon_n: List[float] = []
fwhm: List[float] = []
# DATA_PATH = "../Data/GPO_50_Iterations_Run"
# DATA_PATH = "../Data/GPO_First_Run"
DATA_PATH = "../Data/GPO_Expected_Improvement"
with open(f"{DATA_PATH}/GPO_results.csv", newline="", encoding="utf-8") as csvfile:
    results = csv.reader(csvfile, delimiter=",")
    for i, row in enumerate(results):
        if i == 0 or i == 10:
            # Don't read the headers.
            continue

        eta_g.append(float(row[1]))
        epsilon_n.append(float(row[2]))
        fwhm.append(float(row[3]))

# Convert lists to numpy arrays
eta_g = np.array(eta_g)
epsilon_n = np.array(epsilon_n)
fwhm = np.array(fwhm)

# Convert FWHM data to the objective function value for clarity
objective_fn = fwhm * -1

# Create x-y points to be used in heatmap
xi = np.linspace(eta_g.min(), eta_g.max(), 1000)
yi = np.linspace(epsilon_n.min(), epsilon_n.max(), 1000)

# Interpolate for plotting
zi = griddata(
    (eta_g, epsilon_n), objective_fn, (xi[None, :], yi[:, None]), method="cubic"
)

# Remove any interpolated values that go beyond my chosen minimum (-2 pi),
# or those that go above zero (non-physical).
Z_MAX = 0
Z_MIN = -2 * np.pi
zi[(zi < Z_MIN) | (zi > Z_MAX)] = None

fig = plt.figure(figsize=(10, 6), dpi=300)

# Create the contour plot
levels = np.linspace(Z_MIN, Z_MAX, 100)
CS = plt.contourf(xi, yi, zi, cmap="gnuplot", levels=levels)

ticks = np.linspace(Z_MIN, Z_MAX, 3)
c_bar = plt.colorbar()
c_bar.set_label("Objective Function Value [radians]", labelpad=10)
c_bar.set_ticks(ticks)
c_bar.ax.set_yticklabels([r"$-2\pi$", r"$-\pi$", "0"])

plt.xlim(1, 5)
plt.ylim(0.01, 0.1)
plt.xlabel(
    r"$\eta_g$ [arb. units]",
    fontsize=15.0,
)
plt.ylabel(
    r"$\epsilon_n$ [arb. units]",
    fontsize=15.0,
)
plt.title("30 Iterations - ExpectedImprovement - Without Anomalous Result")
plt.scatter(eta_g, epsilon_n, c="black", marker="x")

fig.savefig(f"{DATA_PATH}/GPO_analysis.png", format="png", dpi=250, bbox_inches="tight")
