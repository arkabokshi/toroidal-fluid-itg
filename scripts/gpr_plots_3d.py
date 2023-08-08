"""
Produces a series of plots, given input data from a GPR run, showing:
    - The GPR predicted mean
    - The GPR prediction uncertainty

The input data is expected to be a csv file, with the column headers:
run #  |  eta_g  |  epsilon_n  |  shear  |  FWHM

The plots produced by this script will be writted to a directory named:
./images
"""
import csv
from typing import List
from math import pi
from pathlib import Path
from numpy import array, linspace, where, zeros
import matplotlib.pyplot as plt
from inference.gp import GpRegressor

eta_g: List[float] = []
epsilon_n: List[float] = []
shear: List[float] = []
fwhm: List[float] = []

with open(
    "../Data/GPO_200_3D/GPO_results.csv",
    newline="",
    encoding="utf-8",
) as csvfile:
    results = csv.reader(csvfile, delimiter=",")
    for i, row in enumerate(results):
        if i == 0:
            # Don't read the headers.
            continue

        eta_g.append(float(row[1]))
        epsilon_n.append(float(row[2]))
        shear.append(float(row[3]))
        fwhm.append(float(row[4]))

# Convert lists to numpy arrays
eta_g = array(eta_g)
epsilon_n = array(epsilon_n)
shear = array(shear)
fwhm = array(fwhm)

# Convert FWHM data to the objective function value for clarity.
# Objective function was simply searching for the negative of the FWHM.
objective_fn = fwhm * -1

# GPR uses spatial data of shape [n_points, n_dimensions]
x_data = array([eta_g, epsilon_n, shear]).T

# train the GP using the data
gpr = GpRegressor(x=x_data, y=objective_fn, cross_val=False)

# create axes for the two variables on which to evaluate the GP
AXIS_SIZE = 128
eta_axis = linspace(eta_g.min(), eta_g.max(), AXIS_SIZE)
eps_axis = linspace(epsilon_n.min(), epsilon_n.max(), AXIS_SIZE)
shear_axis = linspace(shear.min(), shear.max(), AXIS_SIZE)
mu_grid = zeros([AXIS_SIZE, AXIS_SIZE, AXIS_SIZE])
sg_grid = zeros([AXIS_SIZE, AXIS_SIZE, AXIS_SIZE])

# loop over all the points in the grid
for i in range(AXIS_SIZE):
    for j in range(AXIS_SIZE):
        for k in range(AXIS_SIZE):
            mu, sigma = gpr([eta_axis[i], eps_axis[j], shear_axis[k]])
            mu_grid[i, j, k] = mu
            sg_grid[i, j, k] = sigma

# Ensure the images directory exists
Path("./images").mkdir(parents=True, exist_ok=True)

for i in range(AXIS_SIZE):
    if i == 0:
        filter_arr = shear <= shear_axis[i + 1]
    elif i == AXIS_SIZE - 1:
        filter_arr = shear >= shear_axis[i - 1]
    else:
        filter_arr = where((shear_axis[i - 1] < shear) & (shear < shear_axis[i + 1]))

    fig = plt.figure(figsize=(12, 5))
    fig.suptitle(f"Shear ~{shear_axis[i]:4.2f}")
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)

    # plot the gp prediction mean
    levels_1 = linspace(-2 * pi, 0, AXIS_SIZE)
    cs1 = ax1.contourf(
        eta_axis, eps_axis, mu_grid[:, :, i].T, levels=levels_1, extend="min"
    )
    cbar1 = plt.colorbar(cs1, ax=ax1)
    ax1.plot(
        eta_g[filter_arr], epsilon_n[filter_arr], marker="x", color="red", ls="none"
    )
    ax1.set_xlabel(r"$\eta_g$", fontsize=15)
    ax1.set_ylabel(r"$\epsilon_n$", fontsize=15)
    ax1.set_title("GPR prediction mean")

    # plot the gp prediction standard-deviation
    max_val = sg_grid[:, :, i].flatten().max()
    levels_2 = linspace(0, max_val, AXIS_SIZE)
    cs2 = ax2.contourf(eta_axis, eps_axis, sg_grid[:, :, i].T, levels=levels_2)
    cbar2 = plt.colorbar(cs2, ax=ax2)
    ax2.plot(
        eta_g[filter_arr], epsilon_n[filter_arr], marker="x", color="red", ls="none"
    )
    ax2.set_xlabel(r"$\eta_g$", fontsize=15)
    ax2.set_ylabel(r"$\epsilon_n$", fontsize=15)
    ax2.set_title("GPR prediction uncertainty")

    plt.tight_layout()
    plt.savefig(f"./images/{i}.png", format="png", dpi=250)
    plt.close()
