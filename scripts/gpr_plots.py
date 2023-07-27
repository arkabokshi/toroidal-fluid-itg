"""
A script provided by Chris Bowman.
Produces three plots, given input data from a GPR run, showing:
    - The GPR predicted mean
    - The GPR prediction uncertainty
    - The "Leave-one-out" cross-validation
"""
import csv
from typing import List
from numpy import array, linspace, zeros
import matplotlib.pyplot as plt
from inference.gp import GpRegressor

eta_g: List[float] = []
epsilon_n: List[float] = []
fwhm: List[float] = []

with open(
    "../Data/GPO_50_Iterations_Run/GPO_results_appended.csv",
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
        fwhm.append(float(row[3]))

# Convert lists to numpy arrays
eta_g = array(eta_g)
epsilon_n = array(epsilon_n)
fwhm = array(fwhm)

# Convert FWHM data to the objective function value for clarity.
# Objective function was simply searching for the negative of the FWHM.
objective_fn = fwhm * -1

# GPR uses spatial data of shape [n_points, n_dimensions]
x_data = array([eta_g, epsilon_n]).T

# train the GP using the data
gpr = GpRegressor(x=x_data, y=objective_fn, cross_val=False)

# create axes for the two variables on which to evaluate the GP
n, m = 128, 128
eta_axis = linspace(eta_g.min(), eta_g.max(), n)
eps_axis = linspace(epsilon_n.min(), epsilon_n.max(), m)
mu_grid = zeros([n, m])
sg_grid = zeros([n, m])

# loop over all the points in the grid
for i in range(n):
    for j in range(m):
        mu, sigma = gpr([eta_axis[i], eps_axis[j]])
        mu_grid[i, j] = mu
        sg_grid[i, j] = sigma


fig = plt.figure(figsize=(12, 5))
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

# plot the gp prediction mean
cs1 = ax1.contourf(eta_axis, eps_axis, mu_grid.T, 64, vmin=-2.0)
cbar1 = plt.colorbar(cs1, ax=ax1)
ax1.plot(eta_g, epsilon_n, marker="x", color="red", ls="none")
ax1.set_xlabel(r"$\eta_g$", fontsize=15)
ax1.set_ylabel(r"$\epsilon_n$", fontsize=15)
ax1.set_title("GPR prediction mean")

# plot the gp prediction standard-deviation
cs2 = ax2.contourf(eta_axis, eps_axis, sg_grid.T, 64)
cbar2 = plt.colorbar(cs2, ax=ax2)
ax2.plot(eta_g, epsilon_n, marker="x", color="red", ls="none")
ax2.set_xlabel(r"$\eta_g$", fontsize=15)
ax2.set_ylabel(r"$\epsilon_n$", fontsize=15)
ax2.set_title("GPR prediction uncertainty")

plt.tight_layout()
plt.show()


# get the leave-one-out predictions
mu_loo, sig_loo = gpr.loo_predictions()

# plot the loo-predictions against the observed function values
plt.errorbar(
    objective_fn, mu_loo, yerr=sig_loo, marker="o", markerfacecolor="none", ls="none"
)
plt.plot([-4, -0.0], [-4, -0], ls="dashed", color="black")
plt.xlabel("observed LOO objective value")
plt.ylabel("prediction of LOO objective value")
plt.title("Leave-one-out cross validation")
plt.xlim([-2, -0.25])
plt.ylim([-2, -0.25])
plt.grid()
plt.tight_layout()
plt.show()
