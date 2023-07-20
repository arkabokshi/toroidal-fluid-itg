"""
    This script should automate the GPR workflow.
"""
from numpy import linspace, array
from inference.gp import GpOptimiser, UpperConfidenceBound
from models import ToroidalFluidITG


def objective_function(pol_width: float) -> float:
    """
    We want to find the minimum poloidal width, so our
    objective function needs to return its inverse.
    """
    return 1.0 / pol_width


# Initial values for [ eta_g, epsilon_n ]
x = array([[1.0], [0.01]])

# Create an instance of the simulation class
simulation = ToroidalFluidITG()

BASE_RUN_PATH = "GPO_"

# Initial result for the poloidal width
INITIAL_POL_WIDTH = simulation.run(
    eta_g=x[0][0], epsilon_n=x[1][0], run_path=BASE_RUN_PATH + "0"
)
y = array([objective_function(INITIAL_POL_WIDTH)])

# define bounds for the optimisation of [ eta_g, epsilon_n ]
bounds = [(1.0, 5.0), (0.01, 0.1)]

# Create an instance of GpOptimiser
GPO = GpOptimiser(x, y, bounds=bounds, acquisition=UpperConfidenceBound)

# Store the current state of the system for plotting later
resolution = [5, 10]  # for [ eta_g, epsilon_n ]
x_gp = array([linspace(*bounds[0], resolution[0]), linspace(*bounds[1], resolution[1])])
mu, sig = GPO(x_gp)
means = [mu]
sigmas = [sig]
acquis = [array([GPO.acquisition(k) for k in x_gp])]


for i in range(5):
    # request the proposed evaluation
    new_x = GPO.propose_evaluation()
    print("new_x:", new_x)

    # evaluate the new point
    POLOIDAL_WIDTH = simulation.run(
        eta_g=new_x[0][0], epsilon_n=new_x[1][0], run_path=f"{BASE_RUN_PATH}{i + 1}"
    )
    new_y = objective_function(POLOIDAL_WIDTH)
    print("new_y", new_y)

    # update the gaussian process with the new information
    GPO.add_evaluation(new_x, new_y)

    # Store the current state of the system for plotting later
    mu, sig = GPO(x_gp)
    means.append(mu)
    sigmas.append(sig)
    acquis.append(array([GPO.acquisition(k) for k in x_gp]))

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, axes = plt.subplots(2, 3, gridspec_kw={"height_ratios": [1, 1]}, figsize=(15, 10))

for k, i, j in [(i, *divmod(i, 3)) for i in range(6)]:
    aq = acquis[k] - acquis[k].min()
    proposal = x_gp[aq.argmax()]

    divider = make_axes_locatable(axes[i, j])
    acq_ax = divider.append_axes("bottom", size="40%", pad=0.0)

    axes[i, j].plot(
        GPO.x[: k + 4], GPO.y[: k + 4], "o", c="red", label="observations", zorder=5
    )
    axes[i, j].plot(x_gp, means[k], lw=2, c="blue", label="GP prediction")
    axes[i, j].fill_between(
        x_gp,
        (means[k] - 2 * sigmas[k]),
        y2=(means[k] + 2 * sigmas[k]),
        color="blue",
        alpha=0.15,
        label="95% confidence interval",
    )
    axes[i, j].set_ylim([-1.5, 4])
    axes[i, j].set_xlim([-8, 8])
    axes[i, j].set_xticks([])
    axes[i, j].text(
        0.99,
        0.98,
        f"iteration {k + 1}",
        transform=axes[i, j].transAxes,
        verticalalignment="top",
        horizontalalignment="right",
        fontsize=13,
    )

    acq_ax.fill_between(x_gp, 0.9 * aq / aq.max(), color="green", alpha=0.15)
    acq_ax.plot(x_gp, 0.9 * aq / aq.max(), color="green", label="acquisition function")
    acq_ax.plot(
        [proposal] * 2, [0.0, 1.0], c="green", ls="dashed", label="acquisition maximum"
    )
    acq_ax.set_ylim([0, 1])
    acq_ax.set_xlim([-8, 8])
    acq_ax.set_yticks([])
    acq_ax.set_xlabel("x")

    if k == 0:
        axes[i, j].legend()
        acq_ax.legend()

plt.tight_layout()
plt.savefig("GPO.png", format="png")
plt.show()
