"""
    Plots the data created by the GpOptimiser class used in automation.py.
    Example taken from the inference tools demo:
    https://github.com/C-bowman/inference-tools/blob/master/demos/gp_optimisation_demo.ipynb
"""
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import load

with open("means.npy", "rb") as f:
    means = load(f)
with open("sigmas.npy", "rb") as f:
    sigmas = load(f)
with open("acquis.npy", "rb") as f:
    acquis = load(f)
with open("x_gp.npy", "rb") as f:
    x_gp = load(f)

with open("GPO.pickle", "rb") as file:
    GPO = pickle.load(file)

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
