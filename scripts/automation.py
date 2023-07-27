"""
    This script automates the GPO workflow.
    Run this script from the root directory of the toroidal-fluid-itg repository.
    Produces a .csv file that summarises the results.
"""
import csv
import pickle
from numpy import linspace, array, empty, save
from inference.gp import GpOptimiser, UpperConfidenceBound
from models import ToroidalFluidITG
import load_data


def objective_function(pol_width: float) -> float:
    """
    We want to find the minimum poloidal width, but the GPO algorithm searches
    for a maximum, so our objective function needs to return the negative.
    """
    return pol_width * -1.0


# Create an instance of the simulation class
simulation = ToroidalFluidITG()
BASE_RUN_PATH = "GPO_"

LOAD_DATA = True

if LOAD_DATA:
    DATA_COUNT = load_data.get_csv_row_count(
        "../Data/GPO_50_Iterations_Run/GPO_results_appended.csv"
    )

    # Initial values, in pairs, for [ eta_g, epsilon_n ]
    x = empty((DATA_COUNT, 2))
    y = empty((DATA_COUNT))

    (
        previous_eta_g_inputs,
        previous_epsilon_n_inputs,
        previous_objective_fn_results,
    ) = load_data.from_csv("../Data/GPO_50_Iterations_Run/GPO_results_appended.csv")

    for i, value in enumerate(previous_objective_fn_results):
        x[i] = [previous_eta_g_inputs[i], previous_epsilon_n_inputs[i]]
        y[i] = previous_objective_fn_results[i]

else:
    # If not loading data from a previous run, at least two data points
    # have to be created in order for the GPR algorithm to proceed.
    x = array([[1.0, 0.01], [5, 0.1]])

    LOWER_BOUND_RUN_PATH = BASE_RUN_PATH + "0_lower_bound"
    UPPER_BOUND_RUN_PATH = BASE_RUN_PATH + "0_upper_bound"

    y1 = objective_function(
        simulation.run(eta_g=x[0][0], epsilon_n=x[0][1], run_path=LOWER_BOUND_RUN_PATH)
    )

    y2 = objective_function(
        simulation.run(eta_g=x[1][0], epsilon_n=x[1][1], run_path=UPPER_BOUND_RUN_PATH)
    )

    y = array([y1, y2])

# Define bounds for the optimisation of [ eta_g, epsilon_n ]
# These are implied in the initial values of x, but need to be in this form,
# i.e. as tuples, for the GPO class.
bounds = [(1.0, 5.0), (0.01, 0.1)]

# Create an instance of GpOptimiser
GPO = GpOptimiser(x, y, bounds=bounds, acquisition=UpperConfidenceBound)

RESOLUTION = 100
ETA_G_POINTS = linspace(*bounds[0], RESOLUTION)
EPSILON_N_POINTS = linspace(*bounds[1], RESOLUTION)
x_gp = empty((RESOLUTION, 2))

for i in range(RESOLUTION):
    x_gp[i] = [ETA_G_POINTS[i], EPSILON_N_POINTS[i]]

# Store the current state of the system for plotting later
mu, sig = GPO(x_gp)
means = [mu]
sigmas = [sig]
acquis = [array([GPO.acquisition(k) for k in x_gp])]

# Prepare a csv file to output iteration results to.
with open("GPO_results.csv", "w", newline="", encoding="utf-8") as outcsv:
    writer = csv.writer(outcsv)
    writer.writerow([r"Run \#", r"$\eta_g$", r"$\epsilon_n$", "FWHM [radians]"])

    # This is where the magic happens...
    for i in range(30):
        # Request the proposed evaluation
        new_x = GPO.propose_evaluation()

        # Evaluate the new point
        run_path = f"{BASE_RUN_PATH}{i + 1}"
        poloidal_width = simulation.run(
            eta_g=new_x[0], epsilon_n=new_x[1], run_path=run_path
        )
        new_y = objective_function(poloidal_width)

        # Update the gaussian process with the new information
        # Note: this function mutates the shape of new_x
        GPO.add_evaluation(new_x, new_y)

        # Store the current state of the system for plotting later
        mu, sig = GPO(x_gp)
        means.append(mu)
        sigmas.append(sig)
        acquis.append(array([GPO.acquisition(k) for k in x_gp]))
        writer.writerow(
            [str(i + 1), str(new_x[0][0]), str(new_x[0][1]), str(poloidal_width)]
        )

# Write out all data for easy loading in the gpo_plot.py script
with open("means.npy", "wb") as f:
    save(f, means)
with open("sigmas.npy", "wb") as f:
    save(f, sigmas)
with open("acquis.npy", "wb") as f:
    save(f, acquis)
with open("x_gp.npy", "wb") as f:
    save(f, x_gp)
with open("GPO.pickle", "wb") as file:
    pickle.dump(GPO, file)
