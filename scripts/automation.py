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

LOAD_DATA = False

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
    x = array([[2.0, 25.0], [2.0, 1.0], [3.0, 15.0], [4.0, 30.0], [5.0, 20.0]])
    y = empty(len(x))

    for j, params in enumerate(x):
        # Evaluate the new point
        run_path = f"{BASE_RUN_PATH}init_{j + 1}"
        poloidal_width = simulation.run(
            eta_g=params[0], epsilon_n=0.08, shear=params[1], run_path=run_path
        )
        y[j] = objective_function(poloidal_width)

# Define bounds for the optimisation of [ eta_g, epsilon_n ]
# These are implied in the initial values of x, but need to be in this form,
# i.e. as tuples, for the GPO class.
# bounds = [(1.0, 5.0), (0.01, 0.1)]  # eta_g, epsilon_n
bounds = [(2.0, 5.0), (1, 30)]  # eta_g, shear

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
    csv_writer = csv.writer(outcsv)
    csv_writer.writerow([r"Run \#", r"$\eta_g$", "shear", "FWHM [radians]"])

    # Add all of the init values to the output
    for k, params in enumerate(x):
        # The values in the y array have already been passed through the
        # objective function, so I need to convert them back to FWHM.
        csv_writer.writerow(
            [
                f"Init {k + 1}",
                str(params[0]),
                str(params[1]),
                str(objective_function(y[k])),
            ]
        )

    # This is where the magic happens...
    for i in range(95):
        # Request the proposed evaluation
        new_x = GPO.propose_evaluation()

        # Evaluate the new point
        run_path = f"{BASE_RUN_PATH}{i + 1}"
        poloidal_width = simulation.run(
            eta_g=new_x[0], epsilon_n=0.08, shear=new_x[1], run_path=run_path
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
        csv_writer.writerow(
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
