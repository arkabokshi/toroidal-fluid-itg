"""
    This script automates the GPO workflow.
    Run this script from the root directory of the toroidal-fluid-itg repository.
    Produces a .csv file that summarises the results.
"""
import csv
import pickle
import subprocess
from pathlib import Path
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


# Read in a path to store all of the output data
data_path = input("Enter relative path for GPR output storage: ")
Path(data_path).mkdir(parents=True, exist_ok=False)

# Read in a number of iterations
iterations = input(
    "Please enter the number of iterations you would like the GPR algorithm to carry out: "
)
iterations = int(iterations)

# Create an instance of the simulation class
simulation = ToroidalFluidITG()
BASE_RUN_PATH = "GPO_"

LOAD_DATA = False

if LOAD_DATA:
    PREVIOUS_RESULTS = "../Data/GPO_100_3D/GPO_results.csv"
    DATA_COUNT = load_data.get_csv_row_count(PREVIOUS_RESULTS)

    # Initial values, in threes, for [ eta_g, epsilon_n, shear ]
    x = empty((DATA_COUNT, 3))
    y = empty((DATA_COUNT))

    (
        previous_eta_g_inputs,
        previous_epsilon_n_inputs,
        previous_shear_inputs,
        previous_objective_fn_results,
    ) = load_data.from_csv(PREVIOUS_RESULTS)

    for i in range(DATA_COUNT):
        x[i] = [
            previous_eta_g_inputs[i],
            previous_epsilon_n_inputs[i],
            previous_shear_inputs[i],
        ]
        y[i] = previous_objective_fn_results[i]

else:
    # If not loading data from a previous run, at least two data points
    # have to be created in order for the GPR algorithm to proceed.
    x = array(
        [
            # 6 Face centers
            [3.5, 0.01, 15.5],
            [5.0, 0.055, 15.5],
            [3.5, 0.1, 15.5],
            [2.0, 0.055, 15.5],
            [3.5, 0.055, 30.0],
            [3.5, 0.055, 1.0],
            # 8 Corners
            [2.0, 0.01, 1.0],
            [5.0, 0.01, 1.0],
            [2.0, 0.01, 30.0],
            [5.0, 0.01, 30.0],
            [2.0, 0.1, 1.0],
            [5.0, 0.1, 1.0],
            [2.0, 0.1, 30.0],
            [5.0, 0.1, 30.0],
            # 12 Edge centers
            [3.5, 0.01, 1.0],
            [2.0, 0.01, 15.5],
            [3.5, 0.01, 30.0],
            [5.0, 0.01, 15.5],
            [3.5, 0.1, 1.0],
            [2.0, 0.1, 15.5],
            [3.5, 0.1, 30.0],
            [5.0, 0.1, 15.5],
            [2.0, 0.055, 1.0],
            [2.0, 0.055, 30.0],
            [5.0, 0.055, 1.0],
            [5.0, 0.055, 30.0],
            # Centre
            [3.5, 0.055, 15.5],
        ]
    )
    y = empty(len(x))

    for j, params in enumerate(x):
        # Evaluate the new point
        run_path = f"{BASE_RUN_PATH}init_{j + 1}"
        poloidal_width = simulation.run(
            eta_g=params[0], epsilon_n=params[1], shear=params[2], run_path=run_path
        )
        y[j] = objective_function(poloidal_width)

# Define bounds for the optimisation of [ eta_g, epsilon_n, shear ]
# These are implied in the initial values of x, but need to be in this form,
# i.e. as tuples, for the GPO class.
bounds = [(2.0, 5.0), (0.01, 0.1), (1, 30)]

# Create an instance of GpOptimiser
GPO = GpOptimiser(x, y, bounds=bounds, acquisition=UpperConfidenceBound)

RESOLUTION = 100
ETA_G_POINTS = linspace(*bounds[0], RESOLUTION)
EPSILON_N_POINTS = linspace(*bounds[1], RESOLUTION)
SHEAR_POINTS = linspace(*bounds[2], RESOLUTION)

x_gp = empty((RESOLUTION, len(bounds)))

for i in range(RESOLUTION):
    x_gp[i] = [ETA_G_POINTS[i], EPSILON_N_POINTS[i], SHEAR_POINTS[i]]

# Store the current state of the system for plotting later
mu, sig = GPO(x_gp)
means = [mu]
sigmas = [sig]
acquis = [array([GPO.acquisition(k) for k in x_gp])]

# Prepare a csv file to output iteration results to.
with open("GPO_results.csv", "w", newline="", encoding="utf-8") as outcsv:
    csv_writer = csv.writer(outcsv)
    csv_writer.writerow(
        [r"Run \#", r"$\eta_g$", r"$\epsilon_n$", "shear", "FWHM [radians]"]
    )

    # Add all of the init values to the output
    for k, params in enumerate(x):
        # The values in the y array have already been passed through the
        # objective function, so I need to convert them back to FWHM.
        csv_writer.writerow(
            [
                f"init_{k + 1}",
                str(params[0]),
                str(params[1]),
                str(params[2]),
                str(objective_function(y[k])),
            ]
        )

    # This is where the magic happens...
    for i in range(iterations):
        # Request the proposed evaluation
        new_x = GPO.propose_evaluation()

        # Evaluate the new point
        run_path = f"{BASE_RUN_PATH}{i + 1}"
        poloidal_width = simulation.run(
            eta_g=new_x[0], epsilon_n=new_x[1], shear=new_x[2], run_path=run_path
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
            [
                str(i + 1),
                str(new_x[0][0]),
                str(new_x[0][1]),
                str(new_x[0][2]),
                str(poloidal_width),
            ]
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

subprocess.run(f"mv *.pickle *.npy GPO_* {data_path}", shell=True, check=True)