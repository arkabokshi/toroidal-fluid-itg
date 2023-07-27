"""
Contains methods for the importing of output data produced by the
automation.py script. 
"""
import csv
import numpy as np


def from_csv(file_path: str):
    """
    The given path should point to a csv file that contain the column headers:
    Iteration #  |  eta_g  |  epsilon_n  |  FWHM [radians]

    I.e. a csv file produced by the automation.py script.
    """
    eta_g = []
    epsilon_n = []
    fwhm = []

    with open(file_path, newline="", encoding="utf-8") as csvfile:
        results = csv.reader(csvfile, delimiter=",")
        for i, row in enumerate(results):
            if i == 0:
                # Don't read the headers.
                continue

            eta_g.append(float(row[1]))
            epsilon_n.append(float(row[2]))
            fwhm.append(float(row[3]))

    # Convert lists to numpy arrays
    eta_g = np.array(eta_g)
    epsilon_n = np.array(epsilon_n)
    fwhm = np.array(fwhm)

    # Convert FWHM data to the objective function value for clarity.
    # Objective function was simply searching for the negative of the FWHM.
    objective_fn = fwhm * -1

    # Create a list of indices that contain physical results
    physical_results = np.argwhere(objective_fn > (-2 * np.pi))

    return (
        eta_g[physical_results].flatten(),
        epsilon_n[physical_results].flatten(),
        objective_fn[physical_results].flatten(),
    )


def get_csv_row_count(file_path: str) -> int:
    """
    Counts and returns the number of entries in a csv file.
    Assumes that the first row consists of column headers, and ignores it.

    :param file_path: \
        The path to the csv file.

    :return: \
        The number of entries in the csv file (excluding the headers).
    """
    with open(file_path, newline="", encoding="utf-8") as csvfile:
        return sum(1 for line in csvfile) - 1
