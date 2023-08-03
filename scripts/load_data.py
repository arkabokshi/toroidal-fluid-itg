"""
Contains methods for the importing of output data produced by the
automation.py script. 
"""
import csv
import numpy as np


def from_csv(file_path: str):
    """
    The given path should point to a csv file that contain the column headers:
    Iteration #  |  eta_g  |  epsilon_n  |  Shear  |  FWHM [radians]

    I.e. a csv file produced by the automation.py script.
    """
    eta_g = []
    epsilon_n = []
    shear = []
    fwhm = []

    with open(file_path, newline="", encoding="utf-8") as csvfile:
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
    eta_g = np.array(eta_g)
    epsilon_n = np.array(epsilon_n)
    shear = np.array(shear)
    fwhm = np.array(fwhm)

    # Convert FWHM data to the objective function value for clarity.
    # Objective function was simply searching for the negative of the FWHM.
    objective_fn = fwhm * -1

    return (
        eta_g.flatten(),
        epsilon_n.flatten(),
        shear.flatten(),
        objective_fn.flatten(),
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
