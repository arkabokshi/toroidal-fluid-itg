"""
Contains classes pertaining to different simulation models.
"""
import math
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import fit


class ToroidalFluidITG:
    def __init__(
        self,
        run_path="T_F_ITG",
        eta_g=1,
        epsilon_n=0.01,
        shear=25.0,
        input_file_backup="",
    ) -> None:
        self.run_path = run_path
        self.eta_g = eta_g
        self.epsilon_n = epsilon_n
        self.shear = shear
        self.input_file_backup = input_file_backup

    def run(self, eta_g, epsilon_n, shear, run_path) -> float:
        self.eta_g = eta_g
        self.epsilon_n = epsilon_n
        self.shear = shear
        self.run_path = run_path

        self.set_variables()

        subprocess.run("make clean", shell=True, check=True)
        subprocess.run("make itg", shell=True, check=True)

        self.reset_variables()

        return self.get_poloidal_width()

    def set_variables(self) -> None:
        # Read in the file
        with open("inputdata.f90", "r", encoding="utf-8") as file:
            file_data = file.read()
            self.input_file_backup = file_data

        print(self.eta_g, self.epsilon_n, self.run_path)
        # Replace the target string
        file_data = file_data.replace("%ETA_G%", str(self.eta_g))
        file_data = file_data.replace("%EPSILON_N%", str(self.epsilon_n))
        file_data = file_data.replace("%SHEAR%", str(self.shear))
        file_data = file_data.replace("%RUN_PATH%", str(self.run_path))

        # Write the file out again
        with open("inputdata.f90", "w", encoding="utf-8") as file:
            file.write(file_data)

    def reset_variables(self) -> None:
        with open("inputdata.f90", "w", encoding="utf-8") as file:
            file.write(self.input_file_backup)

    def get_poloidal_width(self, plot=False) -> float:
        param = np.loadtxt(f"./{self.run_path}/parameters.txt")

        for file in os.listdir(self.run_path):
            # Loop through the output data to find the time step of the last
            # output file. In most cases this will probably be "20000", but it
            # could be Lower if the code exited via the gamma tolerance.
            if ".dat" in file and file.replace(".dat", "").isnumeric():
                time_step = file.replace(".dat", "")[1:]

        real_field_end = np.fromfile(f"./{self.run_path}/1{time_step}.dat")
        imag_field_end = np.fromfile(f"./{self.run_path}/2{time_step}.dat")

        length = int(param[2])
        num_modes = int(param[3])
        initial_mode = int(param[4])
        final_mode = int(param[5])

        real_final_modes = np.reshape(
            real_field_end, (num_modes, length)
        )  # / np.max(np.abs(real_field_end))
        imag_final_modes = np.reshape(
            imag_field_end, (num_modes, length)
        )  # / np.max(np.abs(imag_field_end))
        complex_final_modes = np.zeros(real_final_modes.shape, dtype=complex)
        complex_final_modes.real = real_final_modes
        complex_final_modes.imag = imag_final_modes

        # --------------------------- #
        # Routine for plotting phi(x) #
        # --------------------------- #

        num_theta = 720

        potential = np.zeros((length, num_theta))
        real_local_potential = np.zeros(length)
        imag_local_potential = np.zeros(length)

        for jj in range(num_theta):
            theta = float(jj) / (num_theta - 1) * 2.0 * math.pi
            real_local_potential = 0.0
            imag_local_potential = 0.0
            for m in range(initial_mode, final_mode + 1):
                real_local_potential += real_final_modes[
                    m - initial_mode, :
                ] * math.cos(m * theta) + imag_final_modes[
                    m - initial_mode, :
                ] * math.sin(
                    m * theta
                )
                imag_local_potential += imag_final_modes[
                    m - initial_mode, :
                ] * math.cos(m * theta) - real_final_modes[
                    m - initial_mode, :
                ] * math.sin(
                    m * theta
                )
            potential[:, jj] = np.sqrt(
                real_local_potential**2 + imag_local_potential**2
            )

        max_index = np.unravel_index(np.argmax(potential, axis=None), potential.shape)

        # This is the value that will be returned if the result is deemed anomalous.
        penalty_value = 2 * np.pi

        # Assume the FWHM is anomalous, until proven otherwise.
        is_anomalous = True
        reason = "Peak not centered"

        if 600 <= max_index[0] <= 1000:
            # An index of 800 would indicate that the peak value occurs in the centre
            # of the poloidal section.
            is_anomalous = False

        pol_section = potential[max_index[0], :]
        pol_section_centred = np.roll(pol_section, int(pol_section.size / 2))
        radians = np.linspace(0, 2 * np.pi, pol_section.size)

        peaks, _ = find_peaks(pol_section_centred, prominence=1)

        if len(peaks) > 1:
            # Finding more than one peak is an indication of an unconverged mode.
            is_anomalous = True
            reason = "Too many peaks"

        # Guess is an array containing [y_0, A, x_0, c]
        gaussian_guess = [0.0, np.amax(pol_section_centred), np.pi, np.pi / 2]
        # pylint: disable=unbalanced-tuple-unpacking
        popt, pcov = curve_fit(
            fit.gaussian, radians, pol_section_centred, p0=gaussian_guess
        )
        fwhm = fit.fwhm(popt[2], popt[3])
        # Write out the FWHM
        with open(f"./{self.run_path}/fwhm.txt", "w", encoding="utf-8") as file:
            file.write(str(fwhm))

        # ------------------------------- #
        # Routine for poloidal width plot #
        # ------------------------------- #

        if plot:
            # Generate the fit using the newly fitted guess parameters
            gaussian_fit = fit.gaussian(radians, *popt)

            # Plot the simulation data
            plt.plot(radians, pol_section_centred, label="toroidal-fluid-itg")

            # Plot the Gaussian fit
            plt.plot(radians, gaussian_fit, lw=2, label="Gaussian Fit")

            # Set x-ticks to be from 0 -> 2 * pi
            tick_pos = [0, np.pi, 2 * np.pi]
            labels = ["0", r"$\pi$", r"$2\pi$"]
            plt.xticks(tick_pos, labels)

            plt.xlabel("Poloidal Angle [radians]")
            plt.ylabel("Intensity [arb. units]")
            plt.legend()

            plt.savefig(f"{self.run_path}_polWidth.png", format="png")
            plt.clf()

        print("Run path:", self.run_path)
        print("Measured FWHM:", fwhm, "radians")
        print("is_anomalous?", is_anomalous)
        if is_anomalous:
            print("  -> why?", reason)
        return penalty_value if is_anomalous else fwhm
