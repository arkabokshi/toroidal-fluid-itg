"""
Contains classes pertaining to different simulation models.
"""
import math
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import fit


class ToroidalFluidITG:
    def __init__(
        self, run_path="T_F_ITG", eta_g=1, epsilon_n=0.01, input_file_backup=""
    ) -> None:
        self.run_path = run_path
        self.eta_g = eta_g
        self.epsilon_n = epsilon_n
        self.input_file_backup = input_file_backup

    def run(self, eta_g, epsilon_n, run_path) -> float:
        self.eta_g = eta_g
        self.epsilon_n = epsilon_n
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

        # Replace the target string
        file_data = file_data.replace("%ETA_G%", self.eta_g)
        file_data = file_data.replace("%EPSILON_N%", self.epsilon_n)
        file_data = file_data.replace("%RUN_PATH%", self.run_path)

        # Write the file out again
        with open("inputdata.f90", "w", encoding="utf-8") as file:
            file.write(file_data)

    def reset_variables(self) -> None:
        with open("inputdata.f90", "w", encoding="utf-8") as file:
            file.write(self.input_file_backup)

    def get_poloidal_width(self) -> float:
        param = np.loadtxt(f"../{self.run_path}/parameters.txt")

        real_field_end = np.loadtxt(f"../{self.run_path}/120000.txt")
        imag_field_end = np.loadtxt(f"../{self.run_path}/220000.txt")

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
        complex_final_modes = np.zeros(real_final_modes.shape, dtype=np.complex)
        complex_final_modes.real = real_final_modes
        complex_final_modes.imag = imag_final_modes

        # --------------------------- #
        # Routine for plotting phi(x) #
        # --------------------------- #

        NUM_THETA = 720

        Potential = np.zeros((length, NUM_THETA))
        RealLocalPotential = np.zeros(length)
        ImagLocalPotential = np.zeros(length)

        for jj in range(NUM_THETA):
            Theta = float(jj) / (NUM_THETA - 1) * 2.0 * math.pi
            RealLocalPotential = 0.0
            ImagLocalPotential = 0.0
            for m in range(initial_mode, final_mode + 1):
                RealLocalPotential += real_final_modes[m - initial_mode, :] * math.cos(
                    m * Theta
                ) + imag_final_modes[m - initial_mode, :] * math.sin(m * Theta)
                ImagLocalPotential += imag_final_modes[m - initial_mode, :] * math.cos(
                    m * Theta
                ) - real_final_modes[m - initial_mode, :] * math.sin(m * Theta)
            Potential[:, jj] = np.sqrt(
                RealLocalPotential**2 + ImagLocalPotential**2
            )

        maxIndex = np.unravel_index(np.argmax(Potential, axis=None), Potential.shape)
        polSection = Potential[maxIndex[0], :]
        polSectionCentred = np.roll(polSection, int(polSection.size / 2))
        radians = np.linspace(0, 2 * np.pi, polSection.size)

        # Guess is an array containing [y_0, A, x_0, c]
        gaussian_guess = [0.0, np.amax(polSectionCentred), np.pi, np.pi / 2]
        # pylint: disable=unbalanced-tuple-unpacking
        popt, pcov = curve_fit(
            fit.gaussian, radians, polSectionCentred, p0=gaussian_guess
        )
        fwhm = fit.fwhm(popt[2], popt[3])
        print("FWHM:", fwhm, "radians")
        # Write out the FWHM
        with open(f"../{self.run_path}/fwhm.txt", "w", encoding="utf-8") as file:
            file.write(fwhm)

        # ------------------------------- #
        # Routine for poloidal width plot #
        # ------------------------------- #
        POLOIDAL_WIDTH_PLOT = False

        if POLOIDAL_WIDTH_PLOT:
            # Generate the fit using the newly fitted guess parameters
            gaussian_fit = fit.gaussian(radians, *popt)

            # Plot the simulation data
            plt.plot(radians, polSectionCentred, label="toroidal-fluid-itg")

            # Plot the Gaussian fit
            plt.plot(radians, gaussian_fit, lw=2, label="Gaussian Fit")

            # Set x-ticks to be from 0 -> 2 * pi
            tick_pos = [0, np.pi, 2 * np.pi]
            labels = ["0", r"$\pi$", r"$2\pi$"]
            plt.xticks(tick_pos, labels)

            plt.xlabel("Poloidal Angle [radians]")
            plt.ylabel("Intensity [arb. units]")
            plt.legend()

            plt.show()

        return fwhm
