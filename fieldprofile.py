from math import pi, sin, cos
import numpy as np
import matplotlib.pyplot as plt

# ------------ #
# Loading data #
# ------------ #

RUN_PATH = "./run_3"

tstamp = input("Enter last 5-digits of time-slice:")
tstamp = int(tstamp)

param = np.loadtxt(RUN_PATH + "/parameters.txt")
xx = np.fromfile(RUN_PATH + "/xx.dat")
realfieldstart = np.fromfile(RUN_PATH + "/100000.dat")
imagfieldstart = np.fromfile(RUN_PATH + "/200000.dat")
realfieldend = np.fromfile(RUN_PATH + "/" + str(100000 + tstamp) + ".dat")
imagfieldend = np.fromfile(RUN_PATH + "/" + str(200000 + tstamp) + ".dat")

ns = param[0]
dt = param[1]
length = int(param[2])
NumModes = int(param[3])
InitialMode = int(param[4])
FinalMode = int(param[5])

RealInitialModes = np.reshape(realfieldstart, (NumModes, length))
RealFinalModes = np.reshape(
    realfieldend, (NumModes, length)
)  # / np.max(np.abs(realfieldend))
ImagFinalModes = np.reshape(
    imagfieldend, (NumModes, length)
)  # / np.max(np.abs(imagfieldend))
ComplexFinalModes = np.zeros(RealFinalModes.shape, dtype=complex)
ComplexFinalModes.real = RealFinalModes
ComplexFinalModes.imag = ImagFinalModes

# --------------------------- #
# Routine for plotting phi(x) #
# --------------------------- #

NUM_THETA = 720

Potential = np.zeros((NUM_THETA, length))
RealLocalPotential = np.zeros(length)
ImagLocalPotential = np.zeros(length)

for jj in range(NUM_THETA):
    Theta = float(jj) / (NUM_THETA - 1) * 2.0 * pi
    RealLocalPotential = 0.0
    ImagLocalPotential = 0.0
    for m in range(InitialMode, FinalMode + 1):
        RealLocalPotential += RealFinalModes[m - InitialMode, :] * cos(
            m * Theta
        ) + ImagFinalModes[m - InitialMode, :] * sin(m * Theta)
        ImagLocalPotential += ImagFinalModes[m - InitialMode, :] * cos(
            m * Theta
        ) - RealFinalModes[m - InitialMode, :] * sin(m * Theta)
    Potential[jj, :] = np.sqrt(RealLocalPotential**2 + ImagLocalPotential**2)
    # Potential[jj,:] = RealLocalPotential

# figure(1)
# for ii in range(NumModes):
# plot(xx,RealInitialModes[ii,:])


modeamp = np.abs(ComplexFinalModes).max(axis=1)
plt.plot(modeamp)

fig2 = plt.figure(2, figsize=(4, 4))
for ii in range(NumModes):
    plt.plot(xx, RealFinalModes[ii, :], linewidth=0.5)
    # plot(xx,np.abs(ComplexFinalModes[ii,:]),linewidth=0.5)
plt.ylabel("$\mathfrak{R}[\phi_m(y)]$", fontsize=25.0)
plt.xlabel("$y$ (arb. units)", fontsize=25.0, fontname="serif")
plt.xticks(fontsize=20.0, fontname="serif")
plt.yticks(fontsize=20.0, fontname="serif")

plt.show()
fig2.savefig("realfieldtoroid.jpg", format="jpeg", dpi=250, bbox_inches="tight")


# ---------------------------- #
# Routine for 2D Poloidal plot #
# ---------------------------- #
POLOIDAL_PLOT = True


if POLOIDAL_PLOT:
    INNER_RADIUS = 0.25
    OUTER_RADIUS = 0.5

    rr = np.linspace(INNER_RADIUS, OUTER_RADIUS, length)
    pa = np.linspace(0.0, 2.0 * pi, NUM_THETA)
    yy = np.zeros((NUM_THETA, length))
    xx = np.zeros((NUM_THETA, length))
    for jj in range(NUM_THETA):
        yy[jj, :] = rr * sin(pa[jj])
        xx[jj, :] = rr * cos(pa[jj])

    levels = np.linspace(np.min(Potential), np.max(Potential), 50)
    fig4 = plt.figure(4, figsize=(4, 4))
    ax = fig4.add_subplot(1, 1, 1)
    ax.set_aspect("equal")
    plt.contourf(xx, yy, Potential, levels=levels, cmap="jet")
    plt.xticks(fontsize=20.0, fontname="serif")
    plt.yticks(fontsize=20.0, fontname="serif")
    plt.ylabel("$Z(m)$", fontsize=25.0)
    plt.xlabel("$R(m)$", fontsize=25.0, fontname="normal")
    fig4.savefig("isolatedpotential.png", format="png", dpi=250, bbox_inches="tight")
    plt.show()
