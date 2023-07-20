from math import pi, cos, sin
import numpy as np
import matplotlib.pyplot as plt

# ------------ #
# Loading data #
# ------------ #

RUN_PATH = "./run_3"

parameters = np.loadtxt(RUN_PATH + "/parameters.txt")
videoparam = np.loadtxt(RUN_PATH + "/videoparam.txt")

xx = np.fromfile(RUN_PATH + "/xx.dat")
growthrate = np.fromfile(RUN_PATH + "/growthrate.dat")
frequency = np.fromfile(RUN_PATH + "/frequency.dat")
Globalgamma = np.fromfile(RUN_PATH + "/globalgamma.dat")
Globalgamma = Globalgamma / 2.0
gammaE_t = np.fromfile(RUN_PATH + "/gammaE_t.dat")

NumSteps = int(parameters[0])
dt = parameters[1]
Length = int(parameters[2])
NumModes = int(parameters[3])
InitialMode = int(parameters[4])
FinalMode = int(parameters[5])

NumFiles = int(videoparam[0])
DelPrint = int(videoparam[1])


time = np.linspace(0.0, NumSteps * dt, NumSteps)
time = time / 1.0e4
mode_gamma = np.reshape(growthrate, (NumModes, NumSteps))
mode_omega = np.reshape(frequency, (NumModes, NumSteps))

RealData = np.zeros((NumModes, Length, NumFiles))
ImagData = np.zeros((NumModes, Length, NumFiles))

for tt in range(NumFiles):
    tstamp = DelPrint * tt
    TempData = np.fromfile(RUN_PATH + "/" + str(100000 + tstamp) + ".dat")
    RealData[:, :, tt] = np.reshape(TempData, (NumModes, Length))
    TempData = np.fromfile(RUN_PATH + "/" + str(200000 + tstamp) + ".dat")
    ImagData[:, :, tt] = np.reshape(TempData, (NumModes, Length))


# --------------------------- #
# Routine for plotting phi(x) #
# --------------------------- #
print("------------------------------------")
print("Reconstructing poloidal potential...")
print("------------------------------------")

NUM_THETA = 720

Potential = np.zeros((NUM_THETA, Length, NumFiles))
RealLocalPotential = np.zeros(Length)
ImagLocalPotential = np.zeros(Length)


D_PLOT = 1

for tt in range(0, NumFiles, D_PLOT):
    print("Timeslice:", tt + 1, "/", NumFiles)

    for jj in range(NUM_THETA):
        Theta = float(jj) / (NUM_THETA - 1) * 2.0 * pi

        RealLocalPotential = 0.0
        ImagLocalPotential = 0.0

        for m in range(InitialMode, FinalMode + 1):
            RealLocalPotential += RealData[m - InitialMode, :, tt] * cos(
                m * Theta
            ) + ImagData[m - InitialMode, :, tt] * sin(m * Theta)
            ImagLocalPotential += ImagData[m - InitialMode, :, tt] * cos(
                m * Theta
            ) - RealData[m - InitialMode, :, tt] * sin(m * Theta)

        Potential[jj, :, tt] = np.sqrt(
            RealLocalPotential**2 + ImagLocalPotential**2
        )
        # Potential[jj,:,tt] = RealLocalPotential


# --------------------------------- #
# Routine for 2D Poloidal animation #
# --------------------------------- #

INNER_RADIUS = 0.25
OUTER_RADIUS = 0.5


rr = np.linspace(INNER_RADIUS, OUTER_RADIUS, Length)
pa = np.linspace(0.0, 2.0 * pi, NUM_THETA)
yy = np.zeros((NUM_THETA, Length))
xx = np.zeros((NUM_THETA, Length))
for jj in range(NUM_THETA):
    yy[jj, :] = rr * sin(pa[jj])
    xx[jj, :] = rr * cos(pa[jj])


fig = plt.figure(figsize=(16, 5), dpi=1000)
# plt.hold(False)
ims = []

for tt in range(0, NumFiles, D_PLOT):
    levels = np.linspace(np.min(Potential[:, :, tt]), np.max(Potential[:, :, tt]), 50)
    timepoint = DelPrint * tt * dt / 1.0e4

    yll = Globalgamma[-1] * 0.90
    yul = Globalgamma[-1] * 1.10

    plt.subplot(1, 3, 2)

    plt.contourf(xx, yy, Potential[:, :, tt], levels=levels, cmap="jet")
    # plt.xticks(fontsize=20.0,fontname='serif')
    # plt.yticks(fontsize=20.0,fontname='serif')
    # plt.imshow(Potential[:,:,tt].transpose())
    plt.ylabel("Z (arb. units)", fontsize=15.0, fontname="serif")
    plt.xlabel("R (arb. units)", fontsize=15.0, fontname="serif")
    plt.title(r"$|\phi_1(x,\Theta)|$", fontsize=20.0, fontname="serif")

    plt.subplot(1, 3, 3)
    for i in range(20, 60):
        plt.plot(time, mode_gamma[i, :], linewidth=0.5)
    plt.plot(time, Globalgamma, "y--", linewidth=2.5)
    # plt.plot([0.0,time[-1]],[0.743984083872,0.743984083872],'k--',linewidth=1.5)
    plt.plot([timepoint, timepoint], [yll, yul], "k", linewidth=1.5)
    # plt.ylabel('$\gamma$',fontsize=15.0)
    plt.xlabel("Time (ms)", fontsize=15.0, fontname="serif", visible=True)
    plt.title(r"$\gamma / \omega_e^*$", fontsize=20.0, fontname="serif")
    plt.xticks(
        np.linspace(0.0, np.max(time), 3), fontsize=15.0, fontname="serif", visible=True
    )
    plt.yticks(fontsize=15.0, fontname="serif")
    plt.ylim([yll, yul])

    plt.subplot(1, 3, 1)
    plt.plot(time, gammaE_t, "b", linewidth=2.0)
    plt.plot([timepoint, timepoint], [-0.004, 0.004], "k", linewidth=1.5)
    # plt.ylabel('$\gamma_E$',fontsize=15.0,visible='False')
    plt.title(r"$\gamma_E = d\Omega/dq $", fontsize=20.0, fontname="serif")
    plt.xlabel("Time (ms)", fontsize=15.0, fontname="serif")
    plt.xticks(np.linspace(0.0, np.max(time), 3), fontsize=15.0, fontname="serif")
    plt.yticks(fontsize=15.0, fontname="serif")
    plt.ylim([-0.004, 0.004])

    # plt.tight_layout()

    fig.savefig(
        "./images/" + str(tt) + ".jpg", format="jpeg", dpi=250, bbox_inches="tight"
    )

    print("Saving:", tt + 1, "/", NumFiles)
    plt.clf()
