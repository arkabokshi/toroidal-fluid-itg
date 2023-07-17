import matplotlib.pyplot as plt
import numpy as np

# ------------ #
# Loading data #
# ------------ #

RUN_PATH = "./run_3"

param = np.loadtxt(RUN_PATH + "/parameters.txt")
xx = np.loadtxt(RUN_PATH + "/xx.txt")
growthrate = np.fromfile(RUN_PATH + "/growthrate.dat")
frequency = np.fromfile(RUN_PATH + "/frequency.dat")
gammaE_t = np.fromfile(RUN_PATH + "/gammaE_t.dat")

ThetaMaxima = np.fromfile(RUN_PATH + "/thetamaxima.dat")
GlobalOmega = np.fromfile(RUN_PATH + "/globalomega.dat")
Globalgamma = np.fromfile(RUN_PATH + "/globalgamma.dat")

Globalgamma = Globalgamma / 2.0
Globalomega = np.sqrt(GlobalOmega**2 - Globalgamma**2)
Globalomega = -1.0 * Globalomega


NumSteps = int(param[0])
dt = param[1]
length = int(param[2])
NumModes = int(param[3])

time = np.linspace(0.0, NumSteps * dt, NumSteps)
mode_gamma = np.reshape(growthrate, (NumModes, NumSteps))
mode_omega = np.reshape(frequency, (NumModes, NumSteps))


# -------- #
# Analysis #
# -------- #

gamma_max = np.max(Globalgamma)
gamma_max_loc = np.where(Globalgamma == gamma_max)
# gamma_max_loc = (array([25125]),)

print("Max(gamma) loc:", gamma_max_loc)
print("Gamma at gamma_max:", Globalgamma[gamma_max_loc])
print("omega at gamma_max:", Globalomega[gamma_max_loc])
print("theta at gamma_max:", ThetaMaxima[gamma_max_loc])
print("gamma_E at gamma_max:", gammaE_t[gamma_max_loc])
print("")

# print Globalomega[theta0]
# print Globalgamma[theta0]

print("Mean [ gamma ] :", np.mean(Globalgamma[55000:60000]))
print("Mean [ omega ] :", np.mean(Globalomega[55000:60000]))
print("Mean [ theta ] :", np.mean(ThetaMaxima[55000:60000]))


# -------- #
# Plotting #
# -------- #

fig1 = plt.figure(1, figsize=(12, 8))

plt.subplot(221)
# for i in range(NumModes):
#  plot(time,mode_gamma[i,:],linewidth=0.5)

plt.plot(time, Globalgamma, "y--", linewidth=2.0)
plt.plot([0.0, time[-1]], [0.743984083872, 0.743984083872], "k--", linewidth=1.5)
plt.scatter(float(gamma_max_loc[-1]) * dt, Globalgamma[gamma_max_loc])
plt.ylabel(r"$\gamma$", fontsize=15.0)
plt.xlabel("Time", fontsize=15.0, fontname="serif", visible=False)
plt.xticks(fontsize=15.0, fontname="serif", visible=False)
plt.yticks(fontsize=15.0, fontname="serif")
# ylim([0.69,0.74])
# xlim([300.,600.])

plt.subplot(222)
# for i in range(NumModes):
#  plot(time,mode_omega[i,:],linewidth=0.5)
plt.plot(time, Globalomega, "r--", linewidth=2.0)
plt.plot([0.0, time[-1]], [-0.757499193109, -0.757499193109], "k--", linewidth=1.5)
plt.scatter(float(gamma_max_loc[-1]) * dt, Globalomega[gamma_max_loc])
plt.ylabel(r"$\omega$", fontsize=15.0)
plt.xlabel("Time", fontsize=15.0, fontname="serif", visible=False)
plt.xticks(fontsize=15.0, fontname="serif", visible=False)
plt.yticks(fontsize=15.0, fontname="serif")
# ylim([-1.15,-0.95])
# xlim([300.,600.])

plt.subplot(223)
plt.plot(time, gammaE_t, linewidth=2.0)
plt.plot([time[0], time[-1]], [0.0, 0.0])
plt.scatter(float(gamma_max_loc[-1]) * dt, gammaE_t[gamma_max_loc])
plt.ylabel(r"$\gamma_E$", fontsize=15.0)
plt.xlabel("Time", fontsize=15.0, fontname="serif")
plt.xticks(fontsize=15.0, fontname="serif")
plt.yticks(fontsize=15.0, fontname="serif")
# tight_layout(h_pad=0.4)
# ylim([-0.0032,0.0032])
# xlim([300.,600.])

plt.subplot(224)
plt.plot(time, ThetaMaxima, linewidth=2.0)
plt.scatter(float(gamma_max_loc[-1]) * dt, ThetaMaxima[gamma_max_loc])
plt.ylabel(r"$max(\Theta)$ [$\pi$ Radians]", fontsize=15.0)
plt.xlabel("Time", fontsize=15.0, fontname="serif")
plt.xticks(fontsize=15.0, fontname="serif")
plt.yticks(fontsize=15.0, fontname="serif")
# tight_layout(h_pad=0.4)
# ylim([-0.0032,0.0032])
# xlim([300.,600.])

fig1.savefig("omega.png", format="png", dpi=250, bbox_inches="tight")

plt.show()
