import numpy as np
from matplotlib.pyplot import *
from numpy import *
from scipy.optimize import *
from numpy.fft import *


# ------------ #
# Loading data #
# ------------ #

runpath = './run_3'

param = np.loadtxt(runpath+'/parameters.txt')
xx = np.fromfile(runpath+'/xx.dat')
growthrate = np.fromfile(runpath+'/growthrate.dat')
frequency  = np.fromfile(runpath+'/frequency.dat')
gammaE_t = np.fromfile(runpath+'/gammaE_t.dat')

ThetaMaxima = np.fromfile(runpath+'/thetamaxima.dat')
GlobalOmega = np.fromfile(runpath+'/globalomega.dat')
Globalgamma = np.fromfile(runpath+'/globalgamma.dat')

Globalgamma = Globalgamma / 2.0
Globalomega = np.sqrt( GlobalOmega**2 - Globalgamma**2 )
Globalomega = -1.0*Globalomega


NumSteps = int(param[0])
dt = param[1]
length = int(param[2])
NumModes = int(param[3])

time = linspace(0.0,NumSteps*dt,NumSteps)
mode_gamma = reshape(growthrate,(NumModes,NumSteps))
mode_omega = reshape(frequency, (NumModes,NumSteps))


# -------- #
# Analysis #
# -------- #

gamma_max = np.max( Globalgamma )
gamma_max_loc = where( Globalgamma==gamma_max )
#gamma_max_loc = (array([25125]),)

print('Max(gamma) loc:',gamma_max_loc)
print('Gamma at gamma_max:', Globalgamma[gamma_max_loc])
print('omega at gamma_max:', Globalomega[gamma_max_loc])
print('theta at gamma_max:', ThetaMaxima[gamma_max_loc])
print('gamma_E at gamma_max:', gammaE_t[gamma_max_loc])
print('')

#print Globalomega[theta0]
#print Globalgamma[theta0]
  
print('Mean [ gamma ] :' , mean(Globalgamma[55000:60000]))
print('Mean [ omega ] :' , mean(Globalomega[55000:60000]))
print('Mean [ theta ] :' , mean(ThetaMaxima[55000:60000]))


# -------- #
# Plotting #
# -------- #

fig1 = figure(1,figsize=(12,8))

subplot(221)
#for i in range(NumModes):
#  plot(time,mode_gamma[i,:],linewidth=0.5)

plot(time,Globalgamma,'y--',linewidth=2.0)
plot([0.0,time[-1]],[0.743984083872,0.743984083872],'k--',linewidth=1.5)
scatter(float(gamma_max_loc[-1])*dt,Globalgamma[gamma_max_loc])
ylabel('$\gamma$',fontsize=15.0)
xlabel('Time',fontsize=15.0,fontname='serif',visible=False)
xticks(fontsize=15.0,fontname='serif',visible=False)
yticks(fontsize=15.0,fontname='serif')
#ylim([0.69,0.74])
#xlim([300.,600.])

subplot(222)
#for i in range(NumModes):
#  plot(time,mode_omega[i,:],linewidth=0.5)
plot(time,Globalomega,'r--',linewidth=2.0)
plot([0.0,time[-1]],[-0.757499193109,-0.757499193109],'k--',linewidth=1.5)
scatter(float(gamma_max_loc[-1])*dt,Globalomega[gamma_max_loc])
ylabel('$\omega$',fontsize=15.0)
xlabel('Time',fontsize=15.0,fontname='serif',visible=False)
xticks(fontsize=15.0,fontname='serif',visible=False)
yticks(fontsize=15.0,fontname='serif')
#ylim([-1.15,-0.95])
#xlim([300.,600.])

subplot(223)
plot(time,gammaE_t,linewidth=2.0)
plot([time[0],time[-1]],[0.0,0.0])
scatter(float(gamma_max_loc[-1])*dt,gammaE_t[gamma_max_loc])
ylabel('$\gamma_E$',fontsize=15.0)
xlabel('Time',fontsize=15.0,fontname='serif')
xticks(fontsize=15.0,fontname='serif')
yticks(fontsize=15.0,fontname='serif')
#tight_layout(h_pad=0.4)
#ylim([-0.0032,0.0032])
#xlim([300.,600.])

subplot(224)
plot(time,ThetaMaxima,linewidth=2.0)
scatter(float(gamma_max_loc[-1])*dt,ThetaMaxima[gamma_max_loc])
ylabel('$max(\Theta)$ [$\pi$ Radians]',fontsize=15.0)
xlabel('Time',fontsize=15.0,fontname='serif')
xticks(fontsize=15.0,fontname='serif')
yticks(fontsize=15.0,fontname='serif')
#tight_layout(h_pad=0.4)
#ylim([-0.0032,0.0032])
#xlim([300.,600.])


fig1.savefig('omega.jpg',format='jpeg',dpi=250,bbox_inches='tight')

show()


