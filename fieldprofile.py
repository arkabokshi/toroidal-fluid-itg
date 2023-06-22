from numpy import *
import numpy as np
from matplotlib.pyplot import *
from scipy.interpolate import griddata
import math

# ------------ #
# Loading data #
# ------------ #

tstamp = input('Enter last 5-digits of time-slice:')
tstamp = int(tstamp)

param = loadtxt('./parameters.txt')
xx = loadtxt('./xx.txt')
realfieldstart = loadtxt('./100000.txt')
imagfieldstart = loadtxt('./200000.txt')
realfieldend = loadtxt('./'+str(100000+tstamp)+'.txt')
imagfieldend = loadtxt('./'+str(200000+tstamp)+'.txt')

ns = param[0]
dt = param[1]
length = int(param[2])
NumModes = int(param[3])
InitialMode = int(param[4])
FinalMode = int(param[5])

RealInitialModes = reshape(realfieldstart,(NumModes,length))
RealFinalModes = reshape(realfieldend,(NumModes,length)) #/ np.max(np.abs(realfieldend))
ImagFinalModes = reshape(imagfieldend,(NumModes,length)) #/ np.max(np.abs(imagfieldend))
ComplexFinalModes=np.zeros(RealFinalModes.shape,dtype=np.complex)
ComplexFinalModes.real=RealFinalModes
ComplexFinalModes.imag=ImagFinalModes

# --------------------------- #
# Routine for plotting phi(x) #
# --------------------------- #

NumTheta = 720

Potential = np.zeros((NumTheta,length))
RealLocalPotential = np.zeros(length)
ImagLocalPotential = np.zeros(length)

for jj in range(NumTheta):
	Theta = float(jj)/(NumTheta-1)*2.0*pi
	RealLocalPotential = 0.0
	ImagLocalPotential = 0.0
	for m in range(InitialMode,FinalMode+1):
		RealLocalPotential += RealFinalModes[m-InitialMode,:]*cos(m*Theta) + ImagFinalModes[m-InitialMode,:]*sin(m*Theta)
		ImagLocalPotential += ImagFinalModes[m-InitialMode,:]*cos(m*Theta) - RealFinalModes[m-InitialMode,:]*sin(m*Theta)
	Potential[jj,:] = sqrt( RealLocalPotential**2 + ImagLocalPotential**2 )
	#Potential[jj,:] = RealLocalPotential

#figure(1)
#for ii in range(NumModes):
  #plot(xx,RealInitialModes[ii,:])


modeamp=np.abs(ComplexFinalModes).max(axis=1)
plot(modeamp)

fig2 = figure(2,figsize=(4,4))
for ii in range(NumModes):
  plot(xx,RealFinalModes[ii,:],linewidth=0.5)
  #plot(xx,np.abs(ComplexFinalModes[ii,:]),linewidth=0.5)
ylabel('$\mathfrak{R}[\phi_m(y)]$',fontsize=25.0)
xlabel('$y$ (arb. units)',fontsize=25.0,fontname='serif')
xticks(fontsize=20.0,fontname='serif')
yticks(fontsize=20.0,fontname='serif')

show()
fig2.savefig('realfieldtoroid.jpg',format='jpeg',dpi=250,bbox_inches='tight')


# ---------------------------- #
# Routine for 2D Poloidal plot #
# ---------------------------- #
PoloidalPlot = True


if PoloidalPlot:
  
    InnerRadius = 0.25
    OuterRadius = 0.5

    rr = np.linspace(InnerRadius,OuterRadius,length) 
    pa = np.linspace(0.0,2.0*pi,NumTheta)
    yy = np.zeros((NumTheta,length))
    xx = np.zeros((NumTheta,length))
    for jj in range(NumTheta):
        yy[jj,:] = rr * sin(pa[jj])
        xx[jj,:] = rr * cos(pa[jj])
  
    levels = np.linspace(np.min(Potential),np.max(Potential),50)
    fig4 = figure(4,figsize=(4,4))
    ax = fig4.add_subplot(1,1,1)
    ax.set_aspect('equal')
    contourf(xx,yy,Potential,levels=levels,cmap=cm.jet)
    xticks(fontsize=20.0,fontname='serif')
    yticks(fontsize=20.0,fontname='serif')
    ylabel('$Z(m)$',fontsize=25.0)
    xlabel('$R(m)$',fontsize=25.0,fontname='normal')
    fig4.savefig('isolatedpotential.jpg',format='jpeg',dpi=250,bbox_inches='tight')
    show()
