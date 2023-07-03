import numpy as np
import matplotlib.pyplot as plt
from numpy import *
from matplotlib.pyplot import *
from math import *
import matplotlib.animation as animation

# ------------ #
# Loading data #
# ------------ #
parameters = np.loadtxt('parameters.txt')
videoparam = np.loadtxt('videoparam.txt')

xx = loadtxt('./xx.txt')
growthrate = loadtxt('./growthrate.txt')
frequency  = loadtxt('./frequency.txt')
#PolPotOmega = loadtxt('./polpotomega.txt')
#PolPotGamma = loadtxt('./polpotgamma.txt')
Globalgamma = loadtxt('./globalgamma.txt')
Globalgamma = Globalgamma / 2.0
gammaE_t = loadtxt('./gammaE_t.txt')

NumSteps = int(parameters[0])
dt = parameters[1]
Length = int(parameters[2])
NumModes = int(parameters[3])
InitialMode = int(parameters[4])
FinalMode = int(parameters[5])

NumFiles = int(videoparam[0])
DelPrint = int(videoparam[1])


time = linspace(0.0,NumSteps*dt,NumSteps)
time = time / 1.0e4
mode_gamma = reshape(growthrate,(NumModes,NumSteps))
mode_omega = reshape(frequency, (NumModes,NumSteps))

RealData = np.zeros((NumModes,Length,NumFiles)) 
ImagData = np.zeros((NumModes,Length,NumFiles))

for tt in range(NumFiles):
  tstamp = DelPrint * tt
  TempData = np.loadtxt(str(100000+tstamp)+'.txt')
  RealData[:,:,tt] = np.reshape(TempData,(NumModes,Length)) 
  TempData = np.loadtxt(str(200000+tstamp)+'.txt')
  ImagData[:,:,tt] = np.reshape(TempData,(NumModes,Length)) 
 
  

# --------------------------- #
# Routine for plotting phi(x) #
# --------------------------- #
print( '------------------------------------' )
print( 'Reconstructing poloidal potential...' )
print( '------------------------------------' )

NumTheta = 720

Potential = np.zeros((NumTheta,Length,NumFiles))
RealLocalPotential = np.zeros(Length)
ImagLocalPotential = np.zeros(Length)


dplot = 1

for tt in range(0,NumFiles,dplot):
  print('Timeslice:',tt+1,'/',NumFiles)

  for jj in range(NumTheta):

	  Theta = float(jj)/(NumTheta-1)*2.0*pi

	  RealLocalPotential = 0.0
	  ImagLocalPotential = 0.0

	  for m in range(InitialMode,FinalMode+1):
	
		  RealLocalPotential += RealData[m-InitialMode,:,tt]*cos(m*Theta) + ImagData[m-InitialMode,:,tt]*sin(m*Theta)
		  ImagLocalPotential += ImagData[m-InitialMode,:,tt]*cos(m*Theta) - RealData[m-InitialMode,:,tt]*sin(m*Theta)
		  
	  Potential[jj,:,tt] = np.sqrt( RealLocalPotential**2 + ImagLocalPotential**2 )
	  #Potential[jj,:,tt] = RealLocalPotential


# --------------------------------- #
# Routine for 2D Poloidal animation #
# --------------------------------- #

InnerRadius = 0.25
OuterRadius = 0.5


rr = np.linspace(InnerRadius,OuterRadius,Length) 
pa = np.linspace(0.0,2.0*pi,NumTheta)
yy = np.zeros((NumTheta,Length))
xx = np.zeros((NumTheta,Length))
for jj in range(NumTheta):
      yy[jj,:] = rr * sin(pa[jj])
      xx[jj,:] = rr * cos(pa[jj])


fig = plt.figure(figsize=(16,5),dpi=1000)
#plt.hold(False)
ims = []

for tt in range(0,NumFiles,dplot):

    
    levels = np.linspace(np.min(Potential[:,:,tt]),np.max(Potential[:,:,tt]),50)
    timepoint = DelPrint * tt * dt / 1.0e4

    yll = Globalgamma[-1] * 0.90
    yul = Globalgamma[-1] * 1.10
    
    plt.subplot( 1,3,2 )
    
    plt.contourf(xx,yy,Potential[:,:,tt],levels=levels,cmap=cm.jet)
    #plt.xticks(fontsize=20.0,fontname='serif')
    #plt.yticks(fontsize=20.0,fontname='serif')
    #plt.imshow(Potential[:,:,tt].transpose())
    plt.ylabel('Z (arb. units)',fontsize=15.0,fontname='serif')
    plt.xlabel('R (arb. units)',fontsize=15.0,fontname='serif')
    plt.title('$|\phi_1(x,\Theta)|$',fontsize=20.0,fontname='serif')
    
    plt.subplot( 1,3,3  )
    for i in range(20,60):
      plot(time,mode_gamma[i,:],linewidth=0.5)
    plot(time,Globalgamma,'y--',linewidth=2.5)
    #plot([0.0,time[-1]],[0.743984083872,0.743984083872],'k--',linewidth=1.5)
    plot([timepoint,timepoint],[yll,yul],'k',linewidth=1.5)
#    ylabel('$\gamma$',fontsize=15.0)
    xlabel('Time (ms)',fontsize=15.0,fontname='serif',visible=True)
    title('$\gamma / \omega_e^*$',fontsize=20.0,fontname='serif')
    xticks(np.linspace(0.0,np.max(time),3), fontsize=15.0,fontname='serif',visible=True)
    yticks(fontsize=15.0,fontname='serif')
    ylim([yll,yul])
    
    plt.subplot( 1,3,1 )
    plot(time,gammaE_t,'b',linewidth=2.0)
    plot([timepoint,timepoint],[-0.004,0.004],'k',linewidth=1.5)
    #ylabel('$\gamma_E$',fontsize=15.0,visible='False')
    title('$\gamma_E = d\Omega/dq $',fontsize=20.0,fontname='serif')
    xlabel('Time (ms)',fontsize=15.0,fontname='serif')
    xticks(np.linspace(0.0,np.max(time),3),fontsize=15.0,fontname='serif')
    yticks(fontsize=15.0,fontname='serif')
    ylim([-0.004,0.004])

    #plt.tight_layout()    

    fig.savefig('./images/'+str(tt)+'.jpg',format='jpeg',dpi=250,bbox_inches='tight')
    
    print( 'Saving:',tt+1,'/',NumFiles )
    plt.clf()
    



