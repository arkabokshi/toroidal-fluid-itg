from numpy import *
import numpy as np
from pylab import *

ModeSpace  = loadtxt('modespace.txt')
GammaSpace = loadtxt('gammaspace.txt')
OmegaSpace = loadtxt('omegaspace.txt')
eta_space = loadtxt('etaspace.txt')
tau_space = loadtxt('tauspace.txt')

etamin = min(eta_space)
etamax = max(eta_space)
taumin = min(tau_space)
taumax = max(tau_space)

Data = transpose(ModeSpace)
LogData = np.log10(Data + 1.0)

X,Y = meshgrid(eta_space,tau_space)

figure(1)
im = imshow(Data,cmap=cm.jet,aspect='auto',origin='lower',extent=(etamin,etamax,taumin,taumax))
#clf()

# Contour lines
Levels = linspace(0,2,2)
cs = contour(X,Y,Data,Levels,colors='gray')
clabel(cs,colors='white',inline=True,fmt='%1.1f',fontsize=15)
print where(ModeSpace==0)

#cf = contourf(X,Y,LogData)
colorbar(im)#ticks = Levels)#spacing='proportional',
xlabel('$\eta_s$',fontsize=20)
ylabel('$T_e/T_i$',fontsize=20)


figure(2)
subplot(121)
#gs = contourf(X,Y,GammaSpace,cmap=cm.jet)
gs = imshow(transpose(GammaSpace),cmap=cm.jet,aspect='auto',origin='lower',extent=(etamin,etamax,taumin,taumax))
colorbar(gs)
title('Gamma for the fastest growing mode')
xlabel('$\eta_s$',fontsize=20)
ylabel('$T_e/T_i$',fontsize=20)
subplot(122)
#os = contourf(X,Y,OmegaSpace,cmap=cm.jet)
os = imshow(transpose(OmegaSpace),cmap=cm.jet,aspect='auto',origin='lower',extent=(etamin,etamax,taumin,taumax))
title('Frequency of the fastest growing mode')
xlabel('$\eta_s$',fontsize=20)
ylabel('$T_e/T_i$',fontsize=20)
colorbar(os)




show()

# Using imshow()	:	http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.imshow
# Chosing colormaps 	:	http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
# Using contour()	:	http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.contour
