import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import *
from matplotlib.axes import Axes


folder = '/home/arkaprava/Project/ITGCylinder/InitialValue_vs_EigenValue/david_param/'

Omega_eps = np.loadtxt( folder+'omegaeps.txt' )
Gamma_eps = np.loadtxt( folder+'gammaeps.txt' )
eps_space = np.loadtxt( folder+'epsspace.txt' )
most_unstable_mode = np.loadtxt( folder+'mostunstablemode.txt' )

print '---------------------------------------'
print '(1) Continuity breaking down at start  '
print 'due to the limit set on MaxOrder       '
print '(2) For ODD modes, EVEN initial guess  '
print '==> convergence to wrong eigenvalue    '
print '---------------------------------------'


plt.figure(1)
plt.subplot(3,1,1)
plt.plot(eps_space,Gamma_eps)
plt.subplot(3,1,2)
plt.plot(eps_space,Omega_eps)
plt.subplot(3,1,3)
plt.plot(eps_space,most_unstable_mode)
#plt.ylim([0,10])
plt.show()


# These values are for tau = 10.0, eta = 1.0, eps = [0.1,1.0]
#NumEps = [0.1,0.2,0.25,0.35,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
#ErrorOmega = [0.7719791396,0.0055739337,0.0097689397,0.0059516492,0.0082779196,
#0.0070676209,0.0079541571,0.0091237362,0.0102405527,0.0112967517,0.0078319938]
#ErrorGamma = [1.6302795805,0.1609515227,0.0035552508,0.025873928,0.0251508976,
#0.0177754993,0.0222275912,0.0284647071,0.0355042122,0.0433763332,0.0164607768]


# PLOT 0  
#fig0,ax0 = plt.subplots(1,figsize=(8,3),sharex=False,sharey=False)
#fig0
#ax0.plot(eps_space,most_unstable_mode,'r',linewidth=3.0)
#ax0.set_ylabel('Most unstable mode',fontsize=20.0,fontname='serif',labelpad=10.0)
#ax0.set_xlabel('$\epsilon_n$',fontsize=25.0,fontname='serif',labelpad=-0.5)
#plt.yticks(fontsize=20.0,fontname='serif')
#plt.xticks(fontsize=20.0,fontname='serif')
#ax0dup = ax0.twinx()
#ax0dup.plot(eps_space,most_unstable_mode,'r',linewidth=1.0)
#ax0dup.set_ylabel('$|\omega_{num}/\omega_{ana} - 1|*100$',color='b',fontsize=25.0,fontname='serif',visible=False)
#plt.yticks(fontsize=20.0,fontname='serif')
#plt.xticks(fontsize=20.0,fontname='serif')

# PLOT 1
#fig1,ax1 = plt.subplots(1,figsize=(8,3),sharex=False)
#ax1.plot(eps_space,Omega_eps,'k-',linewidth=3.0)
#ax1.set_ylabel('$\omega_{ana}$',color='k',fontsize=25.0,fontname='serif')
#ax1.set_xlabel('$\epsilon_n$',fontsize=25.0,fontname='serif',labelpad=-0.5)
#plt.yticks(fontsize=20.0,fontname='serif')
#plt.xticks(fontsize=20.0,fontname='serif')
#for tl in ax1.get_yticklabels():
#    tl.set_color('k')
#ax1dup = ax1.twinx()
#ax1dup.plot(NumEps,ErrorOmega,linewidth=1.5,linestyle='dashed',color='blue',\
#marker='s',markerfacecolor='blue',markeredgecolor='blue',markersize=7.5,markeredgewidth=1.0)
#ax1dup.set_ylabel('$|\omega_{num}/\omega_{ana} - 1|*100$',color='b',fontsize=25.0,fontname='serif')
#ax1dup.set_yscale('log')
#for tl in ax1dup.get_yticklabels():
#    tl.set_color('b')
#plt.yticks(fontsize=20.0,fontname='serif')

# PLOT 2
#fig2,ax2 = plt.subplots(1,figsize=(8,3),sharex=False)
#ax2.plot(eps_space,Gamma_eps,'k-',linewidth=3.0)
#ax2.set_xlabel('$\epsilon_n$',fontsize=25.0,fontname='serif',labelpad=-0.5)
#ax2.set_ylabel('$\gamma_{ana}$',color='k',fontsize=25.0,fontname='serif')
#plt.yticks(fontsize=20.0,fontname='serif')
#plt.xticks(fontsize=20.0,fontname='serif')
#for tl in ax2.get_yticklabels():
#    tl.set_color('k')
#ax2dup = ax2.twinx()
#ax2dup.plot(NumEps,ErrorGamma,linewidth=1.5,linestyle='dashed',color='blue',\
#marker='s',markerfacecolor='blue',markeredgecolor='blue',markersize=7.5,markeredgewidth=1.0)
#ax2dup.set_ylabel('$|\gamma_{num}/\gamma_{ana} - 1|*100$',color='b',fontsize=25.0,fontname='serif')
#ax2dup.set_yscale('log')
#for tl in ax2dup.get_yticklabels():
#    tl.set_color('b')
#plt.yticks(fontsize=20.0,fontname='serif')
#plt.xticks(fontsize=20.0,fontname='serif')


#plt.show()
#fig0.savefig('MostunstableMode.jpg',format='jpeg',dpi=250)#,bbox_inches='tight')
#fig1.savefig('MostunstableOmega.jpg',format='jpeg',dpi=250,bbox_inches='tight')
#fig2.savefig('MostunstableGamma.jpg',format='jpeg',dpi=250,bbox_inches='tight')
  
  
