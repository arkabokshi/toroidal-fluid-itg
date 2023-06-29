from numpy import *
from matplotlib.pyplot import *

xx = loadtxt('xx.txt')
profiles = loadtxt('profiles.txt')
params = loadtxt('parameters.txt')

length = params[2]
tau = params[6]


eta_s = profiles[0*length:1*length]
qprofile = profiles[1*length:2*length]
sigma = profiles[2*length:3*length]
ShearTerm = profiles[3*length:4*length]

eta_i = (eta_s*tau-1.0)/1.5

figure(1)
subplot(221)
plot(xx,eta_s,xx,eta_i)
legend(('$\eta_s$','$\eta_i$'),'best')
title('eta(x)')
subplot(222)
plot(xx,qprofile)
title('q-profile')
subplot(223)
plot(xx,sigma**2)
title('sigma^2')
subplot(224)
plot(xx,ShearTerm)
legend(('$\Omega^{p}$',),'upper left')
title('flow-shear(x)')
show()