import numpy as np
import math
import matplotlib.pyplot as plt


# ---------------------- #
# FUNDAMENTAL PARAMETERS #
# ---------------------- #

a = 0.5
R = 5.0
shear = 2.0
qedge = 3.45
b_rs = 0.1
n0 = 80.0


r_beg = 0.3
r_end = 0.4


# ------------ #
# CALCULATIONS #
# ------------ #

length = 1000
r = np.linspace( r_beg,r_end,length )

q = qedge * np.power( r/a,shear )
ktheta = n0*q/r
dqdr = shear*q/r
nqp = n0*dqdr

m_min = math.ceil(np.min(q)*n0)
m_max = math.floor(np.max(q)*n0)
approx_m0 = (m_max+m_min)/2.0

m0 = np.round(approx_m0)
rs = a * np.power( m0/(qedge*n0),1.0/shear )

ktheta_rs = m0/rs
rhoi_rs = np.sqrt(b_rs)/ktheta_rs
rhoi = rhoi_rs * (R+r)/(R+rs)
b = np.power(ktheta*rhoi,2.0 )

xx = (r-rs)*nqp

# -------- #
# PLOTTING #
# -------- #

print 'm0 = ',m0
print 'rs = ',rs
print 'm_min = ',m_min
print 'm_max = ',m_max
print 'Num. modes = ',m_max-m_min

# NOTE
# r --> xx changes the axis scaling for simulated space

plt.figure(1)
plt.subplot(221)
plt.plot(xx,q)						# q profile
plt.ylabel('$q$',fontsize=15)
plt.subplot(222)
plt.plot(xx,dqdr)					# dq/dr
plt.ylabel('$dq/dr$',fontsize=15)
plt.subplot(223)
plt.plot(xx,b)						# b
plt.ylabel('$b$',fontsize=15)
plt.subplot(224)
plt.plot(xx,nqp,xx,shear*ktheta,'r--')			# nq'
plt.legend(('n*dq/dr','k_theta*shear'),'upper left')
plt.ylabel('$n dq/dr$',fontsize=15)

plt.show()

