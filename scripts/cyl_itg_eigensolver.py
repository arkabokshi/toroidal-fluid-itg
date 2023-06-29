from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import math
from numpy.polynomial.hermite import hermval

#-----------------------#
# SIMULATION PARAMETERS #
#-----------------------#


hermite_mode =  15
tolerance    =  1.0e-6

# SIMULATION DOMAIN

x_beg  =  -45.0
x_end  =  +45.0
PtDen  =  10
length =  int(x_end-x_beg)*PtDen + 1

xx = np.linspace( x_beg,x_end,length )

ic = (length-1)/2
dx = xx[2]-xx[1]

fwd = +1	# Direction for shooting routine
bwd = -1


# DECLARING FIELDS

beta = complex(0.0,0.0)

state_old   =   np.zeros( 2,dtype=np.complex128 )
state_new   =   np.zeros( 2,dtype=np.complex128 )

qprofile  =  np.zeros( length )
eta       =  np.zeros( length )

#--------------------#
# PHYSICAL CONSTANTS #
#--------------------#

m0 = 90
n0 = 50
delM = 0

MinRad = 0.5
MajRad = 5.0

krhoi =  np.sqrt( 0.1 )
shear =  2.0
tau   =  1.0
eps   =  0.03088

eta0 = 5.0
etag = 0.0

b = krhoi**2
c = b * tau

rs       =    MinRad * ( float(m0)/(3.45*float(n0)) )**(1/shear)
ktheta   =    float(m0) / rs
nqp      =    ktheta * shear


eta       +=   eta0 * (1.0 - etag*(xx/nqp)**2)
qprofile  +=   float(m0) / float(n0)
#qprofile  +=   3.45 * ( (xx/nqp + rs)/MinRad )**shear

sigma     =   eps / ( np.sqrt(c)*qprofile )

#---------------------------#
# CALCULATING THE POTENTIAL #
#---------------------------#

def Background_Potential( Omega ):


	# FLR term
	flr_term = -1.0 * c
	
	# ION-SOUND term
	sound_term = ( sigma/Omega * (delM-xx) )**2
	
	# EIGENVALUE term
	eigen_term = -1.0 * ( Omega - 1.0 ) / ( Omega + eta )
	
	# V(x) in ( d2 /d2x + V(x) ) P = 0
	# NOTE : check for any normalisation before the ( d2 /dx2 ) term
	Potential = ( sound_term + eigen_term + flr_term ) / ( c*(shear**2) )

	# returning them allows to update the global variables
	return Potential



#-----------------------------------#
# CALCULATING SLOPE FOR RK4 ROUTINE #
#-----------------------------------#

def df( step, ff, Vr, Vi ):

	dv = np.zeros(2,dtype=np.complex128)

	PotV = complex( Vr(step),Vi(step) )

	dv[0]  =  ff[1] * ( -PotV )			# dv = [ -V*phi, dphi ]
	dv[1]  =  ff[0]  

	return dv



#-------------#
# RK4 ROUTINE #
#-------------#

def shoot( jj, u0, way, Vr, Vi ):			# u0 = [ dphi, phi ]

	dh = dx * way

  	k1  =  df(  xx[jj]        , u0            ,  Vr , Vi  )
	k2  =  df(  xx[jj] + dh/2 , u0 + dh*k1/2  ,  Vr , Vi  )
  	k3  =  df(  xx[jj] + dh/2 , u0 + dh*k2/2  ,  Vr , Vi  )
  	k4  =  df(  xx[jj] + dh   , u0 + dh*k3    ,  Vr , Vi  )

	u1  =  u0 + (dh/6.0)*( k1 + 2.0*k2 + 2.0*k3 + k4 )

	return u1



#--------------------#
# SHOOTING ALGORITHM #
#--------------------#

def f( Omega_guess ):
	
	global phi, phi_evn, phi_odd
	global Vx

	#-----------------------------------#
	# 1.) CREATING A FUNCTION TO RETURN #
	#     POTENTIAL AT ANY SPECIFIED XX #
	#-----------------------------------#	

	# NOTE: it is possible to combine (Vr,Vi) into V in interp1d --> not in documentation

	Vx = Background_Potential( Omega_guess )
	Vr = interp1d( xx, Vx.real, kind='linear', bounds_error=False )
	Vi = interp1d( xx, Vx.imag, kind='linear', bounds_error=False ) 


	# DEFINING FIELDS
	phi         =   np.zeros( length,dtype=np.complex128 )
	dphi        =   np.zeros( length,dtype=np.complex128 )
	phi_odd     =   np.zeros( length,dtype=np.complex128 )
	phi_evn     =   np.zeros( length,dtype=np.complex128 )
	dphi_odd    =   np.zeros( length,dtype=np.complex128 )
	dphi_evn    =   np.zeros( length,dtype=np.complex128 )


	#--------------------------------#
	# 2.) SETTING CONDITIONS AT i=ic #
	#--------------------------------#
	
	# EVEN Function
	phi_evn[ic]   =  complex( 1.0,1.0 ) 
	dphi_evn[ic]  =  complex( 0.0,0.0 )
	
	# ODD Function
	phi_odd[ic]   =  complex( 0.0,0.0 ) 
	dphi_odd[ic]  =  complex( 1.0,1.0 )
	

	#-----------------#
	# 3.) RIGHT SHOOT #
	#-----------------#
	
	# NOTE: state = [ dphi,phi ]


	# 1.) EVEN FUNCTION
	
	for ii in range( ic,length-1 ):			# range(1,5) = [1,2,3,4]
	
		state_old[0]  =  dphi_evn[ii]
		state_old[1]  =  phi_evn[ii]
	
		state_new = shoot( ii,state_old,fwd,Vr,Vi )
	
		dphi_evn[ii+1]  =  state_new[0]
		phi_evn[ii+1]   =  state_new[1] 
		
	
	# 2.) ODD FUNCTION

	for ii in range( ic,length-1 ):
	
		state_old[0]  =  dphi_odd[ii]
		state_old[1]  =  phi_odd[ii]
	
		state_new = shoot( ii,state_old,fwd,Vr,Vi )
	
		dphi_odd[ii+1]  =  state_new[0]	
		phi_odd[ii+1]   =  state_new[1]   

	
	#----------------#
	# 4.) LEFT SHOOT #
	#----------------#
	
	# NOTE: While shooting left, dx -> -dx
	# This has been adjusted in the RK4 routine
	# This also accounts for the change in sign of slope at i = ic
	

	# 1.) EVEN FUNCTION

	for ii in range( ic,0,-1 ):
	
		state_old[0]  =  dphi_evn[ii]
		state_old[1]  =  phi_evn[ii]
	
		state_new = shoot( ii,state_old,bwd,Vr,Vi )
		
		dphi_evn[ii-1]  =  state_new[0]
		phi_evn[ii-1]   =  state_new[1] 
	
	
	# 2.) ODD FUNCTION
	
	for ii in range( ic,0,-1 ):
	
		state_old[0]  =  dphi_odd[ii]
		state_old[1]  =  phi_odd[ii]
	
		state_new = shoot( ii,state_old,bwd,Vr,Vi )
		
		dphi_odd[ii-1]  =  state_new[0] 
		phi_odd[ii-1]   =  state_new[1] 	
		

	#---------------------------------#
	# 5.) FINDING THE TOTAL POTENTIAL #
	#---------------------------------#

	# NOTE:
	# check why the odd and even needs to be divided up
	# depending on even/odd

	if hermite_mode%2 == 0:
		beta = -phi_evn[-1] / phi_odd[-1]
		phi = beta*phi_odd + phi_evn
	else:
		beta = -phi_odd[-1] / phi_evn[-1]
		phi = beta*phi_evn + phi_odd

	
	return phi[0]


#---------------------#
# ANALYTICAL SOLUTION #
#---------------------#

# 1.) EIGEN-VALUE

def analytic_solve( k ):

	
	mode = float( 2*k + 1 )

	ao = complex( 1.0 + c , 0.0 )
	bo = complex( c*eta[ic] - 1.0 , mode * sigma[ic] * shear * np.sqrt(c) )
	co = complex( 0.0 , mode * eta[ic] * sigma[ic] * shear * np.sqrt(c) )

	omega_sub = ( -bo - np.sqrt(bo**2 - 4.0*ao*co) ) / ( 2.0*ao )

	return omega_sub


# 2.) EIGEN-FUNCTION

def analytic_eigenfunction( HermiteMode,Omega ):

	# Defining fields
	GaussianTerm = np.zeros( length,dtype=np.complex128 )
	HermiteTerm  = np.zeros( length,dtype=np.complex128 )

	# Defining coef as required by hermval
	# coef = [ c0,c1,c2... ] where ci multiplies the ith Hermite mode
	ci = np.zeros( 51,dtype=float )
	ci[HermiteMode] = 1


	delta   =  sigma[ic] / ( Omega * shear * np.sqrt(c) )
	
	alpha2  =  complex( 0.0 , delta )
	alpha   =  np.sqrt( alpha2 )

	GaussianTerm  =  np.exp( -1.0*alpha2*(xx**2)/2.0 )
	HermiteTerm   =  hermval( -1.0*alpha*xx,ci )

	Analytical_Phi = GaussianTerm * HermiteTerm

	return Analytical_Phi


#------------------#
# FINDING NEW MODE #
#------------------#

def new_mode( Omega ):

	RHS = -Omega / complex(0.0,sigma[ic]*shear*np.sqrt(c)) * ( (Omega-1.0)/(Omega+eta[ic]) + c )
	kk = ( RHS-1.0 ) / 2.0	

	return kk


#------------------------#
# MAIN ROUTINE           #
# ROOT FINDING ALGORITHM #
#------------------------#

if __name__=="__main__":

	dif=1.0e-3
	mul=1.0+dif

	# Give it an informed 
	# hermite_mode defined at the start
	x0 = analytic_solve( hermite_mode )

	x1 = x0*mul
	x2 = x1*mul
	e0 = f(x0)
	e1 = f(x1)
	e2 = f(x2)

	iter_max  = 100

	print '---------------'
	print 'Finding root...'
	print '---------------'

	for ii in range(iter_max):

		# P(x) = a(x-x2)^2 + b(x-x2) + c
		
		# defining coefficients
		cc   =   e2
		bb   =   ( ((x0-x2)**2)*(e1-e2) - ((x1-x2)**2)*(e0-e2) ) / ( (x0-x2)*(x1-x2)*(x0-x1) ) 
		aa   =   ( (x1-x2)*(e0-e2) - (x0-x2)*(e1-e2) ) / ( (x0-x2)*(x1-x2)*(x0-x1) )
			
		# finding new root xr
		root_1  =  -2.0*cc / ( bb + np.sqrt(bb**2 - 4.0*aa*cc) )
		root_2  =  -2.0*cc / ( bb - np.sqrt(bb**2 - 4.0*aa*cc) )
	
		# root = xr-x2
		# Pick the root closest to x2
		if np.abs(root_1) <= np.abs(root_2):
			xr = root_1 + x2
		else:
			xr = root_2 + x2

		# re-define guesses
		x0 = x1
		x1 = x2
		x2 = xr

		e0 = e1
		e1 = e2
		e2 = f(xr)

		print xr, abs(e2)

		# checking tolerance and exiting
		if abs(e2) <= tolerance:
			print ''
			print '-----------------'
			print 'Converged to root'
			print '-----------------'
			print ''
			print 'Hermite mode to investigate:', hermite_mode
			print 'Analytical Omega :', analytic_solve( hermite_mode )
			print 'Eigenvalue Omega :', xr
			print 'Hermite mode found:', np.abs(new_mode( xr ))
			print ''
			print '-----'
			print 'NOTE:'
			print 'If Hermite mode found NOT integer, mode could be non-physical'
			print 'If Hermite mode found CLOSE to integer, try with increased domain size'
			print '-----'
			print ''
			break

	if not abs(e2) <= tolerance:
		print '----------------------------------------------'
		print 'Maximum iterations exceeded: increase iter_max'
		print '----------------------------------------------'
	
#----------#
# PLOTTING #
#----------#

plt.figure(1)
# plt.subplot(3,1,1)
# plt.plot(xx,phi_odd.real,xx,phi_odd.imag)
# plt.legend(['Real ODD_PHI','Imag ODD_PHI'],'best')
# plt.ylim([-2.0,2.0])
# plt.subplot(3,1,2)
# plt.plot(xx,phi_evn.real,xx,phi_evn.imag)
# plt.legend(['Real EVN_PHI','Imag EVN_PHI'],'best')
# plt.ylim([-2.0,2.0])
# plt.subplot(3,1,3)
plt.plot(xx,phi.real,xx,phi.imag,xx,np.abs(phi))
plt.legend(['Real PHI','Imag PHI','abs PHI'],'best')
plt.ylim([-np.abs(phi).max(),np.abs(phi).max()])

plt.figure(2)
plt.plot( xx,analytic_eigenfunction(hermite_mode,xr) )
plt.legend(['Analytic eigen-function',],'best')
plt.show()
	

	
# Comparing fit to data

#Vr = interp1d( xx, Vx.real, kind='linear', bounds_error=False )
#Vi = interp1d( xx, Vx.imag, kind='linear', bounds_error=False ) 

#plt.figure(3)
#plt.subplot(3,1,1)
#plt.plot(xx,Vx.real,linestyle='none',marker='o',markerfacecolor='white',markeredgecolor='blue')
#plt.plot(xx,Vr(xx),'r')
#plt.legend(['Real Vx','Spline'],'best')
#plt.subplot(3,1,2)
#plt.plot(xx,Vx.imag,linestyle='none',marker='o',markerfacecolor='white',markeredgecolor='blue')
#plt.plot(xx,Vi(xx),'r')
#plt.legend(['Imag Vx','Spline'],'best')
#plt.subplot(3,1,3)
#plt.plot(xx,np.sqrt( Vx.real**2+Vx.imag**2 ),linestyle='none',marker='o',markerfacecolor='white',markeredgecolor='blue')
#plt.legend(['abs(Vx)',],'best')

#plt.show()




