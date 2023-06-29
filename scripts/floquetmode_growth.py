import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#--------------#
# LOADING DATA #
#--------------#
param = np.loadtxt('./parameters.txt')
xx = np.loadtxt('./xx.txt')
Globalgamma = np.loadtxt('./globalgamma.txt')
Globalgamma = Globalgamma / 2.0

NumSteps = int(param[0])
dt = param[1]
length = int(param[2])
NumModes = int(param[3])

time = np.linspace(0.0,NumSteps*dt,NumSteps)
time_si = time / 1.0e4


#----------#
# ANALYSIS #
#----------#

#gamma_max = np.max( Globalgamma[1000:] )
gamma_min = np.min( Globalgamma[5000:] )

#gamma_max_loc = np.where( Globalgamma == gamma_max )
gamma_min_loc = np.where( Globalgamma == gamma_min )


jump = 20000
mask = jump / 2
NumPeaks = 6
peak_loc = np.zeros(NumPeaks,dtype=np.int)


#peak_loc[0] = gamma_max_loc[0]
peak_loc[0] = gamma_min_loc[0]

for pp in range(1,NumPeaks):

	# Next peak guess
	peak_loc[pp] = peak_loc[pp-1] + jump

	# Windowing
	window_min = peak_loc[pp] - mask
	window_max = peak_loc[pp] + mask

	# Picking the local maximum peak	
	data = Globalgamma[window_min:window_max]
	#peak_max = np.max( data )
	peak_min = np.min( data )
	#peak_max_loc = np.where( data == peak_max )
	peak_min_loc = np.where( data == peak_min )

	# Appending to existing list
	#peak_loc[pp] = peak_max_loc[0] + window_min			# indexing in data is local and starts at window_min
	peak_loc[pp] = peak_min_loc[0] + window_min

	#print 'Found peak:', pp+1,'/',NumPeaks


#------------------------#
# FITTING CURVE TO PEAKS #
#------------------------#

peak_num = np.arange(0,NumPeaks)
print 'Peak numbers:', peak_num

gamma_genmode = np.mean(Globalgamma[peak_loc[2]:peak_loc[NumPeaks-1]])
zeroshift = Globalgamma - gamma_genmode								# Creating new growth-rate array with 0 as GM_gamma
												# gamma_GM = avg<gamma_FM>
yy = zeroshift[peak_loc]									# zeroshift is a shift in gamma axis, NOT time
#tt = time_si[peak_loc - gamma_max_loc[0]]
tt = time_si[peak_loc - gamma_min_loc[0]]

highres_tt = np.linspace(0.0,tt[-1],100)
highres_peak_num = np.linspace(0.0,peak_num[-1],100)

def fit_func(t,c0,c1,c2):
	return c0 * np.exp( c1 * (t**c2) )
a, a_error = curve_fit( fit_func, tt, yy, maxfev=10000 )
p, p_error = curve_fit( fit_func, peak_num, yy, maxfev=10000 )

# Calculating time when growth < 0.1% GM growth

#GM_ref_gamma = 1.001 * gamma_genmode			# For decay
#GM_ref_gamma -= gamma_genmode
GM_ref_gamma = -1.001 * gamma_genmode			# For growth
GM_ref_gamma += gamma_genmode



#print 'Reference gamma:', GM_ref_gamma

t_GM = ( np.log(GM_ref_gamma/a[0]) / a[1] )**(1/a[2])				# decay
t_GM_period = ( np.log(GM_ref_gamma/p[0]) / p[1] )**(1/p[2])
#t_GM = ( -np.log(GM_ref_gamma/a[0]) / a[1] )**(1/a[2])				# growth
#t_GM_period = ( -np.log(GM_ref_gamma/p[0]) / p[1] )**(1/p[2])


#--------------#
# PRINT & PLOT #
#--------------#

print '---------------------------------------'
print 'GM formation time:', t_GM
print 'GM formation period:', t_GM_period
print 'Starting growth-rate:', Globalgamma[0]
print 'Average growth rate:', gamma_genmode
print 'Percentage deviation:', np.abs( (Globalgamma[0]-gamma_genmode) / Globalgamma[0] ) * 100.0
print 'Max growth a0:', a[0]
print 'Decay rate a1:', a[1]
print 'Time power a2:', a[2]
print '---------------------------------------'

#plt.plot(time,Globalgamma)
plt.figure(1)
plt.subplot(1,2,1)
#plt.plot( time_si-time_si[gamma_max_loc[0]], zeroshift )
plt.plot( time_si-time_si[gamma_min_loc[0]], zeroshift )
plt.plot( [time_si[0],time_si[-1]],[0.0,0.0] )
plt.ylim([-0.035,0.035])
#plt.plot( tt, yy )
plt.plot( [tt[0],tt[-1]],[GM_ref_gamma,GM_ref_gamma] )
plt.scatter(tt,yy)
#plt.plot( time_si, fit_func(time_si,a[0],a[1],a[2]) )
plt.plot( highres_tt, fit_func(highres_tt,a[0],a[1],a[2]) )

plt.subplot(1,2,2)
plt.scatter( peak_num,yy )
plt.plot( highres_peak_num, fit_func(highres_peak_num,p[0],p[1],p[2]) )

plt.show()



















