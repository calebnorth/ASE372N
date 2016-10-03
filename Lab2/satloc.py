#!/usr/bin/env python
from numpy import *
from navConstants import *


def satloc(gpsWeek, gpsSec, sd):
#+==============================================================================+
# satloc : Return satellite location and velocity expressed in and relative to
# 			the ECEF reference frame.
#
#
# INPUTS
#
# gpsWeek ----- Week of true time at which SV location and velocity are
# 				desired.
#
# gpsSec ------ Seconds of week of true time at which SV location and velocity
# 				are desired.
#
# sd ---------- Ephemeris structure array for a single SV. Let ii be the
# 				numerical identifier (PRN identifier) for the SV whose location
# 				is sought. Then sd = satdata(ii). sd has the following
# 				fields:
#
# 	SVID - satellite number
# 	health - satellite health flag (0 = healthy; otherwise unhealthy)
# 	we - week of ephemeris epoch (GPS week, unambiguous)
# 	te - time of ephemeris epoch (GPS seconds of week)
# 	wc - week of clock epoch (GPS week)
# 	tc - time of clock epoch (GPS seconds of week)
# 	e - eccentricity (unitless)
# 	sqrta - sqrt of orbit semi-major axis (m^1/2)
# 	omega0 - argument of perigee (rad.)
# 	M0 - mean anomaly at epoch (rad.)
# 	L0 - longitude of ascending node at beginning of week (rad.)
# 	i0 - inclination angle at epoch (rad.)
# 	dOdt - longitude rate (rad / sec.)
# 	dn - mean motion difference (rad / sec.)
# 	didt - inclination rate (rad / sec.)
# 	Cuc - cosine correction to argument of perigee (rad.)
# 	Cus - sine correction to argument of perigee (rad.)
# 	Crc - cosine correction to orbital radius (m)
# 	Crs - sine correction to orbital radius (m)
# 	Cic - cosine correction to inclination (rad.)
# 	Cis - sine correction to inclination (rad.)
# 	af0 - 0th order satellite clock correction (s)
# 	af1 - 1st order satellite clock correction (s / s)
# 	af2 - 2nd order satellite clock correction (s / s^2)
# 	TGD - group delay time for the satellite (s)
#
#
# OUTPUTS
#
# rSvEcef ----- Location of SV at desired time expressed in the ECEF reference
# 				frame (m).
#
# vSvEcef ----- Velocity of SV at desired time relative to the ECEF reference
# 				frame and expressed in the ECEF reference frame (m/s). NOTE:
# 				vSvEcef is NOT inertial velocity, e.g., for geostationary SVs
# 				vSvEcef = 0.
#
#+------------------------------------------------------------------------------+
# References:
#
# 1: 'Global positioning system directorate systems engineering & integration
#		interface specification IS-GPS-200', IS-GPS-200H-003, Dec 9, 2015
#
# 2: Remondi, Benjamin W., "Computing satellite velocity using the broadcast 
#		ephemeris," GPS Solutions, vol. 8 no. 3, 2004.
#
#
#
# Author: Caleb North
#+==============================================================================+
	from gps2utc import gps2utc
	
	#~ if sd['health'] != 63:
		#~ print 'Warning: SV is unhealthy!'
		
	t = gps2utc(gpsWeek, gpsSec)
	toe = gps2utc(sd['we'], float(sd['te']))

	E = sd['e']
	A = sd['sqrta']**2 # Semi-major axis
	no = sqrt(GM/A**3) # Computed mean motion
	tk = (t - toe).total_seconds() # Time from ephemeris reference epoch

	n = no + sd['dn'] # Corrected mean motion
	Mk = sd['M0'] + n*tk # Mean anomaly
	while Mk < 0: Mk = Mk + 2*pie
	
	Ek = Mk
	delta = 1
	while delta > 1e-12: 
		Ek_new = Mk + E*sin(Ek) # Eccentric Anomaly: Mk = Ek - E*sin(Ek)
		#~ print Ek_new
		delta = abs(Ek-Ek_new)
		Ek = Ek_new
	
	vk = arctan2( sqrt(1-E**2)*sin(Ek), cos(Ek)-E ) # True Anomaly

		
	# Second Harmonic Perturbations
	Phik = vk + sd['omega0'] # Argument of Latitude
	duk = sd['Cus']*sin(2*Phik) + sd['Cuc']*cos(2*Phik) # Argument of Latitude Correction
	drk = sd['Crs']*sin(2*Phik) + sd['Crc']*cos(2*Phik) # Radius Correction
	dik = sd['Cis']*sin(2*Phik) + sd['Cic']*cos(2*Phik) # Inclination Correction
	
	uk = Phik + duk 				# Corrected Argument of Latitude
	rk = A*(1 - E*cos(Ek)) + drk 	# Corrected Radius
	ik = sd['i0'] + dik + sd['didt']*tk 		# Corrected Inclination
	
	xkp = rk*cos(uk) # 
	ykp = rk*sin(uk) # 
	
	Omegak = sd['L0'] + sd['dOdt']*tk - OmegaE*gpsSec # 
	
	xk = xkp*cos(Omegak) - ykp*cos(ik)*sin(Omegak)
	yk = xkp*sin(Omegak) + ykp*cos(ik)*cos(Omegak)
	zk = ykp*sin(ik)
	
	rSvEcef = array([[xk,yk,zk]]).T
	
	
	# Calculating Velocity #
	MkDot = n
	EkDot = MkDot/(1-E*cos(Ek))
	vkDot = sin(Ek)*EkDot*(1+E*cos(vk)) / ((1-E*cos(Ek))*sin(vk))
	PhikDot = vkDot
	dukDot = 2*(sd['Cus']*cos(2*Phik)-sd['Cuc']*sin(2*Phik))*PhikDot
	drkDot = 2*(sd['Crs']*cos(2*Phik)-sd['Crc']*sin(2*Phik))*PhikDot
	dikDot = 2*(sd['Cis']*cos(2*Phik)-sd['Cic']*sin(2*Phik))*PhikDot
	ukDot = PhikDot + dukDot
	rkDot = A*E*sin(Ek)*EkDot+drkDot
	ikDot = sd['didt']+ dikDot
	xkpDot = rkDot*cos(uk) - ykp*ukDot
	ykpDot = rkDot*sin(uk) + xkp*ukDot
	OmegakDot = sd['dOdt'] - OmegaE
	
	xkDot = xkpDot*cos(Omegak) - ykpDot*cos(ik)*sin(Omegak) + ykp*sin(ik)*sin(Omegak)*ikDot - yk*OmegakDot
	ykDot = xkpDot*sin(Omegak) + ykpDot*cos(ik)*cos(Omegak) - ykp*sin(ik)*cos(Omegak)*ikDot + xk*OmegakDot
	zkDot = ykpDot*sin(ik) + ykp*cos(ik)*ikDot
	
	vSvEcef = array([[xkDot,ykDot,zkDot]]).T
	
	
	return (rSvEcef, vSvEcef)
	
	
def eph2sd(b):
	# Constructs a SV dictionary based off the augmented Epheemris
	# Used for Units Tests
	
	sd= {}
	sd['SVID'] = b[0][0]
	sd['health'] = b[6][1]
	sd['we'] = b[5][2]
	sd['te'] = b[3][0]
	sd['wc'] = 1
	sd['tc'] = 1
	sd['e'] = b[2][1]
	sd['sqrta'] = b[2][3]
	sd['omega0'] = b[4][2]
	sd['M0'] = b[1][3]
	sd['L0'] = b[3][2]
	sd['i0'] = b[4][0]
	sd['dOdt'] = b[4][3]
	sd['dn'] = b[1][2]
	sd['didt'] = b[5][0]
	sd['Cuc'] = b[2][0]
	sd['Cus'] = b[2][2]
	sd['Crc'] = b[4][1]
	sd['Crs'] = b[1][1]
	sd['Cic'] = b[3][1]
	sd['Cis'] = b[3][3]
	sd['af0'] = 1
	sd['af1'] = 1
	sd['af2'] = 1
	sd['TGD'] = b[6][2]
	
	return sd
	
if __name__ == "__main__":
	print 'Performing Unit tests.'
	print ''
	print 'Test 1: Solving for Position and Velocity for SV1 at 1653 570957.338101566'
	
	
	sd = {
'SVID': 1,
'health': 63,
'we': 1653,
'te': 576000,
'wc': 1653,
'tc': 575999.999996647,
'e': 0.000301708350889,
'sqrta': 5153.66840553,
'omega0': 1.77176130104,
'M0': -2.03728834211,
'L0': -2.88543635184,
'i0': 0.960714536657,
'dOdt': -7.8396122656e-009,
'dn': 4.39839749662e-009,
'didt': -4.64305054455e-010,
'Cuc': -5.40167093277e-007,
'Cus': 9.78447496891e-006,
'Crc': 192.0625,
'Crs': -6.90625,
'Cic': 5.40167093277e-008,
'Cis': 1.49011611938e-008,
'af0': -2.90526077151e-006,
'af1': -1.5916157281e-012,
'af2': 0,
'TGD': 7.91624188423e-009}

	gpsWeek = 1653
	gpsSec = 570957.338101566
	(R,V) = satloc(gpsWeek, gpsSec, sd)
	print 'Computed Radius'
	print R
	print 'Computed Velocity'
	print V
	
	print 'Radius Error: '
	print array([[5734379.67991244,-18348853.9562115,-18337913.5840076]]).T - R
	
	print 'Velocity Error: '
	print array([[2075.74428570655,-1061.65685114262,1711.92358273592]]).T - V
	
	print ''
	#~ print 'Test 2: Solving for sd dictionary for SV1 at 1653 570957.338101566'
	#~ b = [[ 1 ,11 , 9 ,17 ,16 , 0 , 0.0,-2.905260771510e-06,-1.591615728100e-12, 0.000000000000e+00],
    #~ [1.050000000000e+02,-6.906250000000e+00, 4.398397496620e-09,-2.037288342110e+00],
   #~ [-5.401670932770e-07, 3.017083508890e-04, 9.784474968910e-06, 5.153668405530e+03],
   #~ [ 5.760000000000e+05, 5.401670932770e-08,-2.885436351840e+00, 1.490116119380e-08],
   #~ [ 9.607145366570e-01, 1.920625000000e+02, 1.771761301040e+00,-7.839612265600e-09],
   #~ [-4.643050544550e-10, 1.000000000000e+00, 1.653000000000e+03, 0.000000000000e+00],
   #~ [ 2.800000000000e+00, 6.300000000000e+01, 7.916241884230e-09, 1.050000000000e+02],
   #~ [ 5.748480000000e+05]]
	#~ sd1 = eph2sd(b)
	#~ print sd1
	#~ print sd
	#~ print ''
	
	print 'Test 3: Solving for Position and Velocity for SV2 at 12:00 Sep 1, 2013'
	from datetime import datetime as dt
	from utc2gps import utc2gps
	time0 = dt(2013,9,1,12,0,0)
	time1 = dt(2013,9,1,12,0,1)
	gpsWeek0, gpsSec0 = utc2gps(time0)
	gpsWeek1, gpsSec1 = utc2gps(time1)
	
	sv2 =  [[ 2, 13,  9,  1, 12,  0,  0.0, 4.492728039620e-04, 1.818989403550e-12, 0.000000000000e+00],
    [4.700000000000e+01, 7.218750000000e+00, 4.824843635730e-09,-1.442125029200e+00],
    [3.185123205190e-07, 1.254422636700e-02, 9.790062904360e-06, 5.153719497680e+03],
    [4.320000000000e+04, 2.309679985050e-07, 3.042975919520e+00,-2.328306436540e-07],
    [9.393591569710e-01, 1.811562500000e+02,-2.563763309250e+00,-8.004619189710e-09],
    [3.557291128330e-10, 1.000000000000e+00, 1.756000000000e+03, 0.000000000000e+00],
    [2.000000000000e+00, 0.000000000000e+00,-1.722946763040e-08, 4.700000000000e+01],
    [4.146000000000e+04]]
	sd2 = eph2sd(sv2)
	(R0,V0) = satloc(gpsWeek0, gpsSec0, sd2)
	(R1,V1) = satloc(gpsWeek1, gpsSec1, sd2)
	Vr = R1-R0
	print 'Computed Radius'
	print R0
	print 'Computed Velocity'
	print V0
	print ''
	print 'Computing solutions for one second later'
	print 'Computed Radius'
	print R1
	print 'Computed Velocity'
	print V1
	print ''
	print 'Estimated Veolocity based of Position differences:'
	print Vr
	print 'Error between Computed Velocity and Estimated Velocity:'
	print V0-Vr
	print ''
	print ''
	
	print 'Test 4: Solving for Position and Velocity for SV5 at 12:00 Sep 1, 2013'
	
	sv5 = [[ 5, 13 , 9 , 1, 12 , 0 , 0.0,-3.996668383480e-04, 0.000000000000e+00, 0.000000000000e+00],
    [3.800000000000e+01, 1.015625000000e+01, 5.044138884360e-09, 8.372052215750e-01],
    [6.239861249920e-07, 3.207757952620e-03, 5.116686224940e-06, 5.153687467570e+03],
    [4.320000000000e+04, 3.725290298460e-09,-2.174601614670e+00, 3.352761268620e-08],
    [9.481020352500e-01, 2.724062500000e+02, 2.569222944390e-01,-8.376420446150e-09],
    [5.821670923110e-11, 1.000000000000e+00, 1.756000000000e+03, 0.000000000000e+00],
    [2.000000000000e+00, 0.000000000000e+00,-9.778887033460e-09, 3.800000000000e+01],
    [3.600000000000e+04]]
	sd5 = eph2sd(sv5)
	(R0,V0) = satloc(gpsWeek0, gpsSec0, sd5)
	(R1,V1) = satloc(gpsWeek1, gpsSec1, sd5)
	Vr = R1-R0
	print 'Computed Radius'
	print R0
	print 'Computed Velocity'
	print V0
	print ''
	print 'Computing solutions for one second later'
	print 'Computed Radius'
	print R1
	print 'Computed Velocity'
	print V1
	print ''
	print 'Estimated Veolocity based of Position differences:'
	print Vr
	print 'Error between Computed Velocity and Estimated Velocity:'
	print V0-Vr
	print ''
	print ''
