#!/usr/bin/env python

from numpy import *


def enu2ecef(lat, lon):
#+==============================================================================+
# ecef2enu : Generate the rotation matrix used to express a vector written in 
#			 local east, north, up (ENU) coordinates as a vector written in ECEF 
#			 coordinates at the position defined by geodetic latitude 
#			 and longitude.
# INPUTS 
# 
# lat (phi)   ----- geodetic latitude in radians 
# 
# lon (lamda) ----- longitude in radians 
# 
# 
# OUTPUTS 
# 
# R ------- 3-by-3 rotation matrix that maps a vector v_enu expressed in the 
# local east, north, up (vertical) reference frame to a vector v_ecef expressed 
# in the ECEF reference frame as follows: 
# v_ecef = R*v_enu. 
# 
#+------------------------------------------------------------------------------+ 
# References: 
#	1:  Fundamentals of inertial Navigation, Satellite-based Positioning and their Integration
#		Noureldin, A; Karamat, T.B.; Gregory, J.
#		2013, XVIII, 314p. Hardcover
#		ISBN: 978-3-642-30465-1
#		Site: http://www.springer.com/978-3-642-30465-1 ### Needs clarification ##
# 
# Author: Caleb North
#+==============================================================================+


	R = array([
		[-sin(lon), -sin(lat)*cos(lon), cos(lat)*cos(lon)],
		[ cos(lon), -sin(lat)*sin(lon), cos(lat)*sin(lon)],
		[ 		0., 		  cos(lat), 	 	 sin(lat)]
		])

	return (R)
	
	
def ecef2enu(lat, lon):
#+==============================================================================+
# ecef2enu : Generate the rotation matrix used to express a vector written in 
#			 ECEF coordinates as a vector written in local east, north, up 
# 			 (ENU) coordinates at the position defined by geodetic latitude 
#			 and longitude.
# INPUTS 
# 
# lat (phi)   ----- geodetic latitude in radians 
# 
# lon (lamda) ----- longitude in radians 
# 
# 
# OUTPUTS 
# 
# R ------- 3-by-3 rotation matrix that maps a vector v_ecef expressed in the 
# ECEF reference frame to a vector v_enu expressed in the local 
# east, north, up (vertical) reference frame as follows: 
# v_enu = R*v_ecef. 
# 
#+------------------------------------------------------------------------------+ 
# References: 
#	1:  Fundamentals of inertial Navigation, Satellite-based Positioning and their Integration
#		Noureldin, A; Karamat, T.B.; Gregory, J.
#		2013, XVIII, 314p. Hardcover
#		ISBN: 978-3-642-30465-1
#		Site: http://www.springer.com/978-3-642-30465-1 ### Needs clarification ##
# 
# Author: Caleb North
#+==============================================================================+


	R = array([
		[		 -sin(lon), 		  cos(lon), 	   0.],
		[-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)],
		[ cos(lat)*cos(lon),  cos(lat)*sin(lon), sin(lat)]
		])



	return (R)


#~ print ecef2enu(0,0)



def ecef2lla(pVec):
#+==============================================================================+
# ecef2lla : Convert from a position vector in the Earth-centered, Earth-fixed 
#			 (ECEF) reference frame to latitude, longitude, and altitude 
#			 (geodetic with respect to the WGS-84 ellipsoid). 
# 
# 
# INPUTS 
# 
# pVec ---- 3-by-1 position coordinate vector in the ECEF reference frame, 
# in meters. 
# 
# OUTPUTS 
# 
# lat ----- latitude in radians 
# 
# lon ----- longitude in radians 
# 
# alt ----- altitude (height) in meters above the ellipsoid 
#
#	1:  Fundamentals of inertial Navigation, Satellite-based Positioning and their Integration
#
# Caveats: Only good for abs(lat) > .0001 ## more clarification
#
# 
# 
#+------------------------------------------------------------------------------+ 
# References: 
#	1:  Fundamentals of inertial Navigation, Satellite-based Positioning and their Integration
#		Noureldin, A; Karamat, T.B.; Gregory, J.
#		2013, XVIII, 314p. Hardcover
#		ISBN: 978-3-642-30465-1
#		Site: http://www.springer.com/978-3-642-30465-1 ### Needs clarification ##
# 
# Author: Caleb North
#+==============================================================================+

	x = pVec[0,0]
	y = pVec[1,0]
	z = pVec[2,0]
	
	
	# WGS84 Values #
	a = 6378137.0			## Semimajor axis (equatorial radius)
	f = 1/298.257223563		## Flattening
	b = a*(1-f)				## Semiminor axis
	e = sqrt(f*(2-f))		## Eccentricity
	E = sqrt(a**2/b**2 -1)	## 
	
	p = sqrt(x**2 + y**2)
	theta = arctan(z*a/(p*b))
	
	# Using closed form solution #
	lat = arctan((z+E**2*b*sin(theta)**3)/(p-e**2*a*cos(theta)**3))
	lon = arctan2(y,x)
	N = a/sqrt(1-e**2*sin(lat)**2)	## Normal Radius
	alt = p/cos(lat) - N
	
	## accounting for numerical instabilities ##
	if p  < 1 and z > 0: alt =  z - b
	if p  < 1 and z < 0: alt = -z - b
	if y == 0 and x > 0: lon = 0
	if y == 0 and x < 0: lon = pi
	if p == 0 and z > 0: lat =  pi/2; lon = 0; alt =  z-b
	if p == 0 and z < 0: lat = -pi/2; lon = 0; alt = -z-b
	
	
	return (lat, lon, alt)




def ecef2lla_test():
#+==============================================================================+
# ecef2lla : Convert from a position vector in the Earth-centered, Earth-fixed 
#			 (ECEF) reference frame to latitude, longitude, and altitude 
#			 (geodetic with respect to the WGS-84 ellipsoid). 
# 
# 
# INPUTS 
# 
# pVec ---- 3-by-1 position coordinate vector in the ECEF reference frame, 
# in meters. 
# 
# OUTPUTS 
# 
# lat ----- latitude in radians 
# 
# lon ----- longitude in radians 
# 
# alt ----- altitude (height) in meters above the ellipsoid 
#
#
# Caveats: Only good for abs(lat) > .0001 ## more clarification
#
# 
# 
#+------------------------------------------------------------------------------+ 
# References: 
#	1:  Fundamentals of inertial Navigation, Satellite-based Positioning and their Integration
#		Noureldin, A; Karamat, T.B.; Gregory, J.
#		2013, XVIII, 314p. Hardcover
#		ISBN: 978-3-642-30465-1
#		Site: http://www.springer.com/978-3-642-30465-1 ### Needs clarification ##
# 
# Author: Caleb North
#+==============================================================================+

	x = linspace(-100,100,10000)
	y = 0.00
	z = 6378137.0	
	
	
	# WGS84 Values #
	a = 6378137.0			## Semimajor axis (equatorial radius)
	f = 1/298.257223563		## Flattening
	b = a*(1-f)				## Semiminor axis
	e = sqrt(f*(2-f))		## Eccentricity
	E = sqrt(a**2/b**2 -1)	## 
	
	p = sqrt(x**2 + y**2)
	theta = arctan(z*a/(p*b))
	#~ lat = arctan((z+E**2*b*sin(theta)**3)/(p-e**2*a*cos(theta)**3))
	#~ lon = 2*arctan(y/(x + p))

	# Using closed form solution #
	lat = arctan2((z+E**2*b*sin(theta)**3)/(p-e**2*a*cos(theta)**3))
	lon = arctan2(y,x)

	N = a/sqrt(1-e**2*sin(lat)**2)	## Normal Radius
	alt = p/cos(lat) - N
	#~ print lon
	for i in range(len(x)):
		## accounting for numerical instabilities ##
		if p[i]  < 10 and z > 0: alt[i] =  z - b
		if p[i]  < 10 and z < 0: alt[i] = -z - b
		#~ if y == 0 and x[i] > 0: lon[i] = 0
		#~ if y == 0 and x[i] < 0: lon[i] = pi
		if p[i] == 0 and z > 0: lat[i] =  pi/2; lon[i] = 0; alt[i] =  z-b
		if p[i] == 0 and z < 0: lat[i] = -pi/2; lon[i] = 0; alt[i] = -z-b
	
	
	return (lat, lon, alt,x)





def lla2ecef(lat,lon,alt):
#+==============================================================================+
# lla2ecef : 	Convert from latitude, longitude, and altitude (geodetic with 
# 				respect to the WGS-84 ellipsoid) to a position vector in the 
# 				Earth-centered, Earth-fixed (ECEF) reference frame. 
# 
# 
# INPUTS 
# 
# lat ----- latitude in radians 
# 
# lon ----- longitude in radians 
# 
# alt ----- altitude (height) in meters above the ellipsoid 
# 
# 
# OUTPUTS 
# 
# pVec ---- 3-by-1 position coordinate vector in the ECEF reference frame, 
# 			in meters. 
# 
#+------------------------------------------------------------------------------+ 
# References: 
#	1:  Fundamentals of inertial Navigation, Satellite-based Positioning and their Integration
#		Noureldin, A; Karamat, T.B.; Gregory, J.
#		2013, XVIII, 314p. Hardcover
#		ISBN: 978-3-642-30465-1
#		Site: http://www.springer.com/978-3-642-30465-1 ### Needs clarification ##
# 
# Author: Caleb North
#+==============================================================================+
	
	# WGS84 Values #
	a = 6378137.0			## Semimajor axis (equatorial radius)
	f = 1/298.257223563		## Flattening
	b = a*(1-f)				## Semiminor axis
	e = sqrt(f*(2-f))		## Eccentricity
	E = sqrt(a**2/b**2 -1)	## 

	N = a/sqrt(1 - e**2*sin(lat)**2)		## Normal radius
	
	pVec = array([
		[(N+alt)*cos(lat)*cos(lon)],
		[(N+alt)*cos(lat)*sin(lon)],
		[(N - N*e**2 + alt)*sin(lat)]
		])
		

	return (pVec)


def Test():
#+==============================================================================+
# Test :	Tests the various definitions found in this script.
# 
# 
# INPUTS 
# 
# N/A
# 
# 
# OUTPUTS 
# 
# Multiple
# 
#+------------------------------------------------------------------------------+ 
# References: 
#	N/A
# 
# Author: Caleb North
#+==============================================================================+
	a = 6378137.0
	import matplotlib.pyplot as plt
	
	a,b,c,x = ecef2lla_test()
	
	plt.plot(x,c)
	
	plt.show()
	
	
	
	
def gpsLeaps():
#+==============================================================================+
# gpsLeaps : Returns a list of dates when leaps seconds were introduced.
#
#
# INPUTS
#
#
# OUTPUTS
#
# Leaps ----------- A list of datetime structured dates which correspond to 
#					leap second insertions.
#
#+------------------------------------------------------------------------------+
# References:
#
#	1: https://confluence.qps.nl/display/KBE/UTC+to+GPS+Time+Correction
#
#	2: https://en.wikipedia.org/wiki/Leap_second
#
#
# Author: Caleb North
#+==============================================================================+
	
	from datetime import datetime as dt
	
	LeapTable= [dt(1981,7,1),
				dt(1982,7,1), 
				dt(1983,7,1),
				dt(1985,7,1), 
				dt(1988,1,1),
				dt(1990,1,1), 
				dt(1991,1,1),
				dt(1992,7,1), 
				dt(1993,7,1),
				dt(1994,7,1), 
				dt(1996,1,1), 
				dt(1997,7,1),
				dt(1999,1,1), 
				dt(2006,1,1),
				dt(2009,1,1), 
				dt(2012,7,1),
				dt(2015,7,1), 
				dt(2017,1,1)]
			
	return LeapTable
	
	

	
def utc2gps(utcTime):
#+==============================================================================+
# utc2gps : Convert UTC time to GPS time expressed in GPS week and GPS
# second of week.
#
#
# INPUTS
#
# n --------------- The UTC time and date expressed as a Matlab datenum. Use
#					the Matlab function datenum() to generate n; use datestr()
#					to render n in a standard format.
#
#
# OUTPUTS
#
# gpsWeek --------- The unambiguous GPS week number where zero corresponds to
# 					midnight on the evening of 5 January/morning of 6 January,
# 					1980. By unambiguous is meant the full week count since
# 					the 1980 reference time (no rollover at 1024 weeks).
#
# gpsSec ---------- The GPS time of week expressed as GPS seconds from midnight
# 					on Saturday.
#
#+------------------------------------------------------------------------------+
# References:
#
#
# Author: Caleb North
#+==============================================================================+
	
	from datetime import datetime as dt
	from datetime import timedelta
	import sys
	
	## Verify that the time is valid
	gpsStart = dt(1980,1,6)
	if utcTime < gpsStart:
		sys.exit('ERROR: The time in question occurs before GPS time begins.')
	
	leaps = 0
	LeapTable = gpsLeaps()
	for leap in LeapTable:
		if utcTime > leap:
			leaps +=1
	
	gpsTime = utcTime + timedelta(seconds=leaps)
	gpsWeek = (gpsTime - gpsStart).days/7
	gpsSec  = (gpsTime - gpsStart - timedelta(days=gpsWeek*7)).total_seconds()
	
	
	## Warn if the LeapTable may be outdated
	if utcTime > LeapTable[-1]:
		print ''
		print 'Warning: The last known leap second occured: ', LeapTable[-1]
		print 'It would be wise to verify there have not been additional leaps.'
		print ''
	
	return gpsWeek, gpsSec
	
	
	
	
	







def gps2utc(gpsWeek, gpsSec):
#+==============================================================================+
# gps2utc : Convert GPS time expressed in GPS week and GPS second of week to
# UTC.
#
#
# INPUTS
#
# gpsWeek --------- The unambiguous GPS week number where zero corresponds to
# 					midnight on the evening of 5 January/morning of 6 January,
# 					1980. By unambiguous is meant the full week count since
# 					the 1980 reference time (no rollover at 1024 weeks).
#
# gpsSec ---------- The GPS time of week expressed as GPS seconds from midnight
# 					on Saturday.
#
#
# OUTPUTS
#
# n --------------- The UTC time and date expressed as a Matlab datenum. Use
# 					the Matlab function datestr() to read n in a standard
# 					format.
#
#+------------------------------------------------------------------------------+
# References:
#
#
# Author: Caleb North
#+==============================================================================+
	
	from datetime import datetime as dt
	from datetime import timedelta
	import sys

	gpsStart = dt(1980,1,6)
	gpsTime = gpsStart + timedelta(days=gpsWeek*7) + timedelta(seconds=gpsSec)
	
	## Verify that the time is valid
	if gpsTime < gpsStart:
		sys.exit('ERROR: The time in question occurs before GPS time begins.')
	
	utcTime = gpsTime
	leaps = 0
	LeapTable = gpsLeaps()
	for leap in LeapTable:
		if utcTime > leap:
			leaps += 1
			utcTime = gpsTime - timedelta(seconds=leaps)
	
	## Warn if the LeapTable may be outdated
	if utcTime > LeapTable[-1]:
		print ''
		print 'Warning: The last known leap second occured: ', LeapTable[-1]
		print 'It would be wise to verify there have not been additional leaps.'
		print ''
	
	
	return utcTime
	
	
	
	
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
#
#
#
# Author:
#+==============================================================================+
	
	
	return 1
	
	
	
	
	
	




