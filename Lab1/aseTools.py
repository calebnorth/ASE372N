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
	
	
