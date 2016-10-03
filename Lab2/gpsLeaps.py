#!/usr/bin/env python
from numpy import *

	
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
	
	

	
if __name__ == "__main__":
	print 'Unit testing needs development'
