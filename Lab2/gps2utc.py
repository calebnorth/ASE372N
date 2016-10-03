#!/usr/bin/env python
from numpy import *


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
# utcTime --------- The UTC time and date expressed as a Python datetime. 
#
#+------------------------------------------------------------------------------+
# References:
#
# 1: 'Global positioning system directorate systems engineering & integration
#		interface specification IS-GPS-200', IS-GPS-200H-003, Dec 9, 2015
#
# Author: Caleb North
#+==============================================================================+
	
	from datetime import datetime as dt
	from datetime import timedelta
	import sys
	from gpsLeaps import gpsLeaps
	
	# Find uncorrected UTC time
	gpsStart = dt(1980,1,6)
	gpsTime = gpsStart + timedelta(days=gpsWeek*7) + timedelta(seconds=gpsSec)
	
	## Verify that the time is valid
	if gpsTime < gpsStart:
		sys.exit('ERROR: The time in question occurs before GPS time begins.')
	
	## Correct for number of Leaps
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
	
	

	
if __name__ == "__main__":
	print 'Unit testing needs development'
