#!/usr/bin/env python

from numpy import *
def utc2gps(utcTime):
#+==============================================================================+
# utc2gps : Convert UTC time to GPS time expressed in GPS week and GPS
# second of week.
#
#
# INPUTS
#
# utcTime --------- The UTC time and date expressed as a Python datetime. 
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
# 1: 'Global positioning system directorate systems engineering & integration
#		interface specification IS-GPS-200', IS-GPS-200H-003, Dec 9, 2015
#
#
# Author: Caleb North
#+==============================================================================+
	
	from datetime import datetime as dt
	from datetime import timedelta
	import sys
	from gpsLeaps import gpsLeaps

	
	## Verify that the time is valid
	gpsStart = dt(1980,1,6)
	if utcTime < gpsStart:
		sys.exit('ERROR: The time in question occurs before GPS time begins.')
	
	## Find number of Leaps
	leaps = 0
	LeapTable = gpsLeaps()
	for leap in LeapTable:
		if utcTime > leap:
			leaps +=1
	
	## Find GPS time
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
	
	
if __name__ == "__main__":
	print 'Unit testing needs development'
