#!/usr/bin/env python

from numpy import *
from pylab import *
from matplotlib import *
#~ from matplotlib import pyplot as plt
import aseTools as w

#(sod, lats, longs, hts, sats)

## Generating a datamap ##
data_files = [
'1a.txt',
'2a.txt',
'3a.txt',
'3b.txt',
'4a.txt',
'4b.txt'
]

out = 0

datamap = {}
for name in data_files:
	data = genfromtxt(name, delimiter=',')
	data[:,1:3] = data[:,1:3]*pi/180.
	datamap[name] = {}
	datamap[name]['lla'] = data

	print name
	
	## Solving for Mean lat, lon, alt, SVs (4.a, 4.b) 
	datamap[name]['llaMeans'] = mean(datamap[name]['lla'], axis=0)
	[SOD, LAT, LON, ALT, SVS] = mean(datamap[name]['lla'], axis=0)
	
	
	## Solving for std lat, lon, alt, SVs 
	datamap[name]['llaStats'] = std(datamap[name]['lla'], axis=0)
	
	
	## Solving for mean ECEF vector (4.d)
	pVec = w.lla2ecef(LAT,LON,ALT)
	datamap[name]['pVec'] = pVec
	#~ print 'ECEF vector is: '
	#~ print pVec
	if out ==1: print ''
	if out ==1: print ''
	print 'Runtime: '
	print (data[-1,0]-data[0,0])/60
	I = len(datamap[name]['lla'][:,0])
	datamap[name]['ecef'] = zeros((I,5))
	datamap[name]['delta'] = zeros((I,5))
	datamap[name]['enu'] = zeros((I,5))
	for i in range(I):
		sod = datamap[name]['lla'][i,0]
		lat = datamap[name]['lla'][i,1]
		lon = datamap[name]['lla'][i,2]
		alt = datamap[name]['lla'][i,3]
		svs = datamap[name]['lla'][i,4]
		pVec = w.lla2ecef(lat,lon,alt)
		R = w.ecef2enu(LAT, LON)
		diffVec = pVec - datamap[name]['pVec']
		datamap[name]['ecef'][i,:]  = column_stack((sod,pVec.T,svs))
		datamap[name]['delta'][i,:] = column_stack((sod,diffVec.T,svs))
		datamap[name]['enu'][i,:]   = column_stack((sod,dot(R,diffVec).T,svs))
	datamap[name]['ecefMeans'] = mean(datamap[name]['ecef'], axis=0)
	datamap[name]['ecefStats'] =  std(datamap[name]['ecef'], axis=0)
	datamap[name]['enuStats']  =  std(datamap[name]['enu'],  axis=0)
	[sod, sigE, sigN, sigU, sigSV] = datamap[name]['enuStats']
	datamap[name]['errors'] = array([sqrt(sigE**2+sigN**2)/1.2, sqrt(sigE**2+sigN**2+sigU**2)/1.3])
	
	#~ print 'LLA means: '
	#~ print datamap[name]['llaMeans']
	#~ print 'LLA sigmas: '
	#~ print datamap[name]['llaStats']
	print 'ECEF means: '
	print datamap[name]['ecefMeans']
	#~ print 'ECEF sigmas: '
	#~ print datamap[name]['ecefStats']
	#~ print 'ENU sigmas: ' 
	#~ print datamap[name]['enuStats']
	#~ print 'Errors: '
	#~ print datamap[name]['errors']
	
	print ''
	print ''
	
## Estimating distance bewteen WAAS 
dVec = datamap['1a.txt']['pVec']-datamap['2a.txt']['pVec']
print 'Distance between WAAS is: ' + str(linalg.norm(dVec))

## Estimating distance bewteen markers (5)
dVec = datamap['3a.txt']['pVec']-datamap['4a.txt']['pVec']
print 'Distance between markers is: ' + str(linalg.norm(dVec))


## Computing residual heights and deviation (7.b)
fname = '3a.txt'
print datamap[fname]['enuStats']
print datamap[fname]['errors']




## Plotting ##

one = '1a'
two = '2a'

# WAAS Comparison
x1 = datamap[one+'.txt']['enu'][:,1]
y1 = datamap[one+'.txt']['enu'][:,2]
x2 = datamap[two+'.txt']['enu'][:,1]
y2 = datamap[two+'.txt']['enu'][:,2]
scatter(x1,y1,color='red',marker='+',s=60,label='WAAS Enabled')
scatter(x2,y2,color='blue',marker='x',s=60,label='WAAS Disabled')
xlabel(r'$\Delta$ East (m)')
ylabel(r'$\Delta$ North (m)')
gca().set_aspect('equal', adjustable='box')
xlim(-10,10)
ylim(-10,10)
gca().set_xticks(arange(-10,11,2))
gca().set_yticks(arange(-10,11,2))
gca().legend(loc='best')
grid()
#~ suptitle('ENU Residuals')
savefig('ENU-residual-WAAS')
close()


# Mustang Comparsions

one = '3a'
two = '3b'


x1 = datamap[one+'.txt']['lla'][:,2]-datamap[one+'.txt']['llaMeans'][2]
y1 = datamap[one+'.txt']['lla'][:,1]-datamap[one+'.txt']['llaMeans'][1]
x2 = datamap[two+'.txt']['lla'][:,2]-datamap[two+'.txt']['llaMeans'][2]
y2 = datamap[two+'.txt']['lla'][:,1]-datamap[two+'.txt']['llaMeans'][1]
scatter(x1*1e9,y1*1e9,color='red',marker='+',s=60,label='Test 1')
scatter(x2*1e9,y2*1e9,color='blue',marker='x',s=60,label='Test 2')
xlabel(r'$\Delta$ Longitude (nrad)')
ylabel(r'$\Delta$ Latitude (nrad)')
gca().set_aspect('equal', adjustable='box')
xlim(-1600,1600)
ylim(-1600,1600)
gca().set_xticks(arange(-1600,1601,400))
gca().set_yticks(arange(-1600,1601,400))
gca().legend(loc='best')
grid()
#~ suptitle('LLA Residuals')
savefig('LLA-residuals-Mustang')
close()



x1 = datamap[one+'.txt']['enu'][:,1]
y1 = datamap[one+'.txt']['enu'][:,2]
x2 = datamap[two+'.txt']['enu'][:,1]
y2 = datamap[two+'.txt']['enu'][:,2]
scatter(x1,y1,color='red',marker='+',s=60,label='Test 1')
scatter(x2,y2,color='blue',marker='x',s=60,label='Test 2')
xlabel(r'$\Delta$ East (m)')
ylabel(r'$\Delta$ North (m)')
gca().set_aspect('equal', adjustable='box')
xlim(-10,10)
ylim(-10,10)
gca().set_xticks(arange(-10,11,2))
gca().set_yticks(arange(-10,11,2))
gca().legend(loc='best')
grid()
#~ suptitle('ENU Residuals')
savefig('ENU-residuals-Mustang')
close()

## Plotting Tests ##

one = '1a'
two = '2a'


x1 = datamap[one+'.txt']['lla'][:,2]-datamap[one+'.txt']['llaMeans'][2]
y1 = datamap[one+'.txt']['lla'][:,1]-datamap[one+'.txt']['llaMeans'][1]
x2 = datamap[two+'.txt']['lla'][:,2]-datamap[two+'.txt']['llaMeans'][2]
y2 = datamap[two+'.txt']['lla'][:,1]-datamap[two+'.txt']['llaMeans'][1]
scatter(x1*1e9,y1*1e9,color='red',marker='+',s=60,label=one)
scatter(x2*1e9,y2*1e9,color='blue',marker='x',s=60,label=two)
xlabel(r'$\Delta$ Longitude ($n rad$)')
ylabel(r'$\Delta$ Latitude ($n rad$)')
gca().set_aspect('equal', adjustable='box')
xlim(-1500,1500)
ylim(-1500,1500)
#~ gca().set_xticks(arange(-1500,1501,500))
#~ gca().set_yticks(arange(-1500,1501,500))
gca().legend(loc='best')
grid()
suptitle('Delta LLA')
#~ show()
close()



x1 = datamap[one+'.txt']['enu'][:,1]
y1 = datamap[one+'.txt']['enu'][:,2]
x2 = datamap[two+'.txt']['enu'][:,1]
y2 = datamap[two+'.txt']['enu'][:,2]
scatter(x1,y1,color='red',marker='+',s=40,label=one)
scatter(x2,y2,color='blue',marker='x',s=40,label=two)
xlabel(r'$\Delta$ East ($m$)')
ylabel(r'$\Delta$ North ($m$)')
gca().set_aspect('equal', adjustable='box')
xlim(-10,10)
ylim(-10,10)
gca().set_xticks(arange(-10,11,2))
gca().set_yticks(arange(-10,11,2))
gca().legend(loc='best')
grid()
suptitle('ENU Difference Vectors')
#~ show()
#~ close()



