#This code aims to make a movie to demonstrate the motions of any specified satellte during a given two-minute observation in a image automatically. The satellite and observation image can be specified by argparse. And the TLE data is going to be downloaded automatically.

#!/usr/bin/env python
import aplpy
import getopt, datetime, pytz
import numpy as np
import ephem
import matplotlib.pyplot as plt
import os
import argparse
import requests

"""
New routines to convert to/from gps time
uses only the list of leap seconds and built-in python utilities

dates of leap seconds
http://en.wikipedia.org/wiki/Leap_second
Updated 2015-05-21
"""
#use in command line
parser = argparse.ArgumentParser(description='Make a movie for a satellite.')
parser.add_argument("id",type=int,help="the number representing the time that the observation started")
parser.add_argument("number",help="the NORAD_CAT_ID of the satellite ")
parser.add_argument("-n",help="the name of the movie that you want to make")
args=parser.parse_args()
'''
#parser.add_argument("-l",help="the longitude of the observation site")
#parser.add_argument("-a",help="the altitude of the observation site")
#parser.add_argument("-e",help="the elevation of the observation site")
#parser.add_argument("-m",help="0 means minus'-'; 1 means positive '+'")
#parser.add_argument("-l1",help='TLE line one')
#parser.add_argument("-l2",help='TLE line two')
#parser.add_argument("-l3",help='TLE line three')
#print args.id,args.l,args.a,args.e,args.l1,args.l2,args.l3
'''


#datetime###########################
Leapseconds=[datetime.datetime(1972, 6, 30,23,59,59,0,pytz.utc),
	datetime.datetime(1972, 12, 31,23,59,59,0,pytz.utc),
	datetime.datetime(1973, 12, 31,23,59,59,0,pytz.utc),
	datetime.datetime(1974, 12, 31,23,59,59,0,pytz.utc),
	datetime.datetime(1975, 12, 31,23,59,59,0,pytz.utc),
	datetime.datetime(1976, 12, 31,23,59,59,0,pytz.utc),
	datetime.datetime(1977, 12, 31,23,59,59,0,pytz.utc),
	datetime.datetime(1978, 12, 31,23,59,59,0,pytz.utc),
	datetime.datetime(1979, 12, 31,23,59,59,0,pytz.utc),
	datetime.datetime(1981, 6, 30,23,59,59,0,pytz.utc),
	datetime.datetime(1982, 6, 30,23,59,59,0,pytz.utc),
	datetime.datetime(1983, 6, 30,23,59,59,0,pytz.utc),
	datetime.datetime(1985, 6, 30,23,59,59,0,pytz.utc),
	datetime.datetime(1987, 12, 31,23,59,59,0,pytz.utc),
	datetime.datetime(1989, 12, 31,23,59,59,0,pytz.utc),
	datetime.datetime(1990, 12, 31,23,59,59,0,pytz.utc),
	datetime.datetime(1992, 6, 30,23,59,59,0,pytz.utc),
	datetime.datetime(1993, 6, 30,23,59,59,0,pytz.utc),
	datetime.datetime(1994, 6, 30,23,59,59,0,pytz.utc),
	datetime.datetime(1995, 12, 31,23,59,59,0,pytz.utc),
	datetime.datetime(1997, 6, 30,23,59,59,0,pytz.utc),
	datetime.datetime(1998, 12, 31,23,59,59,0,pytz.utc),
	datetime.datetime(2005, 12, 31,23,59,59,0,pytz.utc),
	datetime.datetime(2008, 12, 31,23,59,59,0,pytz.utc),
	datetime.datetime(2012, 6, 30,23,59,59,0,pytz.utc),
	datetime.datetime(2015, 6, 30,23,59,59,0,pytz.utc)]
# this is when GPStime=0 is defined
GPSzero=datetime.datetime(1980,1,6,0,0,0,0,pytz.utc)
Leapinterval_Start=[datetime.datetime(1960,1,1,0,0,0,0,pytz.utc)]
Leapinterval_End=[]
SumLeapseconds=[0]
for ls in Leapseconds:
	Leapinterval_End.append(ls)
	SumLeapseconds.append(SumLeapseconds[-1]+1)
	Leapinterval_Start.append(ls)

Leapinterval_End.append(datetime.datetime(2018,1,1,0,0,0,0,pytz.utc))
# number of leap seconds that had already elapsed at GPStime=0
GPSzero_leapseconds=SumLeapseconds[np.where(np.array(Leapinterval_Start)<=GPSzero)[0][-1]]
Leapinterval_Start_GPS=[]
Leapinterval_End_GPS=[]
ls=0
for i in range(len(Leapinterval_Start)):
	dt=Leapinterval_Start[i]-GPSzero+datetime.timedelta(seconds=ls-GPSzero_leapseconds)
	Leapinterval_Start_GPS.append(dt.days*86400 + dt.seconds)
	dt=Leapinterval_End[i]-GPSzero+datetime.timedelta(seconds=ls-GPSzero_leapseconds)
	Leapinterval_End_GPS.append(dt.days*86400 + dt.seconds)
	ls=SumLeapseconds[i]        

##########TIME-CALCULATION############
def gps_datetime(gps_str):
	"""
	datetime=gps_datetime(gps)
	converts from gps seconds to datetime.datetime
	"""
	gps = int(gps_str)
	ls=SumLeapseconds[np.where(np.array(Leapinterval_Start_GPS)<gps)[0][-1]]
	dt=datetime.timedelta(seconds=gps-(ls-GPSzero_leapseconds))
	return GPSzero + dt

########TLE-DATA-Download#############
#FITS ID transfer
obs_id = args.id
'''obs_id = "{}".format(args.id)
obs_id = "1069427592"'''
snapshot = "2_%s_V.fits" %(obs_id)
obs_time = gps_datetime(obs_id)
'''obs_time = datetime.datetime(2013, 11, 25, 15, 14, 2, 0, pytz.utc)'''
end_time = obs_time + datetime.timedelta(seconds=112)
delta_time = datetime.timedelta(seconds=0.5)
#download from website
a=obs_time.date()
b=datetime.datetime(obs_time.year,obs_time.month,obs_time.day+1).date()
os.system("curl -c cookies.txt -b cookies.txt https://www.space-track.org/ajaxauth/login -d 'identity=zsnxxx@sina.com&password=Darkness0924zsn'")
os.system('curl --limit-rate 100K --cookie cookies.txt https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/'+str(a)+'--'+str(b)+'/NORAD_CAT_ID/'+args.number+'/limit/1/format/tle > downloadTLE.txt')
#args.number=37153 for STRELA 3


########TLE-DATA-READ-IN#############
def get_sat_ephem(obs_id):
	# Read the ephemeris for the satellite
	'''data=ephem.readtle(args.l1,args.l2,args.l3)'''
	f=open('downloadTLE.txt','r') #only read
	result=list()
	line0='name'
	for line in f.readlines():			  #read lines one by one
		line=line.strip()
		if not len(line) or line.startswith('#'): #skip the empty and comment lines
			continue
		result.append(line)
	data = ephem.readtle(line0,result[0], result[1])
	return data
	'''line1="1 37153U 10043B   13329.15715830 -.00000006  00000-0 -11935-3 0  9995"
	line2="2 37153 082.4609 055.5481 0006670 321.3736 038.6869 12.40699436145589"'''

########POSITION-CALCULATION#############
def get_sat_radec(sat_ephem, start_time, end_time, delta_time):
	l_ra = []
	l_dec = []
	# Observing from the MRO
	topo=ephem.Observer()
	topo.lon = '116:40:14.93'
	topo.lat = '-26:42:11.95'
	topo.elevation = 377.83
	'''#topo.lon = args.l
	#topo.lat = args.a
	#topo.elevation = float(args.e)'''
	# Work out the ra and dec for each time
	t = start_time
	while t < end_time:
		# Set up the new time
		topo.date = t
		# Work out where the satellite is
		sat_ephem.compute(topo)
		# Get the RA and Dec in radians, convert to degrees and add to the ra/dec list
		l_ra.append(float(sat_ephem.ra) / ephem.degree)
		l_dec.append(float(sat_ephem.dec) / ephem.degree)
		t += delta_time
	return l_ra, l_dec

sat_ephem = get_sat_ephem(obs_id)
l_ra, l_dec = get_sat_radec(sat_ephem, obs_time, end_time, delta_time)

linelist=np.array([l_ra,l_dec])


##Plot the image with the satellite trajectory one by one#################3
for i in range(len(l_ra)):
	fig=aplpy.FITSFigure(snapshot)
	fig.tick_labels.set_xformat("hh:mm")
	fig.tick_labels.set_yformat("dd:mm")
	# Set gray-scale
	fig.show_grayscale(invert=True, aspect='auto')
	# Show colorbar
	fig.add_colorbar()
	fig.colorbar.set_axis_label_text('Flux (Jy/beam)')
	# Add a grid
	fig.add_grid()
	fig.grid.show()
	fig.grid.set_xspacing(5.0)  # degrees
	fig.grid.set_yspacing(5.0)  # degrees
	#circles
	fig.show_circles(l_ra[i], l_dec[i],0.1, color='r',linewidth=1.0, alpha=0.5)
	#fig.show_lines([linelist], color='r', linewidth=2.0, alpha=0.5)
	#plt.show()
	if i<10:
		fig.save('00'+str(i)+'.png')
	elif i<100:
		fig.save('0'+str(i)+'.png')
	else:
		fig.save(str(i)+'.png')
	fig.close()

############make movie###################
os.system('ffmpeg -r 10 -i %03d.png -b:v 50M -r 10 '+args.n+'.mp4')

'''
if __name__ == "__main__":
	print linelist
'''
