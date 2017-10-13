#This code aims to make a GIF to display two transimissions from COSMOS 2438 with a series of 0.5s snapshots.

import os
import aplpy
import glob
import ephem
import argparse
import requests
import imageio
import numpy as np
import getopt, datetime, pytz
import matplotlib.pyplot as plt
from matplotlib import animation 
import astropy.io.fits as pf
from astropy.time import Time


"""
New routines to convert to/from gps time
uses only the list of leap seconds and built-in python utilities

dates of leap seconds
http://en.wikipedia.org/wiki/Leap_second
Updated 2015-05-21
"""
#use in command line
parser = argparse.ArgumentParser(description='Make a movie for a satellite.')
#parser.add_argument("id",type=int,help="the number representing the time that the observation started")
#parser.add_argument("number",help="the NORAD_CAT_ID of the satellite ")
parser.add_argument("-n",help="the name of the movie that you want to make")
args=parser.parse_args()


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
#obs_id = args.id
"""obs_id = "{}".format(args.id)"""
obs_id = "1078397896"
snapshot = "2_%s_0055_I.fits" %(obs_id)
obs_time = gps_datetime(obs_id)+datetime.timedelta(seconds=28)
end_time = obs_time + datetime.timedelta(seconds=125)
delta_time = datetime.timedelta(seconds=0.5)
#download from website
a=obs_time.date()
b=datetime.datetime(obs_time.year,obs_time.month,obs_time.day+1).date()
os.system("curl -c cookies.txt -b cookies.txt https://www.space-track.org/ajaxauth/login -d 'identity=zsnxxx@sina.com&password=Darkness0924zsn'")
os.system('curl --limit-rate 100K --cookie cookies.txt https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/'+str(a)+'--'+str(b)+'/NORAD_CAT_ID/'+args.number+'/limit/1/format/tle > downloadTLE.txt')
"""args.number=37153 for STRELA 3,32955 for cosmos..."""

'''
'''
########TLE-DATA-READ-IN#############
def get_sat_ephem():
	# Read the ephemeris for the satellite
	#data=ephem.readtle(args.l1,args.l2,args.l3)
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
	

########POSITION-CALCULATION#############
def get_sat_radec(sat_ephem,timelist):
	l_ra = []
	l_dec = []
	# Observing from the MRO
	topo=ephem.Observer()
	topo.lon = '116:40:14.93'
	topo.lat = '-26:42:11.95'
	topo.elevation = 377.83
	#topo.lon = args.l
	#topo.lat = args.a
	#topo.elevation = float(args.e)
	# Work out the ra and dec for each time
	for time in timelist:	
		topo.date = time
		# Work out where the satellite is
		sat_ephem.compute(topo)
		# Get the RA and Dec in radians, convert to degrees and add to the ra/dec list
		l_ra.append(float(sat_ephem.a_ra) / ephem.degree)
		l_dec.append(float(sat_ephem.a_dec) / ephem.degree)
	return l_ra, l_dec


#read fits_list in
init=55
fitslist=[]
for i in range(125):
	name='2_1078397896_%04d_I.fits' %init
	fitslist.append(name)
	init+=1

#get seperate time range
timelist=[]
for f in fitslist:
	hdr = pf.getheader(f,0) 
	jd = hdr["JD"] 
	t=Time(jd, format='jd', scale='utc') 
	a=t.iso
	b=datetime.datetime.strptime(a, "%Y-%m-%d %H:%M:%S.%f")
	c=b+datetime.timedelta(seconds=2.5)
	timelist.append(c)

sat_ephem = get_sat_ephem()
l_ra, l_dec = get_sat_radec(sat_ephem,timelist)

linelist=np.array([l_ra,l_dec])

##plot in different frames#################3
j=0
imagelist=[]
for fits in fitslist:
	fig=aplpy.FITSFigure(fits)
	fig.tick_labels.set_xformat("hh:mm")
	fig.tick_labels.set_yformat("dd:mm")
	#add title
	fig.add_label(0.5, 1.05, "time: %s" %(timelist[j]), relative=True, size='xx-large', layer='title')
	# Set gray-scale
	fig.show_grayscale(invert=False, aspect='auto')
	# Show colorbar
	fig.add_colorbar()
	fig.colorbar.set_axis_label_text('Flux (Jy/beam)')
	fig.show_colorscale(stretch="linear",vmin=-8.0,vmax=8.0,cmap="Greys")
	# Add a grid
	fig.add_grid()
	fig.grid.show()
	fig.grid.set_xspacing(5.0)  # degrees
	fig.grid.set_yspacing(5.0)  # degrees
	#circles	
	#fig.show_circles(l_ra, l_dec,0.5, color='cyan',linewidth=1.0, alpha=0.5)
	fig.show_lines([linelist], color='blue', linewidth=2.0, alpha=0.5)
	fig.show_circles(l_ra[2], l_dec[2],0.2, color='cyan',linewidth=1.0, alpha=1.0)
	fig.show_circles(l_ra[122], l_dec[122],0.2, color='cyan',linewidth=1.0, alpha=1.0)
	fig.show_circles(l_ra[j], l_dec[j],0.2, color='r',linewidth=1.0, alpha=1.0)
	#plt.show()
	if j<10:
		fig.save('00'+str(j)+'.png')
		imagelist.append('00'+str(j)+'.png')
	elif j<100:
		fig.save('0'+str(j)+'.png')
		imagelist.append('0'+str(j)+'.png')
	else:
		fig.save(str(j)+'.png')
		imagelist.append(str(j)+'.png')
	fig.close()
        j+=1


############make movie###################
#ffmpeg
#os.system('ffmpeg -r 10 -i %03d.png -q:v 1 -r 10 '+args.n+'.mp4')
#os.system('ffmpeg -r 10 -i %03d.png -b:v 200M -r 10 '+args.n+'.mp4')

#animation
ilist=glob.glob('*.png')
imagelist=np.sort(ilist)
images = []
for filename in imagelist:
    images.append(imageio.imread(filename))
exportname = "fig3v.gif"
imageio.mimsave(exportname, images, 'GIF')
#kargs = { 'duration': 0.5 }
#imageio.mimsave(exportname, images, 'GIF', **kargs)

'''
#with imageio.get_writer('test.mp4', fps=10) as writer: ----need all the images have the same size
with imageio.get_writer('test.gif', mode='I') as writer:
    for imagename in imagelist:
        image = imageio.imread(imagename)
        writer.append_data(image)


if __name__ == "__main__":
	print imagelist
'''
