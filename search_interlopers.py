#This code aims to conduct an automatical search for all the interlopers within a given image and make a movie to show satellites' motions during two minutes.

import os
import numpy as np
import astropy.io.fits as pf
import getopt, datetime, pytz
import ephem
import aplpy


#datetime
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

#TIME-CALCULATION
def gps_datetime(gps_str):
	"""
	datetime=gps_datetime(gps)
	converts from gps seconds to datetime.datetime
	"""
	gps = int(gps_str)
	ls=SumLeapseconds[np.where(np.array(Leapinterval_Start_GPS)<gps)[0][-1]]
	dt=datetime.timedelta(seconds=gps-(ls-GPSzero_leapseconds))
	return GPSzero + dt
#calculate ephem
obs_id = "1069427592"
snapshot = "/home/trina/Desktop/2_%s_V.fits" %(obs_id)
obs_time = gps_datetime(obs_id)
'''obs_time = datetime.datetime(2013, 11, 25, 15, 14, 2, 0, pytz.utc)'''
end_time = obs_time + datetime.timedelta(seconds=112)
delta_time = datetime.timedelta(seconds=0.5)
a=datetime.datetime(obs_time.year,obs_time.month,obs_time.day-2).date()
b=datetime.datetime(obs_time.year,obs_time.month,obs_time.day+2).date()
#TLE-DATA-Download#
'''
#DOWNLOADS
os.system("curl -c cookies.txt -b cookies.txt https://www.space-track.org/ajaxauth/login -d 'identity=zsnxxx@sina.com&password=Darkness0924zsn'")
os.system('curl --limit-rate 100K --cookie cookies.txt https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/'+str(a)+'--'+str(b)+'/orderby/NORAD_CAT_ID/format/3le > /home/trina/Desktop/catalog1.txt')#specific epoch download
#os.system('curl --limit-rate 100K --cookie cookies.txt https://www.space-track.org/basicspacedata/query/class/tle_latest/ORDINAL/1/EPOCH/%3Enow-30/orderby/NORAD_CAT_ID/format/3le > /home/trina/Desktop/catalog.txt')-------full catalog

#https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/2013-10-26--2013-12-21/orderby/NORAD_CAT_ID/format/3le -----repeated
#https://www.space-track.org/basicspacedata/query/class/tle_latest/ORDINAL/1/EPOCH/2013-10-26--2013-12-21/orderby/NORAD_CAT_ID/format/3le -----wrong list without STRELA3
'''

#read_in
f=open('/home/trina/Desktop/catalog1.txt','r') #only read
r=list()
for line in f.readlines():			  #read lines one by one
	line=line.strip()
	if not len(line) or line.startswith('#'): #skip the empty and comment lines
		continue
	r.append(line)
l=len(r)/3
tlelist=np.reshape(r,(l,3))

#delete repeated objects
def dele(tlelist):
	unique,copy=[],[]
	for tle in tlelist:
		temp=tle[1][2:7]
		if not temp in unique:
			copy.append(tle)
			unique.append(temp)
	return copy

finaltlelist=dele(tlelist) #delete repeated ones

#POSITION-CALCULATION#
#def pos_list function
def get_sat_radec(sat_ephem, start_time, end_time, delta_time):
	l_ra = []
	l_dec = []
	# Observing from the MRO
	topo=ephem.Observer()
	topo.lon = '116:40:14.93'
	topo.lat = '-26:42:11.95'
	topo.elevation = 377.83
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


#get cordlist
#loop
cordlist,namelist=[],[]
for tle in finaltlelist:
	sat_ephem = ephem.readtle(tle[0],tle[1],tle[2])
	#l_ra, l_dec = get_sat_radec(sat_ephem, obs_time, end_time, delta_time)
	ob_pos= get_sat_radec(sat_ephem, obs_time, end_time, delta_time)
	cordlist.append(ob_pos)
	namelist.append(tle[0])
	#ralist.append(l_ra)
	#declist.append(l_dec)
#print np.shape(cordlist),namelist

#center_cord
fig=snapshot
hdr = pf.getheader(fig,0) 
ra_c=hdr["RA_PH"]
dec_c=hdr["DEC_PH"] 

#smallest distance
def angular_dist(r1,d1,r2,d2):
	r1=np.radians(r1)
  	r2=np.radians(r2)
  	d1=np.radians(d1)
  	d2=np.radians(d2)
  	a=np.sin(np.abs((d1-d2)/2))**2
  	b=np.cos(d1)*np.cos(d2)*np.sin(np.abs((r1-r2)/2))**2
  	d=2*np.arcsin(np.sqrt(a+b))
  	dd=np.degrees(d)
  	return dd

def find_closest(tup,x,y):
  	min_d=angular_dist(tup[0][0],tup[1][0],x,y)
	minra,mindec=tup[0][0],tup[1][0]
  	for i in range(len(tup[0])):
    		dis=angular_dist(tup[0][i],tup[1][i],x,y)
    		if dis<min_d:
      			min_d=dis
			minra,mindec=tup[0][i],tup[1][i]
  	return min_d,minra,mindec

#judge smallest distance
def judgment(cordlist,max_d):
  	select_namelist,select_tlelist,select_pos_list=[],[],[]
	i=0
    	for ob_pos in cordlist:
  		min_d,minra,mindec=find_closest(ob_pos,ra_c,dec_c)
  	  	if min_d<max_d:
  	    		#t=(i,namelist[i],min_d,minra,mindec)
			t=(i,namelist[i])
			select_namelist.append(t)
			select_tlelist.append(finaltlelist[i])
			select_pos_list.append(ob_pos)
		i+=1
  	return select_namelist,select_tlelist,select_pos_list

select_name,select_tle,select_pos=judgment(cordlist,10)

#get one right position to indicate the existence
def get_stationary_radec(sat_ephem):
	# Observing from the MRO
	topo=ephem.Observer()
	topo.lon = '116:40:14.93'
	topo.lat = '-26:42:11.95'
	topo.elevation = 377.83
	t = '2013-11-25 15:14:04'
	topo.date = t
	sat_ephem.compute(topo)
	# Get the RA and Dec in radians, convert to degrees and add to the ra/dec list
	ra=(float(sat_ephem.ra) / ephem.degree)
	dec=(float(sat_ephem.dec) / ephem.degree)
	return ra, dec
#get ...list
sta_pos_list=[]
sta_ra,sta_dec=[],[]
for tle in select_tle:
	sat_ephem = ephem.readtle(tle[0],tle[1],tle[2])
	pos= get_stationary_radec(sat_ephem)		
	sta_pos_list.append(pos)
	sta_ra.append(pos[0])
	sta_dec.append(pos[1])


###########plot###########
'''
##plot stationary one
fig=aplpy.FITSFigure('/home/trina/Desktop/strela3.fits')
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
fig.show_circles(sta_ra, sta_dec, 0.1, color='r',linewidth=0.5, alpha=1.0)
#fig.show_lines([linelist], color='r', linewidth=2.0, alpha=0.5)
#plt.show()
fig.save('/home/trina/Desktop/interlopers1.png')
fig.close()

##plot stationary path
fig=aplpy.FITSFigure('/home/trina/Desktop/strela3.fits')
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
for pos in select_pos:
	fig.show_circles(pos[0], pos[1], 0.1, color='r',linewidth=1.0, alpha=1.0)
#fig.show_lines([linelist], color='r', linewidth=2.0, alpha=0.5)
#plt.show()
fig.save('/home/trina/Desktop/path.png')
fig.close()
'''

#make move
clist=['pink','deeppink','mediumvioletred','olive','lime','green','salmon','red','orange','yellow','cyan','blue']
for i in range(len(select_pos[0][0])):
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
	j=0
	for pos in select_pos:
		fig.show_circles(pos[0][i], pos[1][i],0.5, color=clist[j],linewidth=1.0, alpha=0.5)
		j+=1
	if i<10:
		fig.save('00'+str(i)+'.png')
	elif i<100:
		fig.save('0'+str(i)+'.png')
	else:
		fig.save(str(i)+'.png')
	fig.close()


#os.system('ffmpeg -r 10 -i %03d.png -b:v 50M -r 10 '+args.n+'.mp4')

'''
if __name__ == "__main__":
	print select_name,select_tle
'''

