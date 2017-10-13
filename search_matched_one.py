#This code aims to conduct a crossmatch for the brightest detection of the given image with bulk catalog of satellites and then return the information about the identified source/satellite.

import os
import numpy as np
import astropy.io.fits as pf
import getopt, datetime, pytz
import ephem
import aplpy
import argparse
import astropy.wcs as wcs
from astropy.coordinates import SkyCoord
from astropy import units as u

#use in command line
parser = argparse.ArgumentParser(description='search interlopers.')
parser.add_argument("id",type=int,help="the number representing the time that the observation started")
args=parser.parse_args()

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
obs_id = args.id
snapshot = "/home/trina/Desktop/2_%s_V.fits" %(obs_id)
obs_time = gps_datetime(obs_id)
end_time = obs_time + datetime.timedelta(seconds=112)
delta_time = datetime.timedelta(seconds=0.1)
a=datetime.datetime(obs_time.year,obs_time.month,obs_time.day-2).date()
b=datetime.datetime(obs_time.year,obs_time.month,obs_time.day+2).date()

#TLE-DATA-Download#
#DOWNLOADS
os.system("curl -c cookies.txt -b cookies.txt https://www.space-track.org/ajaxauth/login -d 'identity=zsnxxx@sina.com&password=Darkness0924zsn'")
os.system('curl --limit-rate 100K --cookie cookies.txt https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/'+str(a)+'--'+str(b)+'/orderby/NORAD_CAT_ID/format/3le > /home/trina/Desktop/catalog.txt')#specific epoch download
#os.system('curl --limit-rate 100K --cookie cookies.txt https://www.space-track.org/basicspacedata/query/class/tle_latest/ORDINAL/1/EPOCH/%3Enow-30/orderby/NORAD_CAT_ID/format/3le > /home/trina/Desktop/catalog.txt')-------full catalog

#https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/2013-10-26--2013-12-21/orderby/NORAD_CAT_ID/format/3le -----repeated
#https://www.space-track.org/basicspacedata/query/class/tle_latest/ORDINAL/1/EPOCH/2013-10-26--2013-12-21/orderby/NORAD_CAT_ID/format/3le -----wrong list without STRELA3


#read_in
f=open('/home/trina/Desktop/catalog.txt','r') #only read
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
	timelist=[]
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
		# Get the RA and Dec in radians, convert to degrees and add to the ra/dec list; get the timelist
		timelist.append(t)
		l_ra.append(float(sat_ephem.a_ra) / ephem.degree)
		l_dec.append(float(sat_ephem.a_dec) / ephem.degree)
		t += delta_time
	return l_ra, l_dec, timelist

#get timelist
atle=finaltlelist[0]
ep = ephem.readtle(atle[0],atle[1],atle[2])
timelist = get_sat_radec(ep, obs_time, end_time, delta_time)[2]
#get cordlist
#loop
cordlist,namelist=[],[]
for tle in finaltlelist:
	sat_ephem = ephem.readtle(tle[0],tle[1],tle[2])
	#l_ra, l_dec = get_sat_radec(sat_ephem, obs_time, end_time, delta_time)
	l_ra,l_dec,tt= get_sat_radec(sat_ephem, obs_time, end_time, delta_time)
	ob_pos=(l_ra,l_dec)
	cordlist.append(ob_pos)
	namelist.append(tle[0])
###the code above can be revised by list comprehension###

#center_cord
f=snapshot
hdr = pf.getheader(f,0) 
ra_c=hdr["RA_PH"]
dec_c=hdr["DEC_PH"] 

#smaller than search_radius
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

def closest_to_central(tup,x,y):
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
  	select_namelist,select_tlelist,select_pos_list,linelist=[],[],[],[]
	i=0
    	for ob_pos in cordlist:
  		min_d,minra,mindec=closest_to_central(ob_pos,ra_c,dec_c)
  	  	if min_d<max_d:
  	    		#t=(i,namelist[i],min_d,minra,mindec)
			t=(i,namelist[i])
			select_namelist.append(t)
			select_tlelist.append(finaltlelist[i])
			select_pos_list.append(ob_pos)
			line=np.array(ob_pos)
	                linelist.append(line)
		i+=1
  	return select_namelist,select_tlelist,select_pos_list,linelist

#get select all lists
select_name,select_tle,select_pos,linelist=judgment(cordlist,10)
select_ephem=[]
for tle in select_tle:
	e = ephem.readtle(tle[0],tle[1],tle[2])
	select_ephem.append(e)
	
#find closest satellite to detection source
#get pixel_cord of detection source---for normal ones
def load_fits(file):
  data=pf.open(file)[0].data
  min=data[0][0]
  for r in range(len(data)):
    for c in range(len(data[0])):
      if data[r][c]<=min:
        min=data[r][c]
        i=r
        j=c
  location=(j,i)#FITS pixel cords need to be reversed
  return location


#get pixel_cord of detection source---for those who have side noise
#method1:set the side noise as zero
def find_peak(f, corner, window):
	v_file = pf.open(f)
	v_hdr = v_file[0].header
	nx = int(v_hdr["NAXIS1"])
	ny = int(v_hdr["NAXIS2"])
	v_wcs = wcs.WCS(v_hdr)
	v_data = np.abs(v_file[0].data) ################absolute value
	# clear out corner
	v_data[0:corner,0:corner] = 0.0
	v_data[ny-corner:ny,0:corner] = 0.0
	v_data[0:corner,nx-corner:nx] = 0.0
	v_data[ny-corner:ny,nx-corner:nx] = 0.0
	# clear out edge window
	v_data[0:ny,0:window] = 0.0
	v_data[0:ny,nx-window:nx] = 0.0
	v_data[0:window,0:nx] = 0.0
	v_data[ny-window:ny,0:nx] = 0.0
	# see where the peak is now
	v_peak = np.argmax(v_data)
	v_max = np.max(v_data)
	maxy, maxx = np.unravel_index(v_peak, v_data.shape)
	#
	sky = v_wcs.wcs_pix2world(maxx, maxy, 0)
	sc = SkyCoord(sky[0], sky[1], frame='fk5', unit=(u.deg, u.deg))
	return maxx, maxy, sc

#method2:by hand
x=[1758,889,1455]
y=[1317,466,1270]

#pixel to world transfer
fig=aplpy.FITSFigure(snapshot)

#normal
#npx,npy = load_fits(snapshot)
#x_world, y_world = fig.pixel2world(npx,npy)

#side noise
x,y,sc=find_peak(snapshot, 250, 70)
x_world, y_world = fig.pixel2world(x,y)

#find pass time index and distance for one sat
def find_pass(tup,x,y):
  	dislist=[]
  	for i in range(len(tup[0])):
    		temp=angular_dist(tup[0][i],tup[1][i],x,y)
    		dislist.append(temp)
	clo_dis=min(dislist)
	clo_time_index=dislist.index(clo_dis)
	return clo_time_index,clo_dis
  	
#find the pass satellite
def find_pass_sat(select_pos):
	mindlist,ptilist=[],[]
    	for ob_pos in select_pos:
		temptimeindex,tempdis=find_pass(ob_pos,x_world,y_world)
		ptilist.append(temptimeindex)
		mindlist.append(tempdis)
	min_dis=min(mindlist)
	min_index=mindlist.index(min_dis)
	return min_index,min_dis,ptilist

sat_index,delta_d,ptilist=find_pass_sat(select_pos)
pass_time_index=ptilist[sat_index]
pass_time=timelist[pass_time_index]

#get right_time cords list	
def get_stationary_radec(selec_pos,index):
	sta_ra, sta_dec=[],[]
	for ob_pos in select_pos:
		r=ob_pos[0][index]
		d=ob_pos[1][index]
		sta_ra.append(r)
		sta_dec.append(d)
	return sta_ra,sta_dec

sta_ra,sta_dec=get_stationary_radec(select_pos,pass_time_index)

###########plot###########

##plot stationary one
fig.tick_labels.set_xformat("hh:mm")
fig.tick_labels.set_yformat("dd:mm")
#add title
fig.add_label(0.5, 1.05, "time: %s" %(pass_time), relative=True, size='xx-large', layer='title')
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
for i in range(len(sta_ra)):
	if i==sat_index:
		fig.show_circles(sta_ra[i], sta_dec[i], 0.2, color='red',linewidth=1.0, alpha=0.5)
	else:
		fig.show_circles(sta_ra[i], sta_dec[i], 0.2, color='blue',linewidth=1.0, alpha=1.0)
#fig.show_circles(sta_ra, sta_dec, 0.1, color='r',linewidth=1.0, alpha=1.0)
for line in linelist:
	fig.show_lines([line], color='c', linewidth=2.0, alpha=0.5)
fig.save('/home/trina/Desktop/2_%s_V.png' %(obs_id))
fig.close()

'''
##plot stationary path
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


#make movive
for i in range(len(select_pos[0][0])):
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
		fig.show_circles(pos[0][i], pos[1][i],0.1, color='r',linewidth=1.0, alpha=0.5)
	if i<10:
		fig.save('00'+str(i)+'.png')
	elif i<100:
		fig.save('0'+str(i)+'.png')
	else:
		fig.save(str(i)+'.png')
	fig.close()


#os.system('ffmpeg -r 10 -i %03d.png -b:v 50M -r 10 '+args.n+'.mp4')


#make animation
ilist=glob.glob('*.png')
imagelist=np.sort(ilist)
images = []
for filename in imagelist:
    images.append(imageio.imread(filename))
exportname = "final.gif"
imageio.mimsave(exportname, images, 'GIF')
#kargs = { 'duration': 0.5 }
#imageio.mimsave(exportname, images, 'GIF', **kargs)

'''
if __name__ == "__main__":
	print obs_id
	print x_world, y_world
        print select_name[sat_index]
	print sta_ra[sat_index], sta_dec[sat_index]
	print delta_d
	print pass_time
