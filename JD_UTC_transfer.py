#This code aims to transfer the Julian Date to UTC time system.

import numpy
import astropy.io.fits as pf
from astropy.time import Time
import glob
date=open('/home/trina/Desktop/doc_date.txt','a+')
files = glob.glob('/home/trina/Desktop/new_jupiter/*.fits')
for f in files:
	hdr = pf.getheader(f,0) 
	jd = hdr["JD"] 
	t=Time(jd, format='jd', scale='utc') 
	date.write(f+' ')
	date.write(t.iso)
	date.write('\n')


'''
#example1
import astropy.io.fits as pf
hdr = pf.getheader("2_1068235576_I.fits”, 0) # Read the first FITS header for the file (first=0)
jd = hdr["JD”] # Get the Julian Date (“JD”) from the header
print(jd)
from astropy.time import Time
t=Time(jd, format='jd', scale='utc’) # Import the Julian Date into astropy
print(t.iso) # Convert the time into a UTC string

#example2
hdr = pf.getheader(f, 0) 
jd = hdr["JD"] 
jd_list=[]
for jd in jd_list:
	t=Time(jd, format='jd', scale='utc’) 
	doc1.write(t.iso)
	doc1.write('\n')
'''


	
