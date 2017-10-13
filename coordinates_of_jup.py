#This code aims to read in a list of time, calculate the coordinates of Jupiter at these moments and finally write down the results into a '.txt' file.

import date_sort
import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.coordinates import get_body
timelist=date_sort.date
record=open('/home/trina/Desktop/doc_loc.txt','a+')
for time in timelist:
	t=Time(time)
	loc = EarthLocation.of_site('greenwich') 
	now = get_body('jupiter', t, loc) 
	record.write(now.to_string('decimal'))
	'''record.write(now.ra.degree)
	record.write(now.dec.degree)'''
	record.write('\n')
record.close()





