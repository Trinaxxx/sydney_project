#This code aims to make a movie demonstrating the retrograde motion of Jupiter with a series of images.

import math
import glob
import aplpy
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import date_sort
import os

data=np.genfromtxt('/home/trina/Desktop/doc_loc.txt',dtype=None,names=True,filling_values=-99)
ra=data['ra']
dec=data['dec']
files=date_sort.name
j=0
for fig in files:
	fig=aplpy.FITSFigure(fig)
	fig.set_xaxis_coord_type('longitude')
	fig.set_yaxis_coord_type('latitude')
	fig.tick_labels.set_xformat("hh:mm:ss.ss")
	fig.tick_labels.set_yformat("dd:mm:ss.ss")
	#grayscale
	fig.show_grayscale(pmin=5.0, pmax=99.5, invert=True, aspect='auto')
	#colorbar
	fig.add_colorbar()
	#fig.colorbar.set_axis_label_text('Flux (Jy/beam)')
	#grid=True
	#grid_setting
	fig.add_grid()
	fig.grid.show()
	fig.grid.set_xspacing(5.0)  # degrees
	fig.grid.set_yspacing(5.0)  # degrees
	#circles
	fig.show_circles(ra, dec, 0.1, edgecolor='cyan', linewidth=1.0, alpha=1.0)
	fig.show_circles(ra[j], dec[j], 0.1, edgecolor='r', linewidth=1.0, alpha=1.0)
	if j<10:
		fig.save(str(j)+'00'+'.png')
	else:
		fig.save(str(j)+'0'+'.png')
	fig.close()
	j+=1
os.system('ffmpeg -r 10 -i %03d.png -b:v 50M -r 10 xxxx.mp4')

