#This code aims to display the calculated positions of Jupiter in an epoch around(before and after) the observation times of two selected image, with obvious marks of the exact positions where Jupiter is supposed to be at the two moments.

import matplotlib.pyplot as plt
#from astropy import units as u
#from astropy.io import fits
import math
import numpy as np
import glob
import aplpy
import os


data=np.genfromtxt('doc_loc.txt',dtype=None,names=True,filling_values=-99)
ra=data['ra']
dec=data['dec']
files = glob.glob('*.fits') 
j=0
for fig in files:
	fig=aplpy.FITSFigure(fig)
	fig.set_xaxis_coord_type('longitude')
	fig.set_yaxis_coord_type('latitude')
	fig.tick_labels.set_xformat("hh:mm:ss.ss")
	fig.tick_labels.set_yformat("dd:mm:ss.ss")
	#grayscale
	fig.show_grayscale(pmin=5.0, pmax=99.5, invert=False, aspect='auto')
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
	if j==0:
		fig.show_circles(ra[6], dec[6], 0.1, edgecolor='r', linewidth=1.0, alpha=1.0)
	else:
		fig.show_circles(ra[14], dec[14], 0.1, edgecolor='r', linewidth=1.0, alpha=1.0)
	fig.save('00'+str(j)+'.png')
	fig.close()
	j+=1
os.system('ffmpeg -r 1 -i %03d.png -b:v 50M -r 10 jup_path_move.mp4')


'''
if j==0:
		for i in range(len(src_ra)):
			if i==6:
				fig.show_circles(src_ra[6], src_dec[6], 0.1, edgecolor='cyan', linewidth=1.0, alpha=1.0)
			else:
				fig.show_circles(src_ra[i], src_dec[i], 0.1, edgecolor='r', linewidth=1.0, alpha=1.0)
	else:
		for i in range(len(src_ra)):
			if i==14:
				fig.show_circles(src_ra[14], src_dec[14], 0.1, edgecolor='cyan', linewidth=1.0, alpha=1.0)
			else:
				fig.show_circles(src_ra[i], src_dec[i], 0.1, edgecolor='r', linewidth=1.0, alpha=1.0)	
	fig.save('00'+str(j)+'.png')
	fig.close()
	j+=1
'''


