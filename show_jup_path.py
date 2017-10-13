#This code aims to display the path where Jupiter moves in a stationary image.

import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits
import math
import numpy as np
import glob
import aplpy


fig=aplpy.FITSFigure('/home/trina/Desktop/the_movie/fits/freq1/2_1068235576_I.fits')
fig.set_xaxis_coord_type('longitude')
fig.set_yaxis_coord_type('latitude')
fig.tick_labels.set_xformat("hh:mm:ss.ss")
fig.tick_labels.set_yformat("dd:mm:ss.ss")
#grayscale
fig.show_grayscale(pmin=5.0, pmax=99.5, invert=False, aspect='auto')
#colorbar
fig.add_colorbar()
#fig.colorbar.set_axis_label_text('Flux (Jy/beam)')
src_ra =[100.89815095,103.80388791,106.90027557,108.99799125,110.69136961,111.71078113,111.94792362,111.31254852,109.94465037,107.8908186,105.7263204,103.58421575,102.07634707,101.33688768,101.17898453]#degrees
src_dec = [22.90967146,22.70110046,22.42367453,22.20612867,22.01800901,21.90860811,21.91100193,22.0323094,22.24774783,22.53380322,22.8003494,23.03412188,23.18469883,23.26024333,23.28716499]#degrees
#grid=True
#grid_setting
fig.add_grid()
fig.grid.show()
fig.grid.set_xspacing(5.0)  # degrees
fig.grid.set_yspacing(5.0)  # degrees
#circles
fig.show_circles(src_ra, src_dec, 0.1, edgecolor='r', linewidth=1.0, alpha=1.0)
fig.save('/home/trina/Desktop/final1.png')
fig.close()
	


