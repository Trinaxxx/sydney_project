#This code aims to make a GIF with two FITS containing Jupiter, demonstrating this variable directly and lively. Particularly, I only focus on the local image where Jupiter are supposed to appear during that epoch in order to achieve a prominent visualization.

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import animation 
import numpy as np

#read_data
d1=fits.getdata('/home/trina/Desktop/the_movie/fits/freq1/2_1068235576_I.fits', ext=0)[1200:1601,800:1201]
d2=fits.getdata('/home/trina/Desktop/the_movie/fits/freq1/2_1078314912_I.fits', ext=0)[1200:1601,800:1201]
#text 
bbox_props = dict(boxstyle="circle,pad=0.3", fc="cyan", ec="b", lw=2)
#loop
ims = []
file=[d1,d2]
cord=[[63,165],[307,282]]
fig = plt.figure("Animation")
ax = fig.add_subplot(111)
for i in range(2):
     img=file[i]
     frame =  ax.imshow(img,origin='lower',cmap='gray')        
     t = ax.text(cord[i][0],cord[i][1], "Jupiter", ha="center", va="center", rotation=0,
            fontsize=2,bbox=bbox_props)
     ims.append([frame,t])
#for imgNum in range(numFiles):
    #fileName= files[imgNum]
    #img = read_image(fileName)
#animation   
anim = animation.ArtistAnimation(fig, ims, interval=350, blit=True, repeat_delay=350)
#anim.save('jupiter.gif',fps=1)
plt.show()
