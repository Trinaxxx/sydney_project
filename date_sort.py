#This code aims to sort a series of images according to the order of time.

import numpy as np
import pandas as pd

d1=np.genfromtxt('/home/trina/Desktop/doc_date.txt',delimiter='  ',dtype=None)
d2=pd.DataFrame(d1)
d3=d2.sort(columns=1,axis=0,ascending=True,inplace=False)
name=np.array(d3[0])
date=np.array(d3[1])

'''
date=open('/home/trina/Desktop/doc_date_sort.txt','a+')
for item in d4:
    date.write(item[1])
    date.write('\n')
date.close()
'''
