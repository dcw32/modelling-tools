# Turns a 3d field into a climatology
from netCDF4 import Dataset
import numpy as N
def extract(oh,nyrs):
	ohj=N.empty([nyrs,85,145,192])
	for i in range (0,nyrs):
		q=12*i
		ohi=oh[q:q+12,:,:,:]
		ohj[i,:,:,:]=N.mean(ohi,axis=0)
	ohstd=N.std(ohj,axis=0)
	oh=N.mean(ohj,axis=0)
	return (oh,ohstd)
