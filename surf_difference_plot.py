# DCW May 2015
# Plots an annual mean difference between two runs
# stippling where the change is NOT statistically significant
# Designed for N96 runs, UMvn8.4
import numpy as N
from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import colorConverter
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from surface_extract import extract
from shift_colormap import shiftedColorMap
# Setup
jobid = ['xlgac','xlgae']
userid = 'dcw32'
path = '/scratch/'+userid+'/netscratch/um/'
diagfile = '_cloud_diags'
diagname = 'cloud_cond_nuc'
colourbar_label = '% CCN increase'
# Number of years in run
nyrs = 5
# Sets the colour bar boundaries
levs=[-10,-5,0,5,15,25,35,45]
# If print_trop_diag == True, calculates and prints the volume weighted diag
print_surf_diag = True
# If percentage_diff == True, calculates the percentage difference
# Else calculates difference
percentage_diff = True
# End Setup

# Gumph to set up arrays etc
# Suitable for troposphere, volume file from N96 L63 low-top version
# vol is the volume of the gridboxes up to L63 in cm^3
volfile=Dataset('/scratch/dcw32/netscratch/um/n96_l63_geovol.nc')
vol=volfile.variables['vol_theta'][:]
vol=vol*1E6
# Set up empty arrays
diagh=N.empty([2,145,192])
stdev_diagh=N.empty([2,145,192])
stdevmask=N.ones([145,192])
# Sets the levels for the drawing of the stippling
levels=[0,0.5,1.01]

for i in range(2):
	# extracts the diagnostic of interest
	file=Dataset(path+jobid[i]+'/'+jobid[i]+diagfile+'.nc')
	diag=file.variables[diagname][:]
	file.close()
	# module turns diag into into climatology
	(diag,stdev_diag)=extract(diag,nyrs)
	# create mask of values above the surface layer
	mask=N.zeros([63,145,192])
	mask[0,:,:]=1
	# Limit to bottom 63 levels, bit of a fudge
	diag=diag[:63,:,:]
	stdev_diag=stdev_diag[:63,:,:]
	if print_surf_diag==True:
		diag_trop="%.2e" %(N.sum(diag * vol * mask)/N.sum(vol*mask))
		print str(diag_trop)+" surface diag"
	# Eliminate values above the surface
	diag=diag[0,:,:]
	stdev_diag=stdev_diag[0,:,:]
#	mass_box=vol*mask
#	diag_diff=oh_box*mass_box
#	ohstdmass=ohstd_box*mass_box
#	diag_diff=N.sum(diag_diff,axis=0)/N.sum(mass_box,axis=0)
#	ohstdmass=N.sum(ohstdmass,axis=0)/N.sum(mass_box,axis=0)
	diagh[i,:,:]=diag[:]
	stdev_diagh[i,:,:]=stdev_diag[:]
diag_diff=diagh[1,:,:]-diagh[0,:,:]
# use the 'base' for the significance calculation
stdev_diff=stdev_diagh[0,:,:]
# stipples if the change is not significant at the 95% level
for i in range(0,145):
	for j in range(0,192):
		if 2.5*stdev_diff[i,j]>N.abs(diag_diff[i,j]):
			stdevmask[i,j]=0
if percentage_diff == True:
	diag_diff=diag_diff/diagh[0,:,:]
	diag_diff=100*diag_diff
map=Basemap(lon_0=0,projection='robin')
map.drawcoastlines()
map.drawmapboundary()

# This does all the horrendous stuff to convert the array into something
# python can actually plot
file=Dataset('/scratch/dcw32/netscratch/um/xlgac/xlgac_dust_aod.nc')
lats=file.variables['latitude'][:]
lons=file.variables['longitude'][:]
lonsi=lons-360
lonsj=lons-360
diag_diff,lonsi=addcyclic(diag_diff,lonsi)
stdevmask,lonsj=addcyclic(stdevmask,lonsj)
diag_diff,lonsi=shiftgrid(-180.0,diag_diff,lonsi,cyclic=360.0)
stdevmask,lonsj=shiftgrid(-180.0,stdevmask,lonsj,cyclic=360.0)
lonsi,latsi=N.meshgrid(lonsi,lats)
a,b=map(lonsi,latsi)
# End horrendous stuff

cm=plt.cm.get_cmap('RdBu')

cm=shiftedColorMap(cm,midpoint=0.20)
cax=map.contourf(a,b,diag_diff,levs,cmap=cm,extend='both')

caz=map.contourf(a,b,stdevmask,hatches=['...',None],cmap=plt.get_cmap('gray'),alpha=0,extend='both')
cbar=plt.colorbar(cax,orientation="horizontal")
cbar.set_label(colourbar_label)
# You CANNOT print to .pdf with the stippling, eps and png fine (Python bug)
plt.savefig('/homes/'+userid+'/figures/cloud/ccn_diffplot_share.png',bbox_inches='tight')
plt.show()
