# Shift the colormap to a defined centre, for example when doing
# a change plot where there's a bias (e.g. changes between -0.1 to 0.3)

# Adapted from somewhere online...
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
import math as m
def latlon_shift(lons,lats,data):
        lonsi=lons-360
        data,lonsi=addcyclic(data,lonsi)
        data,lonsi=shiftgrid(-180.,data,lonsi,cyclic=360.)
        lonsi,latsi=np.meshgrid(lonsi,lats)
        return(lonsi,latsi,data)
def gbox_areas(x,y):
# lats x lons
        area=np.zeros([x,y])
        R=6.371E6
        for j in range(x):
                area[j,:]=(R**2)*m.radians(360./y)*(m.sin(m.radians(90.-(j-0.5)*180./(x-1)))-m.sin(m.radians(90.-(180./(x-1))*(j+0.5))))
        return area
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }
    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)
    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])
    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)
        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))
    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)
    return newcmap
