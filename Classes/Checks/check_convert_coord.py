import numpy as np
from general import general_functions as gf

#This code is meant to check the convert_coord() function found in general.py. It's used when converting from the (x,y)
#coordinates of the simulation and the (z, r) coordinates of the Gaussian Beam source

x = 0
y_min = -10
y_max = 10

origin = [-10, 3]

delta_y = 1
y = np.arange(y_min, y_max, delta_y)

coords = gf.convert_coord(x, y_min, y_max, origin, delta_y)

print(coords)

xy_mag = (x-origin[0])**2 + np.square(y-origin[1])
rz_mag = np.square(coords[0,:]) + np.square(coords[1,:])

print(xy_mag - rz_mag)
