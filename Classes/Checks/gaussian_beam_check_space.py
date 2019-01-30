import numpy as np
import matplotlib.pyplot as plt
from general import general_functions as gf

#This code is meant to display a Gaussian beam over a 2D space and check its characteristics
#Define some variables
delta_y = 1

#Define wave properties
origin = [0, 0]
alpha = 0
l = 30
w_0 = 40
E_0 = 10
omega = 2*np.pi

#Use t = 0
t = 0

#Define bounds for column
x_min = -400
x_max = 400
y_min = -200
y_max = 200

#Create the space
x = np.arange(x_min, x_max, delta_y)
y = np.arange(y_min, y_max, delta_y)

#Declare matrix that will store the I values
I_space = np.zeros((y.shape[0], x.shape[0]))

#Declare empty list that stores w values
wz = []

#Go through each x value and determine the Intensity along the row
for i in range(x.shape[0]):
    # Get (r,z) coordinates
    coords = gf.convert_coord(x[i], y_min, y_max, origin, alpha, delta_y)

    # Calculate and store some beam parameters
    zr = gf.z_R(w_0, l)
    w = gf.w(coords[0, :], w_0, zr, l)
    wz.append(w)
    phi_t = gf.phi_t(coords[0, :], zr, l)
    phi_l = gf.phi_l(coords[1, :], coords[0, :], zr, l)

    # Calculate the spatial factors
    amp = gf.amp(E_0, w_0, w, coords[1, :])

    # Convert these to E
    E = gf.amp2E(amp, phi_t+phi_l, omega, t)

    # Calculate H and multiply E and H to get the intensity (Impedance of free space=377)
    H = E / 377
    I = np.multiply(E, H)

    I_space[:,i] = I

#Plot I
plt.figure(dpi=100)
plt.imshow(I_space, interpolation='spline36', cmap='RdBu', alpha=0.9)
plt.axis('off')
plt.title('Intensity of Gaussian Beam $(\\alpha = 0)$')
plt.colorbar()
plt.show()

#Calculate the divergence of the beam
divergence = l/(np.pi*w_0)

#Calculate the divergence between each x value and the next
wz = np.array(wz)
mydiv = np.arctan((wz - np.roll(wz, 1))/(delta_y))
mydiv = np.delete(mydiv, 0)

#Compare the two (they should converge as x_max increases
print(divergence)
print(mydiv.max())

#Ran program with multiple x_max's and the two do converge