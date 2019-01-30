#This program checks the time averaged irradiance profile to ensure it is Gaussian with the proper parameters

import numpy as np
import matplotlib.pyplot as plt
from general import general_functions as gf

#Define beam parameters
w_0 = 5
l = 5
origin = [-10, 0]
omega = 2*np.pi
E_0 = 1
alpha = 0

#Define where we're looking
x = 0
y_min = -10
y_max = 10

#Define quantization of space/time
delta_y = 0.01
delta_t = 0.01

t = np.arange(0, 1+delta_t, delta_t)
y = np.arange(y_min, y_max+delta_y, delta_y)

#Create coordinate (r,z) matrix and calculate some stuf
coords = gf.convert_coord(x, y_min, y_max+delta_y, origin, alpha, delta_y)
zr = gf.z_R(w_0, l)
w = gf.w(coords[0,:], w_0, zr, l)
phi_t = gf.phi_t(coords[0,:], zr, l)
phi_l = gf.phi_l(coords[1,:], coords[0,:], zr, l)
phase_space = phi_t + phi_l
amp = gf.amp(E_0, w_0, w, coords[1,:])

#Define matrix that stores the intensity
I = np.zeros((y.shape[0], t.shape[0]))

#Go thru time and calculate the intensity at each point (go until omega*t=2*pi)
for i in range(t.shape[0]):
    E_t = gf.amp2E(amp, phase_space, omega, t[i])
    I_t = np.square(E_t)/377
    I[:, i] = I_t

#Calcuate the average intensity over the time period
I_avg = np.mean(I, axis=1)

#Calculate w and compare it to the expected value
#Find out where the intensity has decayed by 1/e^2
a = y[np.argmin(np.abs(I_avg - (I_avg.max()/np.exp(2))))]

#Calculate the FWHM
fwhm = 2*y[np.argmin(np.abs(I_avg-(I_avg.max()/2)))]

#Use the FWHM to estimate the beam waist (FWHM = w(z)/sqrt(2*ln2) )
check = fwhm/1.177

#Print the w(z) and the two methods for calculating the beam width (Should all be the same value ignoring the sign)
print(w.min())
print(a)
print(check)

#Plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(y, I_avg*1000, 'k.')
plt.axis('on')
ax.set_title('Time Averaged Intensity Profile\n of Gaussian Beam (z=10)')
ax.set_xlabel('r')
ax.set_ylabel('Time Averaged Intensity')

textfit1 = 'Expected: w = %.2f\n'  \
           'Model: w =%.2f\n'  \
        % (w.min(), a)
ax.text(0.05, .9, textfit1, transform=ax.transAxes, fontsize=11,
              verticalalignment='top')

plt.show()