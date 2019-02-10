import h5py
import matplotlib.pyplot as plt 
import numpy as np 

f = h5py.File('./n2f_sphere-spectra-15-15000.0-50.h5', 'r')

ezr = np.array(f['ez.r'])
ezi = np.array(f['ez.i'])

s = len(ezr[0])

t = np.arange(1, s+1)
intensity = np.sqrt(ezr[0]**2 + ezi[0]**2)

plt.plot(t, intensity)

plt.savefig('hello.png')