import meep as mp
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

'''
TODO: make sure resolutions, lengthscales, source frequencies make sense
https://meep.readthedocs.io/en/latest/Introduction/#units-in-meep
'''

'''
TODO: make reflectance plots (Fig 5, Fig 6)
https://meep.readthedocs.io/en/latest/Python_Tutorials/Basics/#angular-reflectance-spectrum-of-a-planar-interface
https://meep.readthedocs.io/en/latest/Python_Tutorials/Mode_Decomposition/#reflectance-and-transmittance-spectra-for-planewave-at-oblique-incidence
'''

#soltwisch use finite element to do the simulation, not FDTD
#based on examples/binary_grating.py

#set up geometry
resolution = 10        # pixels/μm

#G2 params in 10^-8m
dpml = 1.             # PML thickness
dsub = 3.             # substrate thickness
dpad = 3.             # air padding between grating and PML at top
gp = 15.             # grating period
gh = 11.950               # grating height
lw = 6.730              #line width
tr = 0.916              #top radius
br = 1.302             #bottom radius
sa = 84.73              #sidewall angle (deg)

a = gh/np.tan(sa*2*np.pi/360) # intermediate for vertex calculation

sx = dpml+gp+dpml
sy = dpml+dsub+gh+dpad+dpml

cell = mp.Vector3(sx,sy,0)
pml_layers = [mp.PML(thickness=dpml)]
#pml_layers = [mp.PML(thickness=dpml,direction=mp.Y)]

#check we have right epsilon
substrate = mp.Medium(epsilon=12)
grating = mp.Medium(epsilon=12)

#one for substrate, one for grating superstructure, using trapezoid prism to approximate right now
vertices = [mp.Vector3(-0.5*lw-0.5*a,-0.5*gh,0),
            mp.Vector3(0.5*lw+0.5*a,-0.5*gh,0),
            mp.Vector3(0.5*lw-0.5*a,0.5*gh,0),
            mp.Vector3(-0.5*lw+0.5*a,0.5*gh,0)]

#mp.Block(material=grating, size=mp.Vector3(lw,gh,mp.inf), center=mp.Vector3(0,-0.5*sy+dpml+dsub+0.5*gh,0)),

geometry = [mp.Block(material=substrate, size=mp.Vector3(sx,dpml+dsub,mp.inf), center=mp.Vector3(0,-0.5*sy+0.5*(dpml+dsub),0)),
            mp.Prism(vertices, height=mp.inf, center=mp.Vector3(0,-0.5*sy+dpml+dsub+0.5*gh,0), material=grating)]
# symmetries=[mp.Mirror(mp.X)]
symmetries = []

#set up source
wvl_min = 0.4           # min wavelength
wvl_max = 0.6           # max wavelength
fmin = 1/wvl_max        # min frequency
fmax = 1/wvl_min        # max frequency
fcen = 0.5*(fmin+fmax)  # center frequency
df = fmax-fmin          # frequency width

theta = 0.86 * 2*np.pi/360 #in fig 5a these are in interval [0.3,1.2]
k = mp.Vector3(0,-np.sin(theta),np.cos(theta)) # complex fields for oblique angles
src_pt = mp.Vector3(0,-0.5*sy+dpml+dsub+gh,0) # incident at top of grating

#in meep, frequency units are in units of 2*pi*c
#read more carefully how they set up simulation of source in paper, this is just a plane hitting the top of grating
sources = [mp.Source(mp.ContinuousSource(frequency=100), component=mp.Ez, center=src_pt, size=mp.Vector3(sx,0,0))]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    k_point=k,
                    sources=sources,
                    symmetries=symmetries)
sim.run(until=20)

eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)

#plot geoms
plt.figure(dpi=100)
plt.imshow(eps_data.transpose()[::-1], interpolation='spline36', cmap='binary')
plt.axis('off')
plt.show()

#plot geoms with field
#is taking real part of field the right thing?
plt.figure(dpi=100)
plt.imshow(eps_data.transpose()[::-1], interpolation='spline36', cmap='binary')
plt.imshow(ez_data.real.transpose()[::-1], interpolation='spline36', cmap='RdBu', alpha=0.9)
plt.axis('off')
plt.show()

plt.figure(dpi=100)
plt.imshow(eps_data.transpose()[::-1], interpolation='spline36', cmap='binary')
plt.imshow(ez_data.imag.transpose()[::-1], interpolation='spline36', cmap='RdBu', alpha=0.9)
plt.axis('off')
plt.show()
