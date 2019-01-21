import meep as mp
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

#soltwisch use finite element to do the simulation, not FDTD
#based on examples/binary_grating.py

#set up geometry
resolution = 10        # pixels/μm

dpml = 1.0             # PML thickness
dsub = 3.0             # substrate thickness
dpad = 3.0             # padding between grating and PML
gp = 10.0              # grating period
gh = 5.0               # grating height
gdc = 0.5              # grating duty cycle

sx = dpml+gp+dpml
sy = dpml+dsub+gh+dpad+dpml

cell = mp.Vector3(sx,sy,0)
pml_layers = [mp.PML(thickness=dpml,direction=mp.Y)]

#make sure materials are right
substrate = mp.Medium(epsilon=12)
grating = mp.Medium(epsilon=12)
#one for substrate, one for grating superstructure, using rectangles to approximate right now
geometry = [mp.Block(material=substrate, size=mp.Vector3(sx,dpml+dsub,mp.inf), center=mp.Vector3(0,-0.5*sy+0.5*(dpml+dsub),0)),
            mp.Block(material=grating, size=mp.Vector3(gdc*gp,gh,mp.inf), center=mp.Vector3(0,-0.5*sy+dpml+dsub+0.5*gh,0))]
symmetries=[mp.Mirror(mp.X)]

#set up source
wvl_min = 0.4           # min wavelength
wvl_max = 0.6           # max wavelength
fmin = 1/wvl_max        # min frequency
fmax = 1/wvl_min        # max frequency
fcen = 0.5*(fmin+fmax)  # center frequency
df = fmax-fmin          # frequency width

theta = 0.3 * 2*np.pi/360 #in fig 5a these are in interval [0.3,1.2]
k = mp.Vector3(0,-np.sin(theta),np.cos(theta)) # complex fields for oblique angles
src_pt = mp.Vector3(0,-0.5*sy+dpml+dsub+gh,0) # incident at top of grating

#in meep, frequency units are in units of 2*pi*c
#read more carefully how they set up simulation of source in paper, this is just a plane hitting the top of grating
sources = [mp.Source(mp.ContinuousSource(frequency=0.15), component=mp.Ez, center=src_pt, size=mp.Vector3(0.5*sx,0,0))]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    k_point=k,
                    sources=sources,
                    symmetries=symmetries)
sim.run(until=200)

eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
print(np.unique(eps_data))

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
