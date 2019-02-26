import math
import cmath
import meep as mp
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

from midterm.materials import Silicon

#project z direction onto xy axis

#set up geometry
ts = 100               # timesteps
resolution = 10.
dpml = 1.             # PML thickness
dsub = 3.             # substrate thickness
dpad = 3.             # air padding between grating and PML at top

#G1 params in 10^-8m
gp = 10.             # grating period
gh = 10.271               # grating height
lw = 5.404              #line width
tr = 1.627              #top radius
br = 1.579             #bottom radius
sa = 90.91              #sidewall angle (deg)

# #G2 params in 10^-8m
# gp = 15.             # grating period
# gh = 11.950               # grating height
# lw = 6.730              #line width
# tr = 0.916              #top radius
# br = 1.302             #bottom radius
# sa = 84.73              #sidewall angle (deg)

tw = gh/np.tan(sa*2*np.pi/360) # part to subtract from top of grating width
a = (gh-tr)/np.tan(sa*2*np.pi/360) # intermediate for vertex calculation

#2 periods of grating
sx = 2*gp
sy = dpml+dsub+gh+dpad+dpml

cell_size = mp.Vector3(sx,sy,0)

pml_layers = [mp.PML(thickness=dpml, direction=mp.Y)]

#set up material properties
si = Silicon()
ind = si.index_refract()
ind = ind.real
cond = ind.imag
wv = 0.016869 #7.35 keV, relate to figure 6
freq = 1/wv

substrate = mp.Medium(index=ind, D_conductivity=2*math.pi*freq*cond/ind)
grating = mp.Medium(index=ind, D_conductivity=2*math.pi*freq*cond/ind)

#one for substrate, one for grating superstructure, using trapezoid prism to approximate right now
vtx = [mp.Vector3(-0.5*lw-0.5*a+gp/2,-0.5*gh,0),
            mp.Vector3(0.5*lw+0.5*a+gp/2,-0.5*gh,0),
            mp.Vector3(0.5*lw-0.5*a+gp/2,0.5*gh-tr,0),
            mp.Vector3(-0.5*lw+0.5*a+gp/2,0.5*gh-tr,0)]

vtx2 = [mp.Vector3(-0.5*lw-0.5*a-gp/2,-0.5*gh,0),
            mp.Vector3(0.5*lw+0.5*a-gp/2,-0.5*gh,0),
            mp.Vector3(0.5*lw-0.5*a-gp/2,0.5*gh-tr,0),
            mp.Vector3(-0.5*lw+0.5*a-gp/2,0.5*gh-tr,0)]

#rounded corners
c1 = mp.Cylinder(radius=tr,
            height=mp.inf,
            axis=mp.Vector3(0,0,1),
            center=mp.Vector3(-gp/2-lw/2+tw/2+tr,-sy/2+dpml+dsub+gh-tr, 0),
            material=grating)
c2 = mp.Cylinder(radius=tr,
            height=mp.inf,
            axis=mp.Vector3(0,0,1),
            center=mp.Vector3(-gp/2+lw/2-tw/2-tr,-sy/2+dpml+dsub+gh-tr, 0),
            material=grating)
c3 = mp.Cylinder(radius=tr,
            height=mp.inf,
            axis=mp.Vector3(0,0,1),
            center=mp.Vector3(gp/2-lw/2+tw/2+tr,-sy/2+dpml+dsub+gh-tr, 0),
            material=grating)
c4 = mp.Cylinder(radius=tr,
            height=mp.inf,
            axis=mp.Vector3(0,0,1),
            center=mp.Vector3(gp/2+lw/2-tw/2-tr,-sy/2+dpml+dsub+gh-tr, 0),
            material=grating)

#top block
b1 = mp.Block(center=mp.Vector3(-gp/2,-sy/2+dpml+dsub+gh-tr/2), size=mp.Vector3(lw-tw-2*tr,tr,mp.inf), material=grating)
b2 = mp.Block(center=mp.Vector3(gp/2,-sy/2+dpml+dsub+gh-tr/2), size=mp.Vector3(lw-tw-2*tr,tr,mp.inf), material=grating)

#ellipsoid cutout
e1 = mp.Ellipsoid(center=mp.Vector3(0,-sy/2+dpml+dsub,0), size=mp.Vector3(gp-lw-a,br,mp.inf), material=mp.Medium(epsilon=1))
e2 = mp.Ellipsoid(center=mp.Vector3(-gp,-sy/2+dpml+dsub,0), size=mp.Vector3(gp-lw-a,br,mp.inf), material=mp.Medium(epsilon=1))
e3 = mp.Ellipsoid(center=mp.Vector3(gp,-sy/2+dpml+dsub,0), size=mp.Vector3(gp-lw-a,br,mp.inf), material=mp.Medium(epsilon=1))


geometry = [mp.Block(material=substrate, size=mp.Vector3(sx,dpml+dsub,mp.inf), center=mp.Vector3(0,-0.5*sy+0.5*(dpml+dsub),0)),
            mp.Prism(vtx, height=mp.inf, center=mp.Vector3(0,-0.5*sy+dpml+dsub+0.5*gh,0), material=grating),
            mp.Prism(vtx2, height=mp.inf, center=mp.Vector3(0,-0.5*sy+dpml+dsub+0.5*gh,0), material=grating),
            c1, c2, c3, c4, b1, b2, e1, e2, e3]
sym = []

#set up source
ng = 1 #index of medium wave is propagating in

theta = 1.09 * 2*np.pi/360 #in fig 6 these are in interval [0.3,1.2]
k = mp.Vector3(0,-np.sin(theta),0).scale(freq*ng) # scaled field incidence
k = mp.Vector3(0,0,0)
src_pt = mp.Vector3(0,sy/2-dpml,0) # plane of incidence top of cell

def pw_amp(k,x0):
    def _pw_amp(x):
        return cmath.exp(1j*2*math.pi*k.dot(x+x0))
    return _pw_amp

#fix this so that it's all the vacuum region
'''Energy isn't allowed to dissipate??? the value of E
should be one but somehow max is 10^10?? (ask jason)'''
sources = [mp.Source(mp.ContinuousSource(freq),
                     component=-mp.Ey,
                     amplitude=np.sin(theta),
                     center=src_pt,
                     size=mp.Vector3(sx,0,0),
                     amp_func=pw_amp(k,src_pt))]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=[],
                    sources=sources,
                    k_point=k,
                    symmetries=sym)

nf_pt = mp.Vector3(0,sy/2-dpml-dpad/2,0)
_nearfield = sim.add_near2far(freq, 0, 1, mp.Near2FarRegion(nf_pt, size=mp.Vector3(sx,0,0), weight=1))

sim.run(until=ts)

_E = np.zeros((1001,1),dtype=np.complex128)
'''todo: figure out what the best distance d is (do the transform described for ewald sphere)'''
d = 1e5*np.sin(theta)
for n in range(-500,501):
    ffy = sim.get_farfield(_nearfield, mp.Vector3(n/10,d,0))
    _E[n] = ffy[1]

_eps_data = sim.get_array(center=mp.Vector3(), size=cell_size, component=mp.Dielectric)
_ez_data = sim.get_array(center=mp.Vector3(), size=cell_size, component=mp.Ey)

sim.reset_meep()

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    k_point=k,
                    symmetries=sym)

nearfield = sim.add_near2far(freq, 0, 1, mp.Near2FarRegion(nf_pt, size=mp.Vector3(sx,0,0), weight=1))

sim.run(until=ts)

E = np.zeros((1001,1),dtype=np.complex128)
for n in range(-500,501):
    ffy = sim.get_farfield(nearfield, mp.Vector3(n/10,d,0))
    E[n] = ffy[1]


'''looks like symmetry is broken because something is ever so slightly off-center (discretization issue)'''
def Iq(x,E):
    c = 1
    n = 1
    e0 = 1
    A = c*n*e0/2
    I = A*(np.abs(E)**2)
    I = I/np.amax(I)
    q = (2*np.pi/wv)*np.arctan(x/d)

    plt.plot(q/10,I)
    plt.xlim(-0.7,0.7)
    plt.ylim(0,1)
    plt.ylabel('Rel. Intensity')
    plt.xlabel('q (1/nm)')
    plt.show()

eps_data = sim.get_array(center=mp.Vector3(), size=cell_size, component=mp.Dielectric)
ez_data = sim.get_array(center=mp.Vector3(), size=cell_size, component=mp.Ey)

#plot geoms with field
plt.figure(dpi=100)
plt.imshow(eps_data.transpose()[::-1], interpolation='spline36', cmap='binary')
plt.imshow((np.abs(ez_data)**2).transpose()[::-1], interpolation='spline36', cmap='RdBu', alpha=0.9)
plt.axis('off')
plt.show()

to_plot = E - _E

#background
Iq(np.linspace(-50,50,1001), _E)
#scattered
Iq(np.linspace(-50,50,1001),to_plot)
