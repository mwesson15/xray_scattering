import numpy as np
import meep as mp
import matplotlib.pyplot as plt

sx = 1500
sy = 1556.4
dpml = sy*0.02
resolution = 4.5/20

wavelength = 2.25*20 #E=5.5keV used in fig 4 in Soltwisch
f = 1/wavelength
df = 0.1
alpha = 0.86 #degrees
alpha = alpha*np.pi/180 #radians
k = mp.Vector3(y=np.sin(alpha)*2*np.pi*f)
cutoff = 1e20
amp_z = np.tan(alpha)

cell = mp.Vector3(sx, sy+2*dpml, 33/resolution)
pml_layers = [mp.PML(dpml, direction=mp.Y), mp.PML(15/resolution, direction=mp.Z)]

real_ind = 1 - 1.1676e-5
imag_ind = 3.9765e-7
real_epsilon = real_ind**2 - imag_ind**2
imag_epsilon = 2*real_ind*imag_ind
conductivity = 2*np.pi*f*imag_epsilon/real_epsilon

material = mp.Medium(epsilon=real_epsilon, D_conductivity=conductivity)

geometry = [mp.Block(mp.Vector3(291.6, 973.2, mp.inf),
                     center=mp.Vector3(-604.2, 0),
                     material=material),
            mp.Block(mp.Vector3(291.6, 973.2, mp.inf),
                     center=mp.Vector3(604.2, 0),
                     material=material),
            mp.Block(mp.Vector3(200, 91.6, mp.inf),
                     center=mp.Vector3(-650, -532.4),
                     material=material),
            mp.Block(mp.Vector3(200, 91.6, mp.inf),
                     center=mp.Vector3(650, -532.4),
                     material=material),
            mp.Cylinder(radius=91.6,
                        height=mp.inf,
                        axis=mp.Vector3(0, 0, 1),
                        center=mp.Vector3(-550, -486.6),
                        material=material),
            mp.Cylinder(radius=91.6,
                        height=mp.inf,
                        axis=mp.Vector3(0, 0, 1),
                        center=mp.Vector3(550, -486.6),
                        material=material),
            mp.Prism(vertices=[mp.Vector3(0, 0, 0), mp.Vector3(89.8, 0, 0), mp.Vector3(0, -973.2, 0)],
                     height=mp.inf,
                     axis=mp.Vector3(0, 0, 1),
                     center=mp.Vector3(-428.5, 162.2, 0),
                     material=material),
            mp.Prism(vertices=[mp.Vector3(0, 0, 0), mp.Vector3(89.8, 0, 0), mp.Vector3(89.8, -973.2, 0)],
                     height=mp.inf,
                     axis=mp.Vector3(0, 0, 1),
                     center=mp.Vector3(428.5, 162.2, 0),
                     material=material),
            mp.Block(mp.Vector3(1500, 291.6, mp.inf),
                     center=mp.Vector3(0, 486.6+291.6/2),
                     material=material),
            mp.Ellipsoid(center=mp.Vector3(0, 486.6, 0),
                         size=mp.Vector3(737.2, 260.4, mp.inf),
                         material=mp.Medium(epsilon=1))]

def pw_amp(k, x0):
    def _pw_amp(x):
        return np.exp(1j*k.dot(x+x0))
    return _pw_amp

sources = [
    mp.Source(
        mp.ContinuousSource(f),
        component=mp.Ex,
        center=mp.Vector3(),
        size=mp.Vector3(sx, sy),
        amp_func=pw_amp(k, mp.Vector3())),
    mp.Source(
        mp.ContinuousSource(f),
        component=mp.Ey,
        center=mp.Vector3(),
        size=mp.Vector3(sx, sy),
        amp_func=pw_amp(k, mp.Vector3())),
    mp.Source(
        mp.ContinuousSource(f),
        component=mp.Ez,
        center=mp.Vector3(),
        size=mp.Vector3(sx, sy),
        amp_func=pw_amp(k, mp.Vector3()),
        amplitude=np.tan(alpha))
]

k_point = mp.Vector3()
symmetry = [mp.Mirror(mp.Z)]

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    k_point=k_point,
                    symmetries=symmetry,
                    sources=sources,
                    resolution=resolution,
                    force_complex_fields=True,
                    Courant=0.49,
                    eps_averaging=True)

sim.run(until=2000)

Ey = sim.get_array(center=mp.Vector3(), size=mp.Vector3(sx, sy), component=mp.Ey, cmplx=True)
EyReal = Ey.real
EyImag = Ey.imag
Ey = np.absolute(Ey)
Ex = sim.get_array(center=mp.Vector3(), size=mp.Vector3(sx, sy), component=mp.Ex, cmplx=True)
ExReal = Ex.real
ExImag = Ex.imag
Ex = np.absolute(Ex)
Ez = sim.get_array(center=mp.Vector3(), size=mp.Vector3(sx, sy), component=mp.Ez, cmplx=True)
EzReal = Ez.real
EzImag = Ez.imag
Ez = np.absolute(Ez)
I = np.sqrt(np.square(Ey) + np.square(Ex) + np.square(Ez))
eps = sim.get_array(center=mp.Vector3(), size=mp.Vector3(sx, sy), component=mp.Dielectric, cmplx=False)

plt.figure(dpi=100)
#plt.imshow(eps.transpose(), cmap='Greys', alpha=0.9)
plt.imshow(I.transpose(), cmap='RdBu', alpha=0.8)
plt.colorbar()
plt.axis('off')
plt.savefig('Grating2Sim.png')
np.savetxt('Grating2SimIntensity.csv', I, delimiter=',')
print('Max Intensity = ' + str(I.max()))