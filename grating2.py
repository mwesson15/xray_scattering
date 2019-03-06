import numpy as np
import meep as mp
import matplotlib.pyplot as plt

sx = 1500
sy = 1556.4
dpml = 50
resolution = 4.5/20
S = 0.49

alpha = 0.86 #degrees
alpha = alpha*np.pi/180 #radians

e = 0
Energies = [5.5, 5.55, 5.6]
wavelengths = [2.25, 2.23, 2.21]
delta = [1.1676e-5, 1.6068e-5, 1.5779e-5]
beta = [3.9765e-7, 7.3314e-7, 7.082e-7]

E = Energies[e]
wavelength = wavelengths[e]
f = 1/wavelength
width = 1000
df = 1/width
alpha = 0.86 #degrees
alpha = alpha*np.pi/180 #radians
k = mp.Vector3(y=np.sin(alpha)*2*np.pi*f)

cell = mp.Vector3(sx, sy+2*dpml)
pml_layers = [mp.PML(dpml, direction=mp.Y)]

real_ind = 1 - delta[e]
imag_ind = beta[e]
real_epsilon = real_ind**2 - imag_ind**2
imag_epsilon = 2*real_ind*imag_ind
omega = 2*np.pi*f*1e10
conductivity = omega*imag_epsilon
factor = (1e-10/3e8)/(real_epsilon*8.85e-12)
D_conductivity = factor*conductivity

material = mp.Medium(epsilon=real_epsilon, D_conductivity=D_conductivity)

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

f = f*np.sin(alpha)

sources = [
    mp.Source(
        mp.GaussianSource(f, df),
        component=mp.Hz,
        center=mp.Vector3(y=-sy/2),
        size=mp.Vector3(sx),
        amplitude=1),
    mp.Source(
        mp.GaussianSource(f, df),
        component=mp.Ez,
        center=mp.Vector3(y=-sy/2),
        size=mp.Vector3(sx),
        amplitude=1)
]

k_point = mp.Vector3(0, 0, 0)

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=[],
                    k_point=k_point,
                    sources=sources,
                    resolution=resolution,
                    force_complex_fields=True,
                    Courant=S,
                    eps_averaging=True)

nearfield = sim.add_near2far(f, 0, 1,
        mp.Near2FarRegion(mp.Vector3(y=-sy/2+100), size=mp.Vector3(sx, 200), weight=1, direction=mp.Z),
        mp.Near2FarRegion(mp.Vector3(y=-46), size=mp.Vector3(737, 1064.4), weight=1, direction=mp.Z))

sim.run(until_after_sources=mp.stop_when_fields_decayed(1000, mp.Hz, mp.Vector3(), 1e-3))
#sim.run(until=5000)

Ey = sim.get_array(center=mp.Vector3(), size=mp.Vector3(sx, sy), component=mp.Ey, cmplx=True)
EyReal = Ey.real
Ey = np.absolute(Ey)
Ex = sim.get_array(center=mp.Vector3(), size=mp.Vector3(sx, sy), component=mp.Ex, cmplx=True)
ExReal = Ex.real
Ex = np.absolute(Ex)
EReal = np.square(ExReal) + np.square(EyReal)
Ez = sim.get_array(center=mp.Vector3(), size=mp.Vector3(sx, sy), component=mp.Ez, cmplx=True)
EzReal = Ez.real
Ez = np.absolute(Ez)
EReal = np.square(ExReal) + np.square(EyReal) + np.square(EzReal)
I = np.sqrt(np.square(Ey) + np.square(Ex) + np.square(Ez))
eps = sim.get_array(center=mp.Vector3(), size=mp.Vector3(sx, sy), component=mp.Dielectric)

Distance2Detector = 1e8

#Generate the x values along the range and then solve for the corresponding y value
Radius = Distance2Detector*np.tan(alpha)
ResolutionX = Radius/500
XValues = np.arange(-Radius, Radius, ResolutionX)
YValues = np.sqrt(Radius**2 - np.square(XValues))
YValues = -YValues
steps = XValues.shape[0]

FarEx = np.zeros(XValues.shape, dtype='Complex128')
FarEy = FarEx
FarEz = FarEx

for i, x in enumerate(XValues):
    if i % 5 == 0:
        print('Near2Far Field on Step ' + str(i) + ' of ' + str(steps))
    y = YValues[i]
    ff = sim.get_farfield(nearfield, mp.Vector3(x, y, Distance2Detector))
    FarEx[i] = ff[0]
    FarEy[i] = ff[1]
    FarEz[i] = ff[2]

FarEx = np.absolute(FarEx)
FarEy = np.absolute(FarEy)
FarEz = np.absolute(FarEz)

FarI = FarEx + FarEy + FarEz
FarI = FarI.astype(np.single)
#FarI = 1.8*FarI/FarI.max()

q = (2*np.pi/wavelength)*(XValues/Radius)*np.cos(np.arcsin(YValues/Radius))/10

#Delete the data points in the center
cutoff = 0.4
#Low = np.argmin(np.abs(q+cutoff))
#Hi = np.argmin(np.abs(q-cutoff))


#plt.plot(q[:Low], FarI[:Low])
plt.plot(q, FarI)
plt.title('Far Field of E = ' + str(E) + 'keV')
plt.xlabel('$\\Delta q_x (nm^{-1})$')
plt.ylabel('Relative Intensity')
plt.axis('on')
plt.show()