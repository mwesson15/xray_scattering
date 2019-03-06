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

e = 2
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
        mp.ContinuousSource(f),
        component=mp.Hz,
        center=mp.Vector3(y=-sy/2),
        size=mp.Vector3(sx),
        amplitude=1),
    mp.Source(
        mp.ContinuousSource(f),
        component=mp.Ez,
        center=mp.Vector3(y=-sy/2),
        size=mp.Vector3(sx),
        amplitude=2*np.tan(alpha))
]

k_point = mp.Vector3()

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    k_point=k_point,
                    sources=sources,
                    resolution=resolution,
                    force_complex_fields=True,
                    Courant=S,
                    eps_averaging=True)

nearfield = sim.add_near2far(f, 0, 1,
        mp.Near2FarRegion(mp.Vector3(y=-sy/2+100), size=mp.Vector3(sx, 200), weight=1, direction=mp.Z),
        mp.Near2FarRegion(mp.Vector3(y=-46), size=mp.Vector3(737, 1064.4), weight=1, direction=mp.Z))

#sim.run(until_after_sources=mp.stop_when_fields_decayed(1000, mp.Hz, mp.Vector3(), 1e-3))
sim.run(until=5000)

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

Distance2Detector = 1e8 #one micron, when did you change it
DetectorHeight = 6e6
DetectorResolutionY = DetectorHeight/50
DetectorWidth = DetectorHeight
DetectorResolutionX = DetectorWidth/50
DetectorDepth = wavelength*1*0
DetectorResolutionZ = wavelength/10

DetectorYValues = np.arange(-(DetectorHeight/2 + DetectorResolutionY), 0, DetectorResolutionY)

DetectorXValues = np.arange(-DetectorWidth/2, DetectorWidth/2 + DetectorResolutionX, DetectorResolutionX)
DetectorZValues = np.arange(-DetectorDepth/2, DetectorDepth/2 + DetectorResolutionZ, DetectorResolutionZ)
steps = DetectorXValues.shape[0]

FarEx = np.zeros((DetectorYValues.shape[0], DetectorXValues.shape[0], DetectorZValues.shape[0]), dtype='Complex128')
FarEy = FarEx
FarEz = FarEx

for k, z in enumerate(DetectorZValues):
    for j, x in enumerate(DetectorXValues):
        if j % 5 == 0:
            print('Near2Far Field on Step ' + str(j) + ' of ' + str(steps))
        for i, y in enumerate(DetectorYValues):
            ff = sim.get_farfield(nearfield, mp.Vector3(x, y, Distance2Detector + z))
            FarEx[i, j, k] = ff[0]
            FarEy[i, j, k] = ff[1]
            FarEz[i, j, k] = ff[2]

FarEx = np.absolute(FarEx)
FarEy = np.absolute(FarEy)
FarEz = np.absolute(FarEz)

FarE = FarEx + FarEy + FarEz
FarE = FarE.real
Detector = np.mean(FarE, axis=2)
Detector = Detector.astype(np.single)

plt.figure(dpi=100)
plt.imshow(eps.transpose(), cmap='binary')
plt.imshow(I.transpose(), cmap='gist_heat', alpha=0.9)
plt.colorbar()
plt.axis('on')
plt.title('Time-Averaged Electric Field\n (E = ' + str(E) + 'keV)')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
plt.savefig('Grating2Sim.png')
#np.savetxt('GratingFarE.csv', Detector, delimiter=',')

def QR(x, y, Detector, d):
    I = (np.square(Detector)) / 2
    I = 1.8 * I / np.amax(I)
    coord = np.ones((y.shape[0], x.shape[0], 2))
    for i in range(x.shape[0]):
        coord[:, i, 0] = x[i]
    for i in range(y.shape[0]):
        coord[i, :, 1] = y[i]
    Radius = (d*np.tan(alpha))**2
    r = np.square(coord[:, :, 0]) + np.square(coord[:, :, 1])
    theta = np.arctan(coord[:, :, 0] / d)
    alpha_f = np.arctan(coord[:, :, 1] / d)
    r = np.absolute(r-Radius)
    qx = np.multiply(np.sin(theta), np.cos(alpha_f))
    I_0 = np.ones(r.shape[1])
    q = np.ones(r.shape[1])
    for i in range(r.shape[1]):
        index = np.argmin(r[:, i])
        I_0[i] = I[index, i]
        q[i] = qx[index, i]
    midpt = int(coord.shape[1]/2)
    LowIndex = np.argmin(np.abs(coord[0, :midpt, 0] - Radius**0.5))
    HiIndex = np.argmin(np.abs(coord[0,midpt:, 0] - Radius**0.5))
    q = q[HiIndex:HiIndex+midpt]
    I_0 = I_0[HiIndex:HiIndex+midpt]
    q = (2*np.pi*q/wavelength)*10
    plt.plot(q, I_0)
    plt.title('Far Field of E = ' + str(E) + 'keV')
    plt.xlabel('$\\Delta q_x (nm^{-1})$')
    plt.ylabel('Relative Intensity')
    plt.axis('on')
    plt.show()
    print(q.shape)

QR(DetectorXValues, DetectorYValues, Detector, Distance2Detector)