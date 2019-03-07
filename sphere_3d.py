import meep as mp
import numpy as np
import matplotlib.pyplot as plt
import math

SphereRadius = 5

sx = 2*SphereRadius+20
sy = 1e3
# sy = SphereRadius*5
# sy = SphereRadius*2+20

sz = sy
# sy = 10
# sz = 10

wavelength = 1.67
f=1/wavelength
df = 0.1
resolution = 4
dpml=3


pml_layers= [mp.PML(dpml)]
cell = mp.Vector3(sx+2*dpml, sy+2*dpml, sz+2*dpml)
# symmetry = [mp.Mirror(mp.Y), mp.Mirror(mp.Z)]
# symmetry = [mp.Rotate4(direction=mp.X)]
symmetry = []



nr = 0.9999432732358271
ni = 6.7430E-06
real_epsilon = nr**2 - ni**2
imag_epsilon = 2*nr*ni
conductivity=2 * np.pi * f *imag_epsilon / (real_epsilon*8.85e-12)

# real_ind = 1 - 1.1676e-5
# imag_ind = 3.9765e-7
# real_epsilon = real_ind**2 - imag_ind**2
# imag_epsilon = 2*real_ind*imag_ind
# omega = 2*np.pi*f*1e10
# conductivity = omega*imag_epsilon
# factor = (1e-10/3e8)/(real_epsilon*8.85e-12)
# D_conductivity = factor*conductivity

mat = mp.Medium(epsilon=real_epsilon, D_conductivity=conductivity)

# geometry = [mp.Sphere(radius=SphereRadius, material=mat)]
geometry = []

# geometry = [mp.Block(mp.Vector3(SphereRadius,SphereRadius,SphereRadius),
#                      center=mp.Vector3(),
#                      material=mat)]



# sources = [mp.Source(mp.GaussianSource(f, fwidth=df),
#                      component=mp.Ez,
#                      center=mp.Vector3(-sx/2),
#                      size=mp.Vector3(y=sy, z=sz)),
#             mp.Source(mp.GaussianSource(f, fwidth=df),
#                      component=mp.Ey,
#                      center=mp.Vector3(-sx/2),
#                      size=mp.Vector3(y=sy, z=sz))]

sources = [mp.Source(mp.GaussianSource(f, fwidth=df),
                     component=mp.Ez,
                     center=mp.Vector3(-sx/2),
                     size=mp.Vector3(y=sy, z=sz)),
        mp.Source(mp.GaussianSource(f, fwidth=df),
                     component=mp.Ey,
                     center=mp.Vector3(-sx/2),
                     size=mp.Vector3(y=sy, z=sz))]
# k_point = mp.Vector3(f*real_ind).rotate(mp.Vector3(z=1), rot_angle)

# sources = [mp.EigenModeSource(src=mp.GaussianSource(f, fwidth=df),
#                               center=mp.Vector3(),
#                               size=mp.Vector3(y=sy, z=sz),
#                               direction=mp.AUTOMATIC if rot_angle == 0 else mp.NO_DIRECTION,
#                               eig_kpoint=k_point,
#                               eig_band=1,
#                               eig_parity=mp.EVEN_Y+mp.ODD_Z if rot_angle == 0 else mp.ODD_Z,
#                               eig_match_freq=True)]

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    symmetries=symmetry,
                    sources=sources,
                    resolution=resolution,
                    force_complex_fields=True,
                    Courant=0.49,
                    eps_averaging=False)

nearfield = sim.add_near2far(f, 0, 1,
        mp.Near2FarRegion(mp.Vector3(sx/2-9), size=mp.Vector3(y=sy, z=sz)))


sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(sx/2 - 3, 0), 1e-3))

Distance2Detector = 1e6
DetectorHeight = 5e5
DetectorResolutionY = DetectorHeight/10
DetectorWidth = wavelength*5
DetectorResolutionX = wavelength/10
DetectorResolutionZ = DetectorResolutionY


DetectorYValues = np.arange(-DetectorHeight/2, DetectorHeight/2 + DetectorResolutionY, DetectorResolutionY)
DetectorZValues = np.arange(-DetectorHeight/2, DetectorHeight/2 + DetectorResolutionZ, DetectorResolutionZ)

DetectorXValues = np.arange(-DetectorWidth/2, DetectorWidth/2 + DetectorResolutionX, DetectorResolutionX)
steps = DetectorXValues.shape[0]

FarErx = np.zeros((DetectorYValues.shape[0], DetectorXValues.shape[0], DetectorZValues.shape[0]))
FarEix = np.zeros((DetectorYValues.shape[0], DetectorXValues.shape[0], DetectorZValues.shape[0]))
FarEry = np.zeros((DetectorYValues.shape[0], DetectorXValues.shape[0], DetectorZValues.shape[0]))
FarEiy = np.zeros((DetectorYValues.shape[0], DetectorXValues.shape[0], DetectorZValues.shape[0]))
FarErz = np.zeros((DetectorYValues.shape[0], DetectorXValues.shape[0], DetectorZValues.shape[0]))
FarEiz = np.zeros((DetectorYValues.shape[0], DetectorXValues.shape[0], DetectorZValues.shape[0]))

for j, x in enumerate(DetectorXValues):
    if j % 5 == 0:
        print('Near2Far Field on step ' + str(j) + ' of ' + str(steps))
    for i, y in enumerate(DetectorYValues):
        for k,z  in enumerate(DetectorZValues):
            ff = sim.get_farfield(nearfield, mp.Vector3(Distance2Detector + x, y, z))
            FarErx[i, j, k] = ff[0].real
            FarEix[i, j, k] = ff[0].imag
            FarEry[i, j, k] = ff[1].real
            FarEiy[i, j, k] = ff[1].imag
            FarErz[i, j, k] = ff[2].real
            FarEiz[i, j, k] = ff[2].imag



FarIx = FarErx**2 + FarEix**2
FarIy = FarEry**2 + FarEiy**2
FarIz = FarErz**2 + FarEiz**2

Detector = FarIx + FarIy + FarIz

wavelengthnm = 0.167
Angles = np.arctan(DetectorYValues/Distance2Detector)
q_values = 4*np.pi*np.sin(Angles)/wavelengthnm
I = np.mean(Detector, axis=1)


plt.figure(dpi=100)

plt.imshow(I,interpolation='spline36', cmap='gist_heat', alpha=0.9, extent=[q_values[0], q_values[-1], q_values[0], q_values[-1]])

np.savetxt('I3dRadius_Ezy_pol%s.csv'%(SphereRadius),I, delimiter=',')

plt.axis('on')
plt.xlabel('qy (1/nm)')
plt.ylabel('qz (1/nm)')

plt.colorbar()
# plt.show()
plt.savefig('3d_image_ezy_002_r10.png')

