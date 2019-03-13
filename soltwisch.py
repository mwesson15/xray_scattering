import numpy as np
import meep as mp
import matplotlib
from tqdm import tqdm
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

'''this code is built to match the publication:
Reconstructing detailed line profiles of lamellar gratings
from GISAXS patterns with a Maxwell solver (Soltwisch et. al. 2017)'''

#main simulation loop
#all distances below are expressed in units of 10^-10m
#in this simulation, x is horizontal, -y is up and z is into the screen
#unlike the paper convention of y is horizontal, z is up and x is into the screen
def run():
    #G1 params in 10^-10m
    gp = 1000.             # grating period
    gh = 1027.1               # grating height
    lw = 540.4              #line width
    tr = 162.7              #top radius
    br = 157.9             #bottom radius
    sa = 90.91              #sidewall angle (deg)
    #G2 params in 10^-10m
    gp = 1500.             # grating period
    gh = 1195.0               # grating height
    lw = 673.0              #line width
    tr = 91.6              #top radius
    br = 130.2             #bottom radius
    sa = 84.73              #sidewall angle (deg)

    #params for spacing above and below grating
    dsub = 300
    dpad = 300
    dpml = 50
    resolution = 4.5/20
    S = 0.49 #courant parameter
    alpha = 0.86*np.pi/180 #from Fig. 4, in radians
    k_point = mp.Vector3(0, 0, 0)

    sx = 2*gp #simulation width
    sy = dsub+gh+dpad #simulation height
    cell = mp.Vector3(sx,sy+2*dpml)

    #only absorb in Y direction so that we can have periodic boundaries in X
    pml_layers = [mp.PML(thickness=dpml,direction=mp.Y)]

    #params for light
    e = 0
    Energies = [5.5, 5.55, 5.6]
    wavelengths = [2.25, 2.23, 2.21]
    #corresponding params for Si at these energies
    delta = [1.1676e-5, 1.6068e-5, 1.5779e-5]
    beta = [3.9765e-7, 7.3314e-7, 7.082e-7]
    E = Energies[e]
    wvl = wavelengths[e]
    f = 1/wvl

    src = build_src(sx=sx, sy=sy, f=f*np.sin(alpha), alpha=alpha) #use frequency in 2d slice of infinite plane
    mat = build_mat(delta=delta[e], beta=beta[e], f=f) #use base frequency to get material params right
    geom = build_geom(sx=sx, sy=sy, dsub=dsub, gp=gp, gh=gh, lw=lw, tr=tr, br=br, sa=sa, material=mat)

    #first run empty sim for the background
    sim = mp.Simulation(cell_size=cell,
                        resolution=resolution,
                        Courant=S,
                        k_point=k_point,
                        boundary_layers=pml_layers,
                        geometry=[],
                        sources=src,
                        force_complex_fields=True,
                        eps_averaging=True)

    #nearfield needs to be specified in a region of homogenous space
    nearfield = sim.add_near2far(f, 0, 1,
            mp.Near2FarRegion(mp.Vector3(y=-sy/2+100), size=mp.Vector3(sx), weight=1, direction=mp.Z))

    #run til things settle
    sim.run(until=5000)

    _eps, _I = localfield(sim=sim, sx=sx, sy=sy)

    _inten = _I/np.amax(_I)
    plt.imshow(_eps.transpose(), cmap='binary')
    plt.imshow(_inten.transpose(), cmap='gist_heat', alpha=0.9)
    plt.colorbar()
    plt.show()

    #calculate far field approximation through meep
    d = 1.7e10 #min detector distance in paper
    steps = 50
    neighbors = 4
    _q, _farI = farfield(sim=sim, nearfield=nearfield, sx=sx, alpha=alpha, wavelength=wvl, d=d, steps=steps, neighbors=neighbors)

    #plot far field (compare to Fig. 6)
    plt.plot(_q, _farI/np.amax(_farI))
    plt.title('Far Field of E = ' + str(E) + 'keV')
    plt.xlabel('$\\Delta q_x (nm^{-1})$')
    plt.ylabel('Relative Intensity')
    plt.axis('on')
    plt.xlim(-0.7,0.7)
    plt.show()

    #now run the same simulation with the gratings
    sim.reset_meep()

    sim = mp.Simulation(cell_size=cell,
                        resolution=resolution,
                        Courant=S,
                        k_point=k_point,
                        boundary_layers=pml_layers,
                        geometry=geom,
                        sources=src,
                        force_complex_fields=True,
                        eps_averaging=True)

    #nearfield needs to be specified in a region of homogenous space
    nearfield = sim.add_near2far(f, 0, 1,
            mp.Near2FarRegion(mp.Vector3(y=-sy/2+100), size=mp.Vector3(sx), weight=1, direction=mp.Z))

    #run til things settle
    sim.run(until=5000)

    #calculate near field
    eps, I = localfield(sim=sim, sx=sx, sy=sy)

    #plot near field (compare to Fig. 4)
    inten = I/np.amax(I)
    plt.imshow(eps.transpose(), cmap='binary')
    plt.imshow(inten.transpose(), cmap='gist_heat', alpha=0.9)
    plt.colorbar()
    plt.show()

    #calculate far field approximation through meep
    d = 1.7e10 #min detector distance in paper
    steps = 50
    neighbors = 4
    q, farI = farfield(sim=sim, nearfield=nearfield, sx=sx, alpha=alpha, wavelength=wvl, d=d, steps=steps, neighbors=neighbors)

    #plot far field (compare to Fig. 6)
    plt.plot(q, farI/np.amax(farI))
    plt.title('Far Field of E = ' + str(E) + 'keV')
    plt.xlabel('$\\Delta q_x (nm^{-1})$')
    plt.ylabel('Relative Intensity')
    plt.axis('on')
    plt.xlim(-0.7,0.7)
    plt.show()

    plt.plot(q, (farI-_farI)/np.amax(farI-_farI))
    plt.title('Far Field of E = ' + str(E) + 'keV')
    plt.xlabel('$\\Delta q_x (nm^{-1})$')
    plt.ylabel('Relative Intensity')
    plt.axis('on')
    plt.xlim(-0.7,0.7)
    plt.show()


#calculates local field intensity for display (Fig. 4)
#accepts:
#   sim - the current simulation object (meep.Simulation)
#   sx - cell size in x (float)
#   sy - cell size in y (float)
#returns:
#   eps - array of epsilon values in the cell (numpy.array of shape (Nx,Ny))
#   I - array of intensity values in the cell (numpy.array of shape (Nx,Ny))
def localfield(sim, sx, sy):
    #intensity calculation from components
    Ex = sim.get_array(center=mp.Vector3(), size=mp.Vector3(sx, sy), component=mp.Ex, cmplx=True)
    Ey = sim.get_array(center=mp.Vector3(), size=mp.Vector3(sx, sy), component=mp.Ey, cmplx=True)
    Ez = sim.get_array(center=mp.Vector3(), size=mp.Vector3(sx, sy), component=mp.Ez, cmplx=True)
    Ex = np.abs(Ex)
    Ey = np.abs(Ey)
    Ez = np.abs(Ez)
    I = np.sqrt(np.square(Ey) + np.square(Ex) + np.square(Ez))
    eps = sim.get_array(center=mp.Vector3(), size=mp.Vector3(sx, sy), component=mp.Dielectric)
    return eps, I


#performs the farfield convolution for our simulation
#accepts:
#   sim - the current simulation object (meep.Simulation)
#   nearfield - a near field region to compute the farfield with (meep.Near2FarRegion)
#   sx - cell size in x (float)
#   alpha - angle of incident light (float)
#   wavelength - wavelength of light (float)
#   d - distance to detector in z direction (float)
#   steps - resolution for far field (int)
#   neighbors - number of grating neighbors to one side (int)
#returns:
#   q - the qy array for generating plots like Fig. 6 (numpy.array of shape (N,1))
#   farI - the farfield intensity at the different qy values (numpy.array (N,1))
def farfield(sim, nearfield, sx, alpha, wavelength, d, steps, neighbors):
    #build the arrays needed for the qz=0 plane (qx=0 in the paper)
    r = d*np.tan(alpha)
    z = 0
    thet = np.linspace(-np.pi/2, np.pi/2, 2*steps+1)
    alph = np.concatenate((np.arange(0,np.pi/2,np.pi/2/steps),np.linspace(np.pi/2,0,steps+1)))

    #here we calculate the far field effect of a single grating on its image straight ahead, and 8 neighbors,
    #we then add in order to get an approximate far field for a single grating periodic in the x direction
    farE = np.zeros((1+2*neighbors,thet.shape[0],3), dtype=np.complex128)

    #far field calculation with meeps built in convolution
    for i in tqdm(range(neighbors+1)):
        for j in tqdm(range(alph.shape[0])):
            if i == 0:
                ff = sim.get_farfield(nearfield, mp.Vector3(r*np.sin(thet[j]),-r*np.sin(alph[j]),z))
                farE[i,j,:] = ff[0:3]
            else:
                ff = sim.get_farfield(nearfield, mp.Vector3(r*np.sin(thet[j])+i*sx,-r*np.sin(alph[j]),z))
                farE[i,j,:] = ff[0:3]
                ff = sim.get_farfield(nearfield, mp.Vector3(r*np.sin(thet[j])-i*sx,-r*np.sin(alph[j]),z))
                farE[i+1,j,:] = ff[0:3]

    #collapse far field neighbors into single array
    farE = np.abs(farE)
    farI = np.sum(farE,axis=2)
    farI = np.sum(farI,axis=0)
    farI = farI.astype(np.single)

    #build qx (qy in paper) array according to equations in paper
    q = (2*np.pi/wavelength)*np.sin(thet)*np.cos(alph)/10

    return q, farI


#builds the source representation in meep
#accepts:
#   sx - cell size in x (float)
#   sy - cell size in y (float)
#   f - frequency of light (float)
#   alpha - angle of incident light (float)
#returns:
#   src - a list of meep.Source classes
def build_src(sx, sy, f, alpha):
    width = 1000
    df = 1/width

    #use plane waves as the incident beam spot is large relative to structure
    #both a Ez and Hz source are used here to generate all fields (unpolarized)
    src = [
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
            amplitude=1)
    ]

    return src


#builds the material representation for a given frequency
#accepts:
#   delta - delta for calculating refractive index (float)
#   beta - beta for calculating refractive index (float)
#   f - frequency of light (float)
#returns:
#   material - a meep representation of material properties (meep.Medium)
def build_mat(delta, beta, f):
    real_ind = 1 - delta
    imag_ind = beta
    real_epsilon = real_ind**2 - imag_ind**2
    imag_epsilon = 2*real_ind*imag_ind
    omega = 2*np.pi*f*1e10
    conductivity = omega*imag_epsilon

    #meep is not entirely unit invariant so this factor is needed
    factor = (1e-10/3e8)/(real_epsilon*8.85e-12)
    D_conductivity = factor*conductivity

    material = mp.Medium(epsilon=real_epsilon, D_conductivity=D_conductivity)

    return material


#builds 2 periods of grating geometries for the simulation according to the paper parameters
#accepts:
#   sx - cell size in x (float)
#   sy - cell size in y (float)
#   dsub - height of substrate in y (float)
#   gp - grating period (float)
#   gh - grating height (float)
#   lw - line width (float)
#   tr - top radius (float)
#   br - bottom radius (float)
#   sa - side angle (float)
#   material - material representation in meep (meep.Medium)
#returns:
#   geometry - a list of meep geometries
def build_geom(sx, sy, dsub, gp, gh, lw, tr, br, sa, material):
    tw = gh/np.tan(sa*2*np.pi/360) # part to subtract from top of grating width
    a = (gh-tr)/np.tan(sa*2*np.pi/360) # intermediate for vertex calculation

    #vertices for trapezoids approximating grating cross-section
    vtx = [mp.Vector3(-0.5*lw-0.5*a+gp/2,-1*(-0.5*gh),0),
                mp.Vector3(0.5*lw+0.5*a+gp/2,-1*(-0.5*gh),0),
                mp.Vector3(0.5*lw-0.5*a+gp/2,-1*(0.5*gh-tr),0),
                mp.Vector3(-0.5*lw+0.5*a+gp/2,-1*(0.5*gh-tr),0)]

    vtx2 = [mp.Vector3(-0.5*lw-0.5*a-gp/2,-1*(-0.5*gh),0),
                mp.Vector3(0.5*lw+0.5*a-gp/2,-1*(-0.5*gh),0),
                mp.Vector3(0.5*lw-0.5*a-gp/2,-1*(0.5*gh-tr),0),
                mp.Vector3(-0.5*lw+0.5*a-gp/2,-1*(0.5*gh-tr),0)]

    #rounded corners for top of trapezoids
    c1 = mp.Cylinder(radius=tr,
                height=mp.inf,
                axis=mp.Vector3(0,0,1),
                center=mp.Vector3(-gp/2-lw/2+tw/2+tr,-1*(-sy/2+dsub+gh-tr), 0),
                material=material)
    c2 = mp.Cylinder(radius=tr,
                height=mp.inf,
                axis=mp.Vector3(0,0,1),
                center=mp.Vector3(-gp/2+lw/2-tw/2-tr,-1*(-sy/2+dsub+gh-tr), 0),
                material=material)
    c3 = mp.Cylinder(radius=tr,
                height=mp.inf,
                axis=mp.Vector3(0,0,1),
                center=mp.Vector3(gp/2-lw/2+tw/2+tr,-1*(-sy/2+dsub+gh-tr), 0),
                material=material)
    c4 = mp.Cylinder(radius=tr,
                height=mp.inf,
                axis=mp.Vector3(0,0,1),
                center=mp.Vector3(gp/2+lw/2-tw/2-tr,-1*(-sy/2+dsub+gh-tr), 0),
                material=material)

    #blocks for top of trapezoids inbetween rounded corners
    b1 = mp.Block(center=mp.Vector3(-gp/2,-1*(-sy/2+dsub+gh-tr/2)), size=mp.Vector3(lw-tw-2*tr,tr,mp.inf), material=material)
    b2 = mp.Block(center=mp.Vector3(gp/2,-1*(-sy/2+dsub+gh-tr/2)), size=mp.Vector3(lw-tw-2*tr,tr,mp.inf), material=material)

    #ellipsoid cutout to make bottom of grating round
    e1 = mp.Ellipsoid(center=mp.Vector3(0,-1*(-sy/2+dsub),0), size=mp.Vector3(gp-lw-a,br,mp.inf), material=mp.Medium(epsilon=1))
    e2 = mp.Ellipsoid(center=mp.Vector3(-gp,-1*(-sy/2+dsub),0), size=mp.Vector3(gp-lw-a,br,mp.inf), material=mp.Medium(epsilon=1))
    e3 = mp.Ellipsoid(center=mp.Vector3(gp,-1*(-sy/2+dsub),0), size=mp.Vector3(gp-lw-a,br,mp.inf), material=mp.Medium(epsilon=1))

    geometry = [mp.Block(material=material, size=mp.Vector3(sx,dsub,mp.inf), center=mp.Vector3(0,-1*(-0.5*sy+0.5*dsub),0)),
                mp.Prism(vtx, height=mp.inf, center=mp.Vector3(0,-1*(-0.5*sy+dsub+0.5*gh),0), material=material),
                mp.Prism(vtx2, height=mp.inf, center=mp.Vector3(0,-1*(-0.5*sy+dsub+0.5*gh),0), material=material),
                c1, c2, c3, c4, b1, b2, e1, e2, e3]

    return geometry


if __name__=='__main__':
    run()
