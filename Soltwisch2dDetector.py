import numpy as np
import meep as mp
import matplotlib
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
            mp.Near2FarRegion(mp.Vector3(y=-sy/2), size=mp.Vector3(sx), weight=-1, direction=mp.Y),
            mp.Near2FarRegion(mp.Vector3(-sx/2,-dsub), size=mp.Vector3(y=sy-dsub), weight=-1, direction=mp.X),
            mp.Near2FarRegion(mp.Vector3(sx/2,-dsub), size=mp.Vector3(y=sy-dsub), weight=1, direction=mp.X))

    #run til things settle
    sim.run(until=5000)

    #calculate far field approximation through meep
    Distance2Detector = 1e10
    QxMin = 0.6
    QxMax = 1.6
    QyMin = -1
    QyMax = 1
    DetectorResolutionX = 200
    DetectorResolutionY = 200
    alpha = 0.86
    wavevector = 2*np.pi/0.167 #Current wave vector
    farI = farfield(sim, nearfield, Distance2Detector, QxMin, QxMax, QyMin, QyMax, DetectorResolutionX,
                    DetectorResolutionY, alpha, wavevector)

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
            mp.Near2FarRegion(mp.Vector3(y=-sy/2), size=mp.Vector3(sx), weight=-1, direction=mp.Y),
            mp.Near2FarRegion(mp.Vector3(-sx/2,-dsub), size=mp.Vector3(y=sy-dsub), weight=-1, direction=mp.X),
            mp.Near2FarRegion(mp.Vector3(sx/2,-dsub), size=mp.Vector3(y=sy-dsub), weight=1, direction=mp.X))

    #run til things settle
    sim.run(until=5000)

    farIVac = farfield(sim, nearfield, Distance2Detector, QxMin, QxMax, QyMin, QyMax, DetectorResolutionX,
                    DetectorResolutionY, alpha, wavevector)

    FarIDiff = farI - farIVac
    plt.imshow(FarIDiff, cmap='gist_heat', extent=[QyMin, QyMax, QxMin, QxMax])
    plt.title('Far Field Difference')
    plt.xlabel('$q_y (nm^{-1})$')
    plt.ylabel('$q_x (nm^{-1})$')
    plt.colorbar()
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


#Perofrms convolution
#Using our definition of x, y, z not the papers
#Q values given in inverse nm (like the paper)
def farfield(sim, nearfield, Distance2Detector, QxMin, QxMax, QyMin, QyMax,
             DetectorResolutionX, DetectorResolutionY, alpha, wavevector):
    #Create matrices of all the Qx and Qy values
    #Divide out the wavevector so that we only have the angular dependence (Need to multiply back in when plotting)
    #Factor of 10 is to convert to our MEEP units (angstroms)
    QxValues = np.linspace(QxMin, QxMax, DetectorResolutionX)/(10*wavevector)
    QyValues = np.linspace(QyMin, QyMax, DetectorResolutionY)/(10*wavevector)

    #We can now solve for all the y values of the far field as qy=k(sin(alpha_f)+sin(alpha_i))
    #So y = tan(arcsin(qy/k - sin(alpha)))*Distance2Detector
    YValues = Distance2Detector*np.tan(np.arcsin(QyValues-np.sin(alpha)))

    #The x values can now be calculated
    #qx=k*sin(theta)*cos(alpha_f)
    #So x=Distance2Detector*tan(arcsin(qx/(k*cos(arctan(y/Distance2Detector)))))
    XValues = Distance2Detector*np.tan(np.arcsin(
        np.divide(QxValues, np.cos(np.arctan(YValues/Distance2Detector)))))

    #We know make an empty array that will store the Far Field
    FarE = np.zeros((YValues.shape[0], XValues.shape[0]))

    #This is used just to track progress for convenience
    steps = YValues.shape[0]

    #We now go thru all of the x and y values and perform the convolution
    for i, y in enumerate(YValues):
        if i % 5 == 0:
            print('Convolution is on step ' + str(i) + ' of ' + str(steps))
        for j, x in enumerate(XValues):
            ff = sim.get_farfield(nearfield, mp.Vector3(x, y, Distance2Detector))
            FarE[i, j] = np.abs(ff[0] + ff[1] + ff[2]) #Storing amplitude of complex field

    #Calculate intensity (ignoring constants since we may renormalize)
    FarI = np.square(FarE)

    #Put back in the factors of q to plot nicely
    QyValues = QyValues*wavevector*10
    QxValues = QxValues*wavevector*10

    #Plot the data
    plt.imshow(FarI, cmap='gist_heat', extent=[QyValues.min(), QyValues.max(), QxValues.min(), QxValues.max()])
    plt.title('Far Field')
    plt.xlabel('$q_y (nm^{-1})$')
    plt.ylabel('$q_x (nm^{-1})$')
    plt.colorbar()
    plt.show()

    return FarI

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