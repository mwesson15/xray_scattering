import meep as mp
import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
import cmath
import math
from materials import Gold, Silicon, Materials

class globalVariables():
    cellXSize = 100
    cellYSize = 100
    epsilonOfMaterial = 1
    indexOfSphere = 0.9982
    # indexOfSphere = 3.4

    sphereRadius = 20
    resolution = 10
    projectName = "sphere_scattering"
    h5Name = "ez"
    imageStorageDirectory = projectName +'-out'

def generateSphere(radius):
    si = Gold()
    ind = si.index_refract()
    ind = ind.real
    cond = ind.imag
    wv = 0.0167
    freq = 1/wv
    mat = mp.Medium(index=ind, D_conductivity=2*math.pi*freq*cond/ind)
    # sphere = mp.Cylinder(radius=radius, material=mp.Medium(index=globalVariables.indexOfSphere))
    sphere = mp.Cylinder(radius=radius, material=mat)

    return sphere

def createSimulationSpace():
    cell = mp.Vector3(globalVariables.cellXSize, globalVariables.cellYSize, 0)

    geometry = [generateSphere(globalVariables.sphereRadius)]


    sources = generateSource()

    pml_layers = [mp.PML(1.0)]

    return cell,geometry,sources,pml_layers

def pw_amp(k,x0):
    def _pw_amp(x):
        return cmath.exp(1j*2*math.pi*k.dot(x+x0))
    return _pw_amp

def generateSource():

    use_cw_solver = False
    src_pt = mp.Vector3(-49,0,0)
    ng = 1
    fcen = 0.15
    wv = 0.0167
    fcen = 1/wv
    df = 0.1
    theta_in = 0
    sy = 100
    k = mp.Vector3(math.cos(theta_in),math.sin(theta_in),0).scale(fcen*ng)
    sources = [mp.Source(mp.ContinuousSource(fcen,fwidth=df) if use_cw_solver else mp.GaussianSource(fcen,fwidth=df),
                     component=mp.Ez,
                     center=src_pt,
                     size=mp.Vector3(0,sy,0),
                     amp_func=pw_amp(k,src_pt))]

    # sources = [mp.Source(mp.ContinuousSource(frequency=0.15),
    #                  component=mp.Ez,
    #                  center=mp.Vector3(0, 0),
    #                  size=mp.Vector3(50, 10, 0),
    #                  amp_func=amp),
    #        mp.Source(mp.CustomSource(src_func=time_dependence),
    #                  component=mp.Ez,
    #                  center=mp.Vector3(0, 0),
    #                  size=mp.Vector3(50, 10, 0),
    #                  amp_func=amp)]
    # sources = [mp.Source(mp.ContinuousSource(frequency=0.15),
    #                  component=mp.Ez,
    #                  center=mp.Vector3(0, 0),
    #                  size=mp.Vector3(50, 10, 0),
    #                  amp_func=amp)]
    return sources


def amp(y):
    z = y[0]
    r = y[1]
    z_R = 3*np.pi*100/20
    if z != 0:
        R = z*(1 + (z_R/z)**2)
    else:
        R=1e20
    A = (1 + (z/z_R)**2)**0.5
    w = 10*A
    k = 3*np.pi/10
    amp = np.exp(-1*(r/w)**2)/A
    phase = (k*z) + (k*(r**2)/(2*R)) - np.arctan(z/z_R)
    return amp*np.exp(1j * phase)

def time_dependence(t):
    freq = 1/0.0167
    return np.exp(1j*2*np.pi*freq*t)



def plot_data(sim,cell):
    eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
    #plt.figure(dpi=100)
    #plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
    #plt.axis('off')
    #plt.show()

    ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
    plt.figure(dpi=100)
    plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
    plt.imshow(ez_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
    plt.axis('off')
    plt.savefig("spheresimulation.png")


def pngToGIF():
    gifCreationString = "convert " + globalVariables.imageStorageDirectory + "/" + globalVariables.projectName + "-" + globalVariables.h5Name + "-*.png " + globalVariables.imageStorageDirectory + "/" + globalVariables.projectName + ".gif"
    print(gifCreationString)
    os.system(gifCreationString)

def compressGIF(scalingFactor):
    print("Did I get here for GIF compression?")
    gifCompressionString = "gifsicle -O3 --colors=128 --scale=" + str(scalingFactor) + " -i " + globalVariables.imageStorageDirectory + "/" + globalVariables.projectName + ".gif" + " -o " + globalVariables.imageStorageDirectory + "/" + globalVariables.projectName + "compressed" + ".gif"
    print(gifCompressionString)
    os.system(gifCompressionString)

def countNumberOfPNG():
    pngCountString = "ls " + globalVariables.imageStorageDirectory + "/" + "*.png | wc -l"
    count = os.system(pngCountString)
    return int(count)

def deletePNG():
    print("Number of PNG Files in directory before is :" + str(countNumberOfPNG()))
    deletePNGString = "rm " + globalVariables.imageStorageDirectory + "/" + globalVariables.projectName + "*.png"
    os.system(deletePNGString)
    print("Number of PNG Files in directory after is :" + str(countNumberOfPNG()))

def deleteCompressedGIF():
    pass


'''
Creates a gif from the hdf5 file.
'''
def inPlaceGifCreation(compressGIFBool=True, compressGIFScale=0.5 ,deletePNGBool=True,generateGIF=True,deleteUncompressedGIF=False):
    pngToGIF()

    if(compressGIFBool==True):
        compressGIF(scalingFactor=compressGIFScale)

    if(deletePNGBool==True):
        deletePNG()


def removeH5Files():
    pass

if __name__=="__main__":
    #cell, geometry, sources, pml_layers = createDiffractionSlit()
    cell , geometry , sources , pml_layers = createSimulationSpace()
    resolution = globalVariables.resolution
    sym = [mp.Mirror(mp.Y)]
    sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    symmetries=sym,
                    sources=sources,
                    resolution=resolution)
    sim.use_output_directory()
    #sim.run(mp.at_beginning(mp.output_epsilon),
       # mp.to_appended("ez", mp.at_every(0.6, mp.output_png(mp.Ez,"-Zc dkbluered"))),
       # until=200)
    sim.run(mp.at_every(1 , mp.output_png(mp.Ez, "-Zc dkbluered")), until_after_sources=70)
    # sim.run(mp.at_every(0.6 , mp.output_png(mp.Ez, "-Zc dkbluered")), until=20)

    plot_data(sim,cell)
    inPlaceGifCreation()