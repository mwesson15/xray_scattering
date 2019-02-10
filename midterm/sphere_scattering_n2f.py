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
    # indexOfSphere = 3.4 #visible
    sphereRadius = 20
    resolution = 10
    projectName = "sphere_scattering"
    h5Name = "ez"
    imageStorageDirectory = projectName +'-out'

def generateSphere(radius):
    si = Gold()
    ind = si.index_refract()
    ind = ind.real
    # ind = 0.8 #we saw scattering here
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
    return sources

# def create_NearFieldFlux():



def plot_data(sim,cell):
    eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
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
    cell , geometry , sources , pml_layers = createSimulationSpace()
    resolution = globalVariables.resolution
    sym = [mp.Mirror(mp.Y)]
    sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    symmetries=sym,
                    sources=sources,
                    resolution=resolution)
    d1 = 15
    fcen = 1/0.0167
    nearfield = sim.add_near2far(
    fcen, 0, 1,
    mp.Near2FarRegion(mp.Vector3(0, 10), size=mp.Vector3(d1)),
    mp.Near2FarRegion(mp.Vector3(d1, 10), size=mp.Vector3(0, -20), weight=-1.0),
    mp.Near2FarRegion(mp.Vector3(0, -10), size=mp.Vector3(d1))
    )

    d2 = 15*1e3
    h = 50
    sim.use_output_directory()
    sim.run(mp.at_every(1 , mp.output_png(mp.Ez, "-Zc dkbluered")), until_after_sources=70)
    sim.output_farfields(nearfield, "spectra-{}-{}-{}".format(d1, d2, h),
                     mp.Volume(mp.Vector3(0,  d2 + (0.5 * h)), size=mp.Vector3(1e4, h)),
                     resolution)

    plot_data(sim,cell)
    inPlaceGifCreation()