import meep as mp
import matplotlib.pyplot as plt
import numpy as np

dpml = 10
sx = 100
sy = 100
resolution = 3
f = 0.15

cell = mp.Vector3(sx+2*dpml, sy)

k_point = mp.Vector3()

pml_layers = [mp.PML(dpml, direction=mp.X)]

def pw_amp(k, x0):
    def _pw_amp(x):
        return np.exp(1j * k.dot(x + x0))
    return _pw_amp

alpha = 0
fcen = 0.15  # pulse center frequency
kdir = mp.Vector3(np.cos(alpha), np.sin(alpha))  # direction of k (length is irrelevant)
n = 1 # refractive index of material containing the source
k = kdir.unit().scale(2 * np.pi * fcen * n)  # k with correct length

geometry = [mp.Cylinder()]

sources = [
    mp.Source(
        mp.ContinuousSource(fcen),
        component=mp.Ez,
        center=mp.Vector3(-0.5 * sx, 0),
        size=mp.Vector3(0, sy),
        amp_func=pw_amp(k, mp.Vector3(x=-0.5 * sx)))]

sim = mp.Simulation(
    cell_size=cell,
    sources=sources,
    resolution=resolution,
    k_point=k_point,
    boundary_layers=pml_layers)

sim.run(until=350)

data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
data=data.real
plt.figure(dpi=100)
plt.imshow(data.transpose(), cmap='RdBu', alpha=0.9)
plt.axis('off')
plt.show()