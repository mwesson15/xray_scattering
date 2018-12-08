
from scipy import constants
import numpy as np 
# import meep as mp
c = constants.value(u'speed of light in vacuum')
r_0 = constants.value(u'classical electron radius') #meters


class Materials(object):
    
    def __init__(self, f1, f2, density, mass, element):
        print('Initializing %s parameters'%element)
        self.element = element
        self.f2 = f2
        self.f1 = f1
        self.density = density*1e6
        self.na = self.density*constants.value(u'Avogadro constant')/mass

    def index_refract(self, energy=7350):
        n = 1 - self.dispersive() + 1j*self.absorptive()
        if self.element is 'vacuum':
            n = 1
        return n

    def dispersive(self, energy=7350):
        delta = self.na*r_0*(1240/energy)**2*1e-18*self.f1/(2*np.pi)
        return delta
        
    def absorptive(self, energy=7350):
        beta = self.na*r_0*(1240/energy)**2*1e-18*self.f2/(2*np.pi)
        return beta