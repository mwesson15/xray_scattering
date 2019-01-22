from scipy import constants
import numpy as np
# import meep as mp
c = constants.value(u'speed of light in vacuum')
r_0 = constants.value(u'classical electron radius') #meters


class Material(object):

    def __init__(self, properties):
        element = properties[5]
        print('Initializing %s parameters'%element)
        self.element = element
        self.f1 = properties[0]
        self.f2 = properties[1]
        self.density = properties[2]*1e6
        mass = properties[3]
        self.conductivity = properties[4]
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
    def conductivity(self):
        return self.conductivity

class Gold(Material):
    #right now all these values are at 7350 eV
    def __init__(self):
        f1 = 75.31
        f2 = 8.945
        density = 19.3 #g/cm^3
        mass = 196.967
        name = 'gold'
        # resistivity = 22.14 #ohms*nm
        conductivity = 4.5e7 #siemens/m
        properties = [f1, f2, density, mass, conductivity, name]

        super().__init__(properties)

class Silicon:
    #right now all these values are at 7350 eV
    def __init__(self):
        f1 = 14.28
        f2 = 0.3864
        density = 2.33 #g/cm^3
        mass = 28.086
        conductivity = 1000 #siemens/m
        name = 'silicon'

        properties = [f1, f2, density, mass, conductivity, name]
        super().__init__(properties)
