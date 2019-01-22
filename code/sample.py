import numpy as np
from scipy import constants as cst
import gf


#This class is currently on used to build a 2D environment

    #Define function that creates a block with a vacuum on top
    #Give the dim as an array of (x, y) in nm
    #Give the dimensions of the entire space
    #Currently assumes the relative mu=1
class Sample:

    def __init__(self, dim, vacuum_height, rel_eps, dx, cond, mcond):
        #Calculate epsilon not in SI units with nm and ns
        #There is not need to convert the conductivity to SI units with nm and ns because
        #you would multiply by 10^0 = 1
        self.eps0 = cst.value('electric constant')*(10**9)

        self.dx = dx
        self.dim = (dim/self.dx).astype(int)
        self.vh = vacuum_height
        self.cond = cond
        self.mcond = mcond


        #all these initial parameters are defined for our substrate, no nanostructures are taken into account in the initialization
        self._init(dim, vacuum_height, rel_eps)

    def _init(self, dim, vacuum_height, rel_eps):
        #Convert relative epsilon to epsilon
        eps = self.eps0*rel_eps

        #Create the space for the vacuum that stores epsilon as epsilon not
        eps_space_vacuum = np.ones((self.vh, self.dim[0]))*self.eps0

        #Create the space for the block that uses the given relative epsilon
        eps_space_substrate = np.ones((self.dim[1]-self.vh, self.dim[0]))*eps

        #Append the two epsilon spaces to create the total space
        self.eps_space = np.append(eps_space_vacuum, eps_space_substrate, axis = 0)

        #Create the conductivity space in the vacuum
        cond_space_vacuum = np.zeros((self.vh, self.dim[0]))

        #Create the conductivity space for the block using the given conductivty
        cond_space_substrate = np.ones((self.dim[1]-self.vh, self.dim[0]))*self.cond

        #Append the two conductivity spaces
        self.cond_space = np.append(cond_space_vacuum, cond_space_substrate, axis = 0)

        #The magnetic conductivity of the vacuum is the same as its electrical one (=0)

        #Create magnetic conductivity of the vacuum for the block
        mcond_space_substrate = np.ones((self.dim[1]-self.vh, self.dim[0]))*self.mcond

        #Append the two magnetic conductivity spaces
        self.mcond_space = np.append(cond_space_vacuum, mcond_space_substrate, axis = 0)

    #Define function that allows you to add a nanostructure
    #Currently just a rectangle
    #Give the center in terms of x
    def add_block(self, block_dim, center, rel_eps, cond, mcond):
        #Break up the space list into its individual parts
        eps_space = self.eps_space
        cond_space = self.cond_space
        mcond_space = self.mcond_space

        #Convert relative epsilon to epsilon
        eps = rel_eps*self.eps0

        #Convert given dimensions from nm to index
        ind_dim = (block_dim/self.dx).astype(int)

        #Create the block space
        eps_block = np.ones((ind_dim[1], ind_dim[0]))*eps

        #Get the dimensions of the entire space
        total_dim = np.shape(self.eps_space)

        #Convert the given center to an index
        ind_center = gf.coord2index(center, total_dim, self.dx)
        # print('center (idx)', ind_center)
        # print('block dim', ind_dim)
        # print('total dim', total_dim)

        #Determine the range of the added submatrix
        row_low = np.int(ind_center[0] - 1 - (ind_dim[0]/2))
        row_hi = np.int(row_low + ind_dim[0])

        col_low = np.int(ind_center[1] - 1 - (ind_dim[1]/2))
        col_hi = np.int(col_low + ind_dim[1])

        # print(row_low, row_hi)
        # print(col_low, col_hi)

        #Redefine epsilon space
        eps_space[row_low:row_hi, col_low:col_hi] = eps_block

        #Create the block's matrix for conductivity
        cond_block = np.ones((ind_dim[1], ind_dim[0]))*cond

        #Redefine conductivity space
        cond_space[row_low:row_hi, col_low:col_hi] = cond_block

        #Create the block's matrix for magentic conductivity
        mcond_block = np.ones((ind_dim[1], ind_dim[0]))*mcond

        #Redefine magnetic conductivity space
        mcond_space[row_low:row_hi, col_low:col_hi] = mcond_block

        self.eps_space = eps_space
        self.cond_space = cond_space
        self.mcond_space = mcond_space
