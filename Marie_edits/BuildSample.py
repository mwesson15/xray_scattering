import numpy as np
from scipy import constants as cst
from GeneralFunctions import general_functions as gf

class build_sample:

    def __init__(self, dim, vacuum_height, relative_epsilon, delta_x, conductivity, conductivity_m):
        #all these initial parameters are defined for our substrate, no nanostructures are taken into account in the initialization
        self.create_substrate(dim, vacuum_height, relative_epsilon, delta_x, conductivity, conductivity_m)

#This class is currently on used to build a 2D environment

    #Define function that creates a block with a vacuum on top
    #Give the dim as an array of (x, y) in nm
    #Give the dimensions of the entire space
    #Currently assumes the relative mu=1
    def create_substrate(self, dim, vacuum_height, relative_epsilon, delta_x, conductivity, conductivity_m):
        #Calculate epsilon not in SI units with nm and ns
        eps_0 = cst.value('electric constant')*10^9
        
        #There is not need to convert the conductivity to SI units with nm and ns because
        #you would multiply by 10^0 = 1
        
        #Convert relative epsilon to epsilon
        eps = relative_epsilon*eps_0
        
        #Create the space for the vacuum that stores epsilon as epsilon not
        eps_space_vacuum = np.ones((vacuum_height, dim[0]))*eps_0
        
        #Create the space for the block that uses the given relative epsilon
        eps_space_substrate = np.ones((dim[1]-vacuum_height, dim[0]))*eps
        
        #Append the two epsilon spaces to create the total space
        eps_space = np.append(eps_space_vacuum, eps_space_substrate, axis = 0)
        
        #Create the conductivity space in the vacuum
        conductivity_space_vacuum = np.zeros((vacuum_height, dim[0]))
        
        #Create the conductivity space for the block using the given conductivty
        conductivity_space_substrate = np.ones((dim[1]-vacuum_height, dim[0]))*conductivity
        
        #Append the two conductivity spaces
        conductivity_space = np.append(conductivity_space_vacuum, conductivity_space_substrate, axis = 0)
        
        #The magnetic conductivity of the vacuum is the same as its electrical one (=0)
        
        #Create magnetic conductivity of the vacuum for the block
        conduct_m_space_substrate = np.ones((dim[1]-vacuum_height, dim[0]))*conductivity_m
        
        #Append the two magnetic conductivity spaces
        conduct_m_space = conductivity_space = np.append(conductivity_space_vacuum, conduct_m_space_substrate, axis = 0)
        
        return [eps_space, conductivity_space, conduct_m_space]
    
    #Define function that allows you to add a nanostructure
    #Currently just a rectangle
    #Give the center in terms of x
    def add_block(self, space, dim, center, block_height, rel_eps, conduct, conduct_m, delta_x):
        #Break up the space list into its individual parts
        eps_space = space[0]
        conduct_space = space[1]
        conduct_m_space = space[2]
        
        #Calculate epsilon not in SI units with nm and ns
        eps_0 = cst.value('electric constant')*10^9
        
        #Convert relative epsilon to epsilon
        eps = rel_eps*eps_0
        
        #Convert given dimensions from nm to index
        ind_dim = dim/delta_x
        
        #Calculate height of block in terms of index
        block_ind = block_height/delta_x
        
        #Create the block space
        eps_block = np.ones((ind_dim[1], ind_dim[0]))*eps
        
        #Get the dimensions of the entire space
        total_dim = np.shape(eps_space)
        
        #Convert the given center to an index
        ind_center = coord2index(center, total_dim, delta_x)
        
        #Determine the range of the added submatrix
        row_low = total_dim[0] - 1 - block_ind
        row_hi = row_low - ind_dim[0]
        
        column_low = ind_center - (ind_dim[1]/2)
        column_hi = colum_low + ind_dim[1]
        
        #Redefine epsilon space
        eps_space[row_hi:row_low, column_low:column_hi] = eps_block
        
        #Create the block's matrix for conductivity
        conduct_block = np.ones((ind_dim[1], ind_dim[0]))*conduct
        
        #Redefine conductivity space
        conduct_space[row_hi:row_low, column_low:column_hi] = conduct_block
        
        #Create the block's matrix for magentic conductivity
        conduct_m_block = np.ones((ind_dim[1], ind_dim[0]))*conduct_m
        
        #Redefine magnetic conductivity space
        conduct_m_space[row_hi:row_low, column_low:column_hi] = conduct_m_block
        
        return [eps_space, conduct_space, conduct_m_space]