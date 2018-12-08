import numpy as np
from scipy import constants as cst

class general_functions:
    #Define function that normalizes vector
    def normalize(self, vector):
        #Calculate magnitude of vector
        amp = (np.dot(vector, vector))**0.5
        
        #Normalize
        normalized_vector = vector/amp
        
        return normalized_vector
    
    #Define function that converts E to E-tilda which has a similar magnitude as H
    #Unsure if the is needed at this time
    #Will need to multiply D by c in order to do the same
    def reduce_E(self, E):
        #Calculate factor to multiply it by
        factor = (cst.value('electric constant')/cst.value('mag. constant'))**0.5
        
        #Multiply E by the factor
        new_E = np.multiply(factor, E)
        
        return new_E
    
    #Define function that converts from index to x, y, or z coordinate
    #Should work with an array index (x,y,z), but I've only tested with 1 value at a time
    def index2coord(self, index, space, delta_x):
        #Get the dimensions of space
        dim = space.shape
        
        #Calculate index of the origin
        origin_index = np.divide(np.subtract(dim, 1), 2)
        
        #Calculate distance between given index and origin
        index_spacing = np.subtract(index, origin_index)
    
        #Convert index_spacing into nanometers
        coord = np.multiply(delta_x, index_spacing)

        return coord, index_spacing
    
    #Define function that converts from (x,y,z) coordinates to index
    #Should work with an array index (x,y,z), but I've only tested with 1 value at a time
    def coord2index(self, coordinates, space, delta_x):
        #Get the dimensions of space
        dim = space.shape
        
        #Convert from nm to index distance
        index_spacing = np.divide(coordinates, delta_x)
        
        #Calculate index of the origin
        origin_index = np.divide(np.subtract(dim, 1), 2)
        
        #Calculate the index
        index = np.add(origin_index, index_spacing)
        
        return index
    
    #Define function that calculates the curl of a vector field at all points
    #Will we need to use different function when calculating curl at boundary?
    #Might have to check equations (ensured proper shifts)
    #Different shift for H?
    def curl(self, field_space, delta_x):
        #Assumes cubic Yee Cells
        #Used for 2D simulation
        
        #Get dimensions of space
        dim = field_space.shape
        
        #Create matrices the are shifted +1 indices in the x or y direction
        #This will allow us to calculate the curl of all the points at once in a matrix
        
        #Split off the last column of the array and store both subarrays
        split_column = np.hsplit(field_space, [dim[1]-1])
        
        #Append the two back together to get the x-shifted matrix
        shifted_x = np.append(split_column[1], split_column[0], axis = 1)
        
        #Repeat the above process to make the y-shifted matrix
        split_row = np.vsplit(field_space, [dim[0]-1])
        shifted_y = np.append(split_row[1], split_row[0], axis = 0)
        
        #Calculate x component
        curl_x = (shifted_y[:, :, 2] - field_space[:, :, 2])/delta_x
        
        #Calculate the y component
        curl_y = (field_space[:, :, 2] - shifted_x[:, :, 2])/delta_x
        
        #Calculate the z component by first calculating the 2 derivatives needed to calculate it
        dEy_dx = (shifted_x[:, :, 1] - field_space[:, :, 1])/delta_x
        dEx_dy = (shifted_y[:, :, 0] - field_space[:, :, 0])/delta_x
        curl_z = dEy_dx - dEx_dy
        
        #Collect all three components into one 3D matrix
        curl = np.concatenate((curl_x, curl_y, curl_z), axis = 2)
        
        return curl
    
    #Calculate the curl of a vector field at 1 point
    #Need to check indices
    def curl_point(self, field_space, index, delta_x):
        #Calculate x component
        curl_x = (field_space[index[0], index[1]-1, 2] - field_space[index[0], index[1], 2])/delta_x
        
        #Calculate the y component
        curl_y = (field_space[index[0], index[1], 2] - field_space[index[0]-1, index[1], 2])/delta_x
        
        #Calculate the z component by first calculating the 2 derivatives needed to calculate it
        dHy_dx = (field_space[index[0]-1, index[1], 1] - field_space[index[0], index[1], 1])/delta_x
        dHx_dy = (field_space[index[0], index[1]-1, 0] - field_space[index[0], index[1], 0])/delta_x
        curl_z = dHy_dx - dHx_dy
        
        #The +1 in the indices causes you to reference the next lower x, y, or z value
        
        #Collect all three components into one vector of (x, y, z)
        curl = np.array([curl_x, curl_y, curl_z])
        
        return curl
    
    #Due to the way we have defined the curl, the collection of H values with the same index as a
    #set of E values that in space are shifted in the +x, -y direction
    
    #Define function that calculates only the phase of the wave at a certain point
    def phase(self, k, r, r_0, omega, t):
        #r_0 = origin of the wave, ideally (x,y) = (-big number, xtan(alpha_i))
        #Calculate distance between r and r_0
        c = r - r_0 #Gives difference between points as a vector
        R = (c[0]**2 + c[1]**2)**0.5
        
        #Calculate the phase
        phi = k*R - omega*t
        
        return phi
    
    #Define function that calculates k from the wavelength (nm)
    #The equation allows for this one function to convert between the 2 in either direction
    def lambda2k(self, wavelength):
        #Plug into equation
        k = 2*np.pi/wavelength
        
        return k
    
    #Define function that converts between angular frequency and frequency
    def omega2f(self, omega):
        #Plug into equation
        f = 2*np.pi/omega
        
        return f
    
    #Define function that converts between frequency and angular frequency
    def f2omega(self, f):
        #Plug into equation
        omega = 2*np.pi*f
        
        return omega
    
    #Define function that returns generic wave equation (E=E_0*sin[k(r-r_0)-omega*t])
    #Ignore initial phase
    #Accepts r and r_0 in nm
    def wave(self, amplitude, k, r, r_0, omega, t):
        #Calculate distance between r and r_0
        c = r - r_0 #Gives difference between points as a vector
        R = (c[0]**2 + c[1]**2)**0.5
        
        #Plug into equation
        wave = amplitude*sin(k*R - omega*t)
        
        return wave
