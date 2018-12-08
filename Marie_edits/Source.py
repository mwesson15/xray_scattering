import numpy as np
from scipy import constants as cst
from GeneralFunctions import general_functions as gf


class source: 
    #Define function that creates the incident wave space for the electric field
    #and calculate its values throughout time
    def incident_E(self, amp, k, omega, delta_x, alpha, delta_t, iterations,
                      vacuum_height, width):
        #Calculate and store an origin point that is far away
        r_x = -10000
        r_y = np.absolute(r_x*np.tan(alpha))

        #Create matrix to store electric field both incoming and out going
        E_in = np.array([])
        E_out = np.array([])
        
        #Create array of y values
        y_in = np.arange(0, vacuum_height, delta_x)
        y_min = (-1*width/2)*np.tan(alpha)
        y_out = np.arange(y_min, vacuum_height, delta_x)
        
        #Create a matrix that has all of the r_0-r values (R)
        dif_y_in = np.subtract(r_y, y_in)
        dif_y_out = np.subtract(r_y, y_out)
        dif_x_out = np.add(r_x, width)
        R_in = np.sqrt(np.add(r_x**2, np.square(dif_y_in)))
        R_out = np.sqrt(np.add(np.square(dif_x_out), np.square(dif_y_out)))
        
        #Initialize time
        t = 0
        
        #Go thru time and calculate incoming incident wave
        for i in range(iterations):
            #Calculate the current E in and E out
            current_E_in = gf.wave(amp, k, R_in, omega, t)
            current_E_out = gf.wave(amp, k, R_out, omega, t)
            
            #Append the current data to their respective matrices
            E_in = np.vstack(E_in, current_E_in)
            E_out = np.vstack(E_out, current_E_out)
            
            #Increment time
            t += delta_t
        
        #Transpose the two matrices to mke indexing more intuitive
        E_in = np.transpose(E_in)
        E_out = np.transpose(E_out)
        
        return [E_in, E_out]
    #The indexing of the E_in and E_out go [y, t]
    
    #Define function that does the same for H
    def incident_H(self, amp, k, omega, delta_x, alpha, delta_t, iterations,
                      vacuum_height, width):
        #Calculate and store an origin point that is far away
        r_x = -10000
        r_y = np.absolute(r_x*np.tan(alpha))

        #Create matrix to store electric field both incoming and out going
        H_in = np.array([])
        H_out = np.array([])
        
        #Create array of y values
        y_in = np.arange(0, vacuum_height, delta_x)
        y_in = np.add(y_in, delta_t/2)
        y_min = (-1*width/2)*np.tan(alpha)
        y_out = np.arange(y_min, vacuum_height, delta_x)
        y_out = np.add(y_out, delta_t/2)
        
        #Create a matrix that has all of the r_0-r values (R)
        dif_y_in = np.subtract(r_y, y_in)
        dif_y_out = np.subtract(r_y, y_out)
        dif_x_out = np.add(r_x, width)
        R_in = np.sqrt(np.add(r_x**2, np.square(dif_y_in)))
        R_out = np.sqrt(np.add(np.square(dif_x_out), np.square(dif_y_out)))
        
        #Initialize time
        t = delta_t/2
        
        #Go thru time and calculate incoming incident wave
        for i in range(iterations):
            #Calculate the current E in and E out
            current_H_in = gf.wave(amp, k, R_in, omega, t)
            current_H_out = gf.wave(amp, k, R_out, omega, t)
            
            #Append the current data to their respective matrices
            H_in = np.vstack(H_in, current_H_in)
            H_out = np.vstack(H_out, current_H_out)
            
            #Increment time
            t += delta_t
        
        #Transpose the two matrices to mke indexing more intuitive
        H_in = np.transpose(H_in)
        H_out = np.transpose(H_out)
        
        return [H_in, H_out]
    
    #Define function that gives the source throughout space
    def source_mtrx(self, E, H, dim_init, space):
        #Get some dimensions
        dim_tot = space.shape
        dim_E_in = E[0].shape
        dim_E_out = E[1].shape
        dim_H_in = H[0].shape
        dim_H_out = H[1].shape
        
