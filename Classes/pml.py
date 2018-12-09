import numpy as np
from general_functions import general_functions as gf

class pml:
    def add_pml(space, d, m, R, delta_x):
        #Define epsilon_0
        eps_0 = cst.value('electric constant')*(10**9)
        
        #Get the dimensions of space
        dim = space.shape
        
        #Calculate sigma_max
        sigma_max = gf.calc_sigma_max(m, eps_0, R, d, delta_x)
        
        #Calculate how sigma evolves throughout the PML
        pml_1d_e = np.arange(1, d+1)
        pml_1d_e = gf.calc_sigma(pml_1d, d, m, sigma_max)
        pml_vert_e = pml_1d_e
        
        #Create the 1D array for the magnetic conductivity
        pml_1d_m = gf.cE2cM(pml_1D_e, eps_0)
        pml_vert_m = pml_1d_m
        
        #Create the 2D array of eps_0
        pml_eps = np.ones((dim[0], d))*eps_0
        
        #Convert these 1D arrays into 2D arrays
        for i in range(dim[0]-1):
            pml_vert_e = np.vstack(pml_vert_e, pml_1d_e)
            pml_vert_m = np.vstack(pml_vert_m, pml_1d_m)
        
        #Append the three matrices together
        pml_vert = np.dstack(pml_eps, pml_vert_e, pml_vert_m)
        
        #Append this matrix to either side of the space
        pml = np.hstack(pml_vert, space, pml_vert)
        
        #Make the horizontal PML's
        
        #Append zero matrices to the vertical PML's
        cap = np.zeros((d,d))
        pml_horz_e = np.hstack(cap, pml_vert_e, cap)
        pml_horz_m = np.hstack(cap, pml_vert_m, cap)
        horz_eps_0 = np.hstack(cap, eps_0, cap)
        
        #You still need to transpose the 3 horziontal PML's
        pml_horz_e = np.transpose(pml_horz_e)
        pml_horz_m = np.transpose(pml_horz_m)
        horz_eps_0 = np.transpose(horz_eps_0)
        
        #Append the 3 horizontal PML together
        pml_horz = np.dstack(horz_eps_0, pml_horz_e, pml_horz_m)
        
        #Append the horizontal PML to the space
        space = np.vstack(pml_horz, space, pml_horz)
        
        return space