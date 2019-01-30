import numpy as np
from general import general_functions as gf

class sources:
    #Create function that models a Gaussian Beam along a 1D array and propogates it through time
    def gaussian_src(E_0, l, omega, w_0, focal_origin, alpha, x, y_min, y_max, delta_y, t_0, t_f, delta_t):
        #Declare some matrices/variables
        y = np.arange(y_min, y_max+delta_y, delta_y)
        T = np.arange(t_0, t_f+delta_t, delta_t/2)
        E = np.zeros((y.shape[0], T.shape[0]/2))
        H = E
        t = t_0
        zr = np.pi*(w_0**2)/l
        coords = gf.convert_coord(x, y_min, y_max, focal_origin, alpha, delta_y)
        w = gf.w(coords[0, :], w_0, zr, l)
        phi_t = gf.phi_t(coords[0, :], zr, l)
        amp = gf.amp(E_0, w_0, w, coords[1, :])

        #Go through time and evaluate the beam along the column
        for i in range(T.shape[0]):
            #Calculate the E field and store it
            E_t = gf.amp2E(amp, phi_t, omega, t)
            E[:,i] = E_t

            #Half increment time
            t += delta_t/2

            #Calculate the H field (377 = impedance of free space)
            H_t = gf.amp2E(amp, phi_t, omega, t)/377
            H[:,i] = H_t

            #Half increment time
            t += delta_t/2

        return [E, H]

    def gaussian_src_1_time(E_0, l, omega, w_0, focal_origin, alpha, x, y_min, y_max, delta_y, t):
        # Declare some matrices/variables
        y = np.arange(y_min, y_max + delta_y, delta_y)
        zr = np.pi * (w_0 ** 2) / l
        coords = gf.convert_coord(x, y_min, y_max, focal_origin, alpha, delta_y)
        w = gf.w(coords[0, :], w_0, zr, l)
        phi_t = gf.phi_t(coords[0, :], zr, l)
        amp = gf.amp(E_0, w_0, w, coords[1, :])

        # Calculate the E field and store it
        E = gf.amp2E(amp, phi_t, omega, t)

        # Calculate the H field (377 = impedance of free space)
        H = E / 377

        return [E, H]