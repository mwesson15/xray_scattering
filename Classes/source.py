import numpy as np
from general import general_functions as gf

class sources:
    #Create function that models a Gaussian Beam along a 1D array and propogates it through time
    def gaussian_src(E_0, l, omega, w_0, focal_origin, alpha, x, y_min, y_max, delta_y, t_0, t_f, delta_t):
        #Declare some matrices/variables
        y_steps = ((y_max-y_min)/delta_y) + 1
        y_steps = int(y_steps)
        n_steps = ((t_f-t_0)/delta_t)+1
        n_steps = int(n_steps)
        E = np.zeros((y_steps, n_steps))
        H = E
        t = t_0
        zr = np.pi*(w_0**2)/l
        coords = gf.convert_coord(x, y_min, y_max+delta_y, focal_origin, alpha, delta_y)
        w = gf.w(coords[0, :], w_0, zr, l)
        phi_t = gf.phi_t(coords[0, :], zr, l)
        phi_l = gf.phi_l(coords[1,:], coords[0,:], zr, l)
        phase = phi_l + phi_t
        amp = gf.amp(E_0, w_0, w, coords[1, :])

        #Go through time and evaluate the beam along the column
        for i in range(n_steps):
            #Calculate the E field and store it
            E_t = gf.amp2E(amp, phase, omega, t)
            E[:,i] = E_t

            #Half increment time
            t += delta_t/2

            #Calculate the H field (377 = impedance of free space)
            H_t = gf.amp2E(amp, phase, omega, t)/377
            H[:,i] = H_t

            #Half increment time
            t += delta_t/2

        return [E, H]

    def gaussian_src_1_time(E_0, l, omega, w_0, focal_origin, alpha, x, y_min, y_max, delta_y, t):
        # Declare some matrices/variables
        zr = np.pi * (w_0 ** 2) / l
        coords = gf.convert_coord(x, y_min, y_max+delta_y, focal_origin, alpha, delta_y)
        w = gf.w(coords[0, :], w_0, zr, l)
        phi_t = gf.phi_t(coords[0, :], zr, l)
        amp = gf.amp(E_0, w_0, w, coords[1, :])

        # Calculate the E field and store it
        E = gf.amp2E(amp, phi_t, omega, t)

        # Calculate the H field (377 = impedance of free space)
        H = E / 377

        return [E, H]