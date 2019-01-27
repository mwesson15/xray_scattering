import numpy as np
from scipy import constants as cst

class general_functions:
    #Define function that normalizes vector
    def normalize(vector):
        #Calculate magnitude of vector
        amp = (np.dot(vector, vector))**0.5

        return np.divide(vector, amp)

    #Define function that converts E to E-tilda which has a similar magnitude as H
    #Some resources use this reduction while others don't
    def reduce_E(E):
        #Calculate factor to multiply it by
        factor = (cst.value('electric constant')/cst.value('mag. constant'))**0.5

        return np.multiply(factor, E)

    #Define function that converts from index to x, y, or z coordinate
    def index2coord(index, space, delta_x):
        #Get the dimensions of space
        dim = space.shape

        #Calculate index of the origin
        origin_index = np.divide(np.subtract(dim, 1), 2)

        #Calculate distance between given index and origin
        index_spacing = np.subtract(index, origin_index)

        #Convert index_spacing into nanometers
        coord = np.multiply(delta_x, index_spacing)

        return coord

    #Define function that converts from (x,y,z) coordinates to index
    def coord2index(coordinates, space, delta_x):
        #Get the dimensions of space
        dim = space.shape

        #Convert from nm to index distance
        index_spacing = np.divide(coordinates, delta_x)

        #Calculate index of the origin
        origin_index = np.divide(np.subtract(dim, 1), 2)

        #Calculate the index
        index = np.add(origin_index, index_spacing)

        return index

    #We'll now define functions needed for the update equations

    #Define function that calculates the curl of the H field at all points
    def curl_H(H, delta_x):
        #Assumes and H values is stored in the (+x,+y) compared to the E values with the same index
        #Need to calculate {dHz_dx, dHz_dy, dHx_dy, dHy_dx}

        #Shift the z values in the x direction then calculate the x component of the curl
        xshifted_z = np.roll(H[:,:,2], 1, axis = 1)
        curl_x = np.divide(np.subtract(xshifted_z, H[:,:,2]), delta_x)

        #Shift the z values in the y direction then calculate the y component of the curl
        yshifted_z = np.roll(H[:,:,2], 1, axis = 0)
        curl_y = np.divide(np.subtract(yshifted_z, H[:,:,2]), delta_x)

        #Calculate the z component of the curl
        yshifted_x = np.roll(H[:,:,0], 1, axis = 0)
        xshifted_y = np.roll(H[:,:,1], 1, axis = 1)
        dHy_dx = np.divide(np.subtract(xshifted_y, H[:,:,1]), delta_x)
        dHx_dy = np.divide(np.subtract(yshifted_x, H[:,:,0]), delta_x)
        curl_z = np.subtract(dHy_dx, dHx_dy)

        return [curl_x, curl_y, curl_z]

    #Define function that calculates the curl of the E field at all points
    def curl_E(E, delta_x):
        #Assumes and H values is stored in the (+x,+y) compared to the E values with the same index

        #Shift the z values in the x direction then calculate the x component of the curl
        xshifted_z = np.roll(E[:,:,2], -1, axis = 1)
        curl_x = np.divide(np.subtract(xshifted_z, E[:,:,2]), delta_x)

        #Shift the z values in the y direction then calculate the y component of the curl
        yshifted_z = np.roll(E[:,:,2], 1, axis = 0)
        curl_y = np.divide(np.subtract(yshifted_z, E[:,:,2]), delta_x)

        #Calculate the z component of the curl
        yshifted_x = np.roll(E[:,:,0], 1, axis = 0)
        xshifted_y = np.roll(E[:,:,1], 1, axis = 1)
        dEy_dx = np.divide(np.subtract(xshifted_y, E[:,:,1]), delta_x)
        dEx_dy = np.divide(np.subtract(yshifted_x, E[:,:,0]), delta_x)
        curl_z = np.subtract(dEy_dx, dEx_dy)

        return [curl_x, curl_y, curl_z]

    #Calculate the curl of the H field at 1 point
    def curl_point_H(H, index, delta_x):
        #Calculate x component
        curl_x = (H[index[0]+1, index[1], 2] - H[index[0], index[1], 2])/delta_x

        #Calculate the y component
        curl_y = (H[index[0], index[1], 2] - H[index[0], index[1]+1, 2])/delta_x

        #Calculate the z component by first calculating the 2 derivatives needed to calculate it
        dHy_dx = (H[index[0], index[1]+1, 1] - H[index[0], index[1], 1])/delta_x
        dHx_dy = (H[index[0]+1, index[1], 0] - H[index[0], index[1], 0])/delta_x
        curl_z = dHy_dx - dHx_dy

        return [curl_x, curl_y, curl_z]

    #Calculate the curl of the E field at 1 point
    def curl_point_E(E, index,delta_x):
        #Calculate x component
        curl_x = (E[index[0], index[1], 2] - E[index[0]-1, index[1], 2])/delta_x

        #Calculate the y component
        curl_y = (E[index[0], index[1]-1, 2] - E[index[0], index[1], 2])/delta_x

        #Calculate the z component by first calculating the 2 derivatives needed to calculate it
        dEy_dx = (E[index[0], index[1], 1] - E[index[0], index[1]-1, 1])/delta_x
        dEx_dy = (E[index[0], index[1], 0] - E[index[0]-1, index[1], 0])/delta_x
        curl_z = dEy_dx - dEx_dy

        return [curl_x, curl_y, curl_z]

    #Define function that calculates the two constants for updating the E values
    def calc_const_E(eps, cond, delta_t):
        #Calculate the one for the previous E value
        numerator_1 = np.subtract(np.multiply(2, eps), np.multiply(cond, delta_t))
        denominator = np.add(np.multiply(2, eps), np.multiply(cond, delta_t))
        c1 = np.divide(numerator_1, denominator)

        #Calculate the one for the curl component
        numerator_2 = np.multiply(2, delta_t)
        c2 = np.divide(numerator_2, denominator)

        return [c1, c2]

    #Define functions that calculates the two constants for updating the H values
    def calc_const_H(mu, cond_m, delta_t):
        #Calculate the one for the previous E value
        numerator_1 = np.subtract(np.multiply(2, mu), np.multiply(cond_m, delta_t))
        denominator = np.add(np.multiply(2, mu), np.multiply(cond_m, delta_t))
        c1 = np.divide(numerator_1, denominator)

        #Calculate the one for the curl component
        numerator_2 = np.multiply(2, delta_t)
        c2 = np.divide(numerator_2, denominator)

        return [c1, c2]

    #Define function that calculates only the phase of the wave at a certain point
    def phase(k, r, r_0, omega, t):
        #Calculate distance between r and r_0
        c = r - r_0 #Gives difference between points as a vector
        R = (c[0]**2 + c[1]**2)**0.5

        #Calculate the phase
        phi = k*R - omega*t

        return phi

    #Define function that calculates k from the wavelength (nm)
    #The equation allows for this one function to convert between the 2 in either direction
    def lambda2k(wavelength):
        #Plug into equation
        k = 2*np.pi/wavelength

        return k

    #Define function that converts between angular frequency and frequency
    def omega2f(omega):
        #Plug into equation
        f = 2*np.pi/omega

        return f

    #Define function that converts between frequency and angular frequency
    def f2omega(f):
        #Plug into equation
        omega = 2*np.pi*f

        return omega

    #Define function that returns generic wave equation (E=E_0*sin[k(r-r_0)-omega*t])
    def wave(amplitude, k, R, omega, t):
        #Plug into equation
        wave = amplitude*sin(k*R - omega*t)

        return wave

    #We'll now define functions needed for the PML production

    #Define a function that converts between the conductivity and magnetic conductivity
    #sigma/epsilon = sigma_m/mu
    #Still assuming mu = mu_0
    def cE2cM(cond, eps):
        #Plug into eq'n
        cond_m = (cond/eps)*cst.value('mag. constant')

        return cond_m

    #Define function that goes from magnetic conductivity to conductivity
    def cM2cE(cond_m, eps):
        #Plug into eq'n
        cond = (cond_m*eps)/cst.value('mag. constant')

        return cond

    #Define function that calculates sigma_max
    def calc_sigma_max(m, eps, R, d, delta_x):
        #Calculate the thickness in nm
        d = d/delta_x

        #Calculate the impedance
        imp = (eps/cst.value('mag. constant'))**0.5

        #Calculate sigma_max
        sigma_max = (-1*(m+1)*np.log(R))/(2*imp*d)

        return sigma_max

    #Define function that is sigma(position)
    #Sigma(x) = (x/d)^m *sigma_max
    #d is the thickness of the PML and m is a choosen number normally between 3 and 4
    def calc_sigma(x, d, m, sigma_max):
        #Plug into eq'n
        sigma = ((x/d)**m)*sigma_max

        return sigma

    #Define functions that are used when constructing the Gaussian Beam

    #Define function that returns the Rayleigh Length
    def z_R(w_0, l):
        return np.pi*(w_0**2)/l

    #Define function that returns the beam width at 1 or many z values
    def w(z, zr, l):
        return w_0 * np.sqrt(1 + np.square(z/zr))

    #Define function that returns the radius of curvature along and array of z values
    def R(z, zr):
        return np.multiply(z, 1 + np.square(np.divide(zr, z)))

    #Define function that returns the longitudinal phase
    def phi_l(r, z, zr, l):
        return np.divide(np.square(r)*(np.pi/(l*z)), 1 + np.square(np.divide(zr, z)))

    #Define function that returns the transversal phase
    def phi_t(z, zr, l):
        return ((2*np.pi/l)*z) - np.arctan(z/zr)

    #Define function that returns the inverted complex beam parameter (circular)
    def q_inv(R, w):
        real = np.power(R, -1)
        imag = np.divide((-1*l/np.pi), np.square(w))
        return real + (1j *imag)

    #Define function that returns the spatial factors
    def spatial_beam(E_0, w_0, w, r, phi_t, phi_l):
        A = E_0 * np.divide(w_0, w) * np.exp(-1*np.square(np.divide(r, w))) * np.cos(phi_l+phi_t)
        B = E_0 * np.divide(w_0, w) * np.exp(-1*np.square(np.divide(r, w))) * np.sin(phi_l+phi_t)
        return np.vstack((A, B))

    #Define function that returns E from the spatial factors
    #Spatial matrix is 2D
    def spatial2E(spatial, omega, t):
        return spatial[0,:]*np.cos(omega*t) + spatial[1,:]*np.sin(omega*t)

    #Define function that returns E from the spatial factors
    #Spatial matrix is 1D
    def spatial2E_pt(selfspatial, omega, t):
        return spatial[0]*np.cos(omega*t) + spatial[1]*np.sin(omega*t)

    #Define function that uses complex parameter to return E
    def complexE(E_0, q_inv, r, l):
        u = q_inv * np.exp(-1j * q_inv * np.square(r) * np.pi / l)
        return E * u.real

    # Define function that takes a constant x value and y_min, y_max and converts it to (r,z) coord
    # Origin is the (x,y) coordinates of r=z=0
    def convert_coord(x, y_min, y_max, origin, delta_y):
        y = np.arange(y_min, y_max, delta_y)

        delta_x = origin[0] - x
        delta_y = np.subtract(origin[1], y)

        #Assuming origin[0] =/= 0
        alpha = np.arctan(np.abs(origin[1] / origin[0]))

        z = np.sqrt(np.add((delta_x * np.cos(alpha)) ** 2,
                           np.square((np.multiply(delta_y, np.sin(alpha))))))
        r = np.sqrt(np.add(np.square(np.multiply(delta_y, np.cos(alpha))),
                           ((delta_x * np.sin(alpha)) ** 2)))

        return np.vstack((r, z))
