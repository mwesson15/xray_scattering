import numpy as np
from scipy import constants as cst
import gf

#These functions don't current take into account conductivities
#Define function that updates the E field
def update_E(self, prev_E, H, delta_x, delta_t, source):
    #Calculate the curl of the H field
    curl_H = gf.curl(H, delta_x)

    #Calculate the constant that the curl is multiplied by in the update eq'n
    const = delta_t/(cst.constant('electric constant')*(10**9))

    #Calculate the new E values
    E = np.add(prev_E, np.multiply(const, curl_H))

    #Go thru boundaries and implement boundary conditions

    return E

#Define function that updates the H field
def update_H(self, prev_H, E, delta_x, delta_t, source):
    #Calculate the curl of the E field
    curl_E = gf.curl(E, delta_x)

    #Get the dimension of the space
    dim = E.shape

    #Due to the way we defined where the H values are stored (a +delta_x/2 step in the
    # x and y direction) we actually calculated the curl at the (i+1, j+1) index
    #We will now realign this by moving the bottom row to the top then the last column
    #to the front
    fix_row = np.vsplit(curl_E, [dim[0]-1])
    curl_E = np.append(fix_row[1], fix_row[0], axis = 0)

    fix_column = np.hsplit(curl_E, [dim[1]-1])
    curl_E = np.append(fix_column[1], fix_column[0], axis = 1)


    #Calculate the constant that the curl is multiplied by in the update eq'n
    const = delta_t/(cst.constants('mag. constant')*(10**-9))

    #Calculate the new H values
    H = np.add(prev_H, np.multiply(const, curl_E))

    #Go thru boundaries and implement boundary conditions

    return H
