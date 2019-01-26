import numpy as np
from scipy import constants as cst
from general_functions.py import general_functions as gf


class build_sample:
    # This class is currently on used to build a 2D environment

    # Define function that creates a block with a vacuum on top
    # Give the dim as an array of (x, y)
    # Give the dimensions of the entire space
    def create_block(dim, vacuum_height, delta_x, delta_t, relative_epsilon,
                     relative_mu, conductivity, conductivity_m):
        # Convert relative epsilon to epsilon
        eps = relative_epsilon * cst.value('electric constant')

        # Convert relative mu to mu
        mu = relative_mu * cst.value('mag. constant')

        # Calculate the two constants needed when updating the E-field in a vacuum
        const_E_vac = gf.calc_const_E(cst.value('electric constant'), 0, delta_t)

        # Calculate the two constants needed when updating the E-field in the block
        const_E_block = gf.calc_const_E(eps, conductivity, delta_t)

        # Calculate the two constants needed when updating the H-field in a vacuum
        const_H_vac = gf.calc_const_H(cst.value('mag. constant'), 0, delta_t)

        # Calculate the two constants needed when updating the H-field in the block
        const_H_block = gf.calc_const_H(mu, conductivity_m, delta_t)

        # Condense the vacuum constants into 1 matrix
        const_vac = np.concatenate((const_E_vac, const_H_vac))

        # Condense the block constants into 1 matrix
        const_block = np.concatenate((const_E_block, const_H_block))

        # Create the space for the vacuum
        # I'm assuming there's a better way of doing this step, but it isn't a big deal
        space_vac = np.ones((vacuum_height / delta_x, dim[0] / delta_x, 4))
        for i in range(4):
            space_vac[:, :, i] = np.multiply(space_vac[:, :, i], const_vac[i])

        # Create the space for the block
        space_block = np.ones(((dim[1] - vacuum_height) / delta_x, dim[0] / delta_x, 4))
        for i in range(4):
            space_block[:, :, i] = np.multiply(space_block[:, :, i], const_block[i])

        # Append the two spaces to create the total space
        space = np.append(space_vac, space_block, axis=0)

        return space

    # Define function that allows you to add a substrate
    # Currently just a rectangle
    # Give the center in terms of x
    def add_substrate(space, center, block_height, rel_eps, rel_mu,
                      conduct, conduct_m, delta_x, delta_t):
        # Convert relative epsilon to epsilon
        eps = rel_eps * cst.value('electric constant')

        # Convert relative mu to mu
        mu = rel_mu * cst.value('mag. constant')

        # Calculate the two constants needed when updating the E-field
        const_E = gf.calc_const_E(eps, cond, delta_t)

        # Calculate the two constants needed when updating the H-field
        const_H = gf.calc_const_H(mu, cond_m, delta_t)

        # Condense constants into 1 array
        const = np.concatenate((const_E, const_H))

        # Get the dimensions of space
        dim = space.shape

        # Convert given dimensions from length to index
        ind_dim = dim / delta_x

        # Calculate height of block in terms of index
        block_ind = block_height / delta_x

        # Create the substrate space
        substrate = np.ones((ind_dim[1], ind_dim[0], 4))
        for i in range(4):
            substrate[:, :, i] = np.mutiply(substrate[:, :, i], const[i])

        # Convert the given center to an index
        # ind_center = coord2index(center, total_dim, delta_x)
        ind_center = center

        # Determine the range of the added submatrix
        row_low = total_dim[0] - 1 - block_ind
        row_hi = row_low - ind_dim[0]

        column_low = ind_center - (ind_dim[1] / 2)
        column_hi = colum_low + ind_dim[1]

        # Redefine space
        space[row_hi:row_low, column_low:column_hi] = substrate

        return space

    # Unsure of the best way to do the process below, but we found the book on how to
    # do it. So I'll get back to it after I have time to read up on it

    # Create functions that add buffers for the scattered field regions
    # Give width in terms of index
    def add_scatter_reg(space, width, height):
        # Get current dimensions of space
        dim = space.shape

        # Generate the three individual matrices for the left
        eps_side = np.ones((dim[0], width_left)) * cst.value('electric constant') * (10 ** 9)
        cond_side = np.zeros((dim[0], width_left))
        # No need for 3rd matrix as it is the same as the conductivity (0)

        # Combine them
        side = np.concatenate((eps_side, cond_side,
                               cond_side), axis=2)

        # Append them to the original space
        space = np.concatenate((side, space, side), axis=0)

        # Get the new dimensions of the space
        dim = space.shape

        # Build the top and bottom regions
        eps_vert = np.ones((height, dim[1])) * cst.value('electric constant') * (10 ** 9)
        cond_vert = np.zeros((height, dim[1]))

        # Combine them
        vert = np.dstack(eps_vert, cond_vert, cond_vert)

        # Add them to the space
        space = np.vstack(vert, space, vert)

        return space