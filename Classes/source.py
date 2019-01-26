import numpy as np


class gaussian_beam:
    # Define function that takes a constant x value and y_min, y_max and converts it to (r,z) coord
    # Origin is the (x,y) coordinates of r=z=0
    def convert_coord_column(x, y_min, y_max, origin, delta_y):
        y = np.arange(y_min, y_max, delta_y)

        delta_x = origin[0] - x
        delta_y = np.subtract(origin[1], y)

        alpha = np.arctan(np.abs(origin[1] / origin[0]))

        z = np.sqrt(np.add((delta_x * np.cos(alpha)) ** 2,
                           np.square((np.multiply(delta_y, np.sin(alpha))))))
        r = np.sqrt(np.add(np.square(np.multiply(delta_y, np.cos(alpha))),
                           ((delta_x * np.sin(alpha)) ** 2)))

        return np.vstack((r, z))

    # Define function that takes coordinates matrix and returns spatial factors
    def spatial_factors_column(E_0, coord, l, w_0):
        m = np.add(np.square(np.divide(np.multiply(coord[:, 1], l), w_0 * np.pi)), 1)
        A = np.multiply(np.divide(
            E_0, np.sqrt(m)), np.exp(np.multiply(
            np.square(np.divide(coord[:, 0], np.multiply(w_0, np.sqrt(m)))), -1)))
        phase = np.subtract(np.add(np.divide(
            np.multiply(2 * np.pi * l, np.square(
                np.divide(coord[:, 0], coord[:, 1]))), m), np.multiply(
            2 * np.pi / l, coord[:, 1])), np.arctan(np.multiply(l / (np.pi * w_0), coord[:, 1])))
        B = np.multiply(A, np.cos(phase))
        C = np.multiply(A, np.sin(phase))
        return np.vstack((B, C))

    # Define function that converts matrix of spatial factors into E values
    def spatial_E_column(spatial, omega, t):
        return np.add(np.multiply(
            spatial[:, 0], np.cos(omega * t)), np.multiply(spatial[:, 1], np.sin(omega * t)))

    # Define function that inputs the spatial factors and temporal info to output E
    def spatial_E(spatial, omega, t):
        return (spatial[0] * np.cos(omega * t)) + (spatial[1] * np.sin(omega * t))