import numpy as np
import matplotlib.pyplot as plt
from source import sources as src

#Define some variables
delta_y = 1

#Define wave properties
origin = [0, 0]
alpha = 0
l = 30
w_0 = 40
E_0 = 10
omega = 2*np.pi

t_i = 0
t_f = 1
delta_t = 1

x_min = -400
x_max = 400
x = np.arange(x_min, x_max+delta_y, delta_y)
y_min = -200
y_max = 200
y = np.arange(y_min, y_max+delta_y, delta_y)


I = np.zeros((int(y.shape[0]), int(x.shape[0])))

for i in range(int(x.shape[0])):
    fields = src.gaussian_src(E_0, l, omega, w_0, origin, alpha, x[i], y_min, y_max, delta_y, t_i, t_f, delta_t)
    E = fields[0]
    E = E[:,-1]
    H = fields[1]
    H = H[:, -1]
    I_x = np.multiply(E, H)

    I[:,i] = I_x

#Plot I
plt.figure(dpi=100)
plt.imshow(I, interpolation='spline36', cmap='RdBu', alpha=0.9)
plt.axis('off')
plt.colorbar()
plt.show()