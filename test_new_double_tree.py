import new_double_tree
import numpy as np 


npts1, npts2 = 1e6, 1e6
Lbox = 250.

x1in = np.random.uniform(0, Lbox, npts1)
y1in = np.random.uniform(0, Lbox, npts1)
z1in = np.random.uniform(0, Lbox, npts1)

x2in = np.random.uniform(0, Lbox, npts2)
y2in = np.random.uniform(0, Lbox, npts2)
z2in = np.random.uniform(0, Lbox, npts2)

rmax = 20.
approx_cell_size = 15.

result = new_double_tree.sorted_points_and_ids(
    x1in, y1in, z1in, x2in, y2in, z2in, 
    approx_cell_size, approx_cell_size, approx_cell_size, 
    approx_cell_size, approx_cell_size, approx_cell_size, 
    rmax, rmax, rmax, Lbox, Lbox, Lbox)
