#!/usr/bin/env python
import numpy as np 
from npairs import npairs
from pure_python_pair_counts import python_npairs
from time import time 

npts1, npts2 = 1e5, 1e5
Lbox = 1.
rbins = np.linspace(0.001, 0.05, 10)

x1 = np.random.uniform(0, Lbox, npts1)
y1 = np.random.uniform(0, Lbox, npts1)
z1 = np.random.uniform(0, Lbox, npts1)
x2 = np.random.uniform(0, Lbox, npts2)
y2 = np.random.uniform(0, Lbox, npts2)
z2 = np.random.uniform(0, Lbox, npts2)

start = time()
result = npairs(x1, y1, z1, x2, y2, z2, rbins, Lbox, 
	approx_cell_size = 0.1)

runtime = time() - start
print("Total runtime = %f seconds" % runtime)




# pure_python_result = python_npairs(x1, y1, z1, x2, y2, z2, rbins, Lbox)
# assert np.all(pure_python_result == result)
# print(pure_python_result)
# print(result)
