import numpy as np
cimport numpy as cnp
import tree_module 

def cython_npairs(double[:,:] data1, double[:, :] data2, 
	double[:] rbins, double[:] period):
	
	cdef double xperiod = period[0]
	cdef double yperiod = period[1]
	cdef double zperiod = period[2]
	cdef int num_rbins = len(rbins)
	cdef double rmax = rbins[num_rbins[num_bins-1]]


