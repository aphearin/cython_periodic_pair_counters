cimport cython
import numpy as np
cimport numpy as np

__all__ = ('npairs_no_pbc', )

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def npairs_no_pbc(np.ndarray[np.float64_t, ndim=1] x_icell1,
                  np.ndarray[np.float64_t, ndim=1] y_icell1,
                  np.ndarray[np.float64_t, ndim=1] z_icell1,
                  np.ndarray[np.float64_t, ndim=1] x_icell2,
                  np.ndarray[np.float64_t, ndim=1] y_icell2,
                  np.ndarray[np.float64_t, ndim=1] z_icell2,
                  np.ndarray[np.float64_t, ndim=1] rbins):
    
    #c definitions
    cdef int nbins = len(rbins)
    cdef np.ndarray[np.int_t, ndim=1] counts = np.zeros((nbins,), dtype=np.int)
    cdef double dsq
    cdef int i, j, k, last_bin_to_increment
    cdef int Ni = len(x_icell1)
    cdef int Nj = len(x_icell2)
    cdef double dx, dy, dz
    
    #square the distance bins to avoid taking a square root in a tight loop
    rbins = rbins*rbins
    #loop over points in grid1's cells
    for i in range(0,Ni):
        #loop over points in grid2's cells
        for j in range(0,Nj):

            #calculate the square distance
            dx = x_icell1[i] - x_icell2[j]
            dy = y_icell1[i] - y_icell2[j]
            dz = z_icell1[i] - z_icell2[j]
            dsq = dx*dx + dy*dy + dz*dz

            k = nbins-1
            while dsq <= rbins[k]:
                counts[k] += 1
                k=k-1
                if k<0: break
        
    return counts

