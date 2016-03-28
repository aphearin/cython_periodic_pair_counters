import numpy as np
cimport numpy as cnp
import new_double_tree 

def cython_npairs(double[:,:] data1, double[:, :] data2, 
    double[:] rbins, double[:] period):
    
    cdef double xperiod = period[0]
    cdef double yperiod = period[1]
    cdef double zperiod = period[2]
    cdef int num_rbins = len(rbins)
    cdef double rmax = rbins[num_rbins[num_bins-1]]

    result = new_double_tree.sorted_points_and_ids(data1[:,0], data1[:,1], data1[:,2], 
        data2[:,0], data2[:,1], data2[:,2], 
        rmax, rmax, rmax, rmax, rmax, rmax, rmax, rmax, rmax, 
        period[0], period[1], period[2])
    x1, y1, z1, cell1_id_indices, x1cell_size, y1cell_size, z1cell_size = result[0:7]
    x2, y2, z2, cell2_id_indices, x2cell_size, y2cell_size, z2cell_size = result[7:]

    cdef int num_cell1 = len(cell1_id_indices)
    cdef int num_cell2 = len(cell2_id_indices)
    cdef int ix1, iy1, iz1, icell1, icell2 

    cdef long[:] counts = np.zeros(len(rbins), dtype=long)

    for icell1 in range(num_cell1):

        x_icell1 = x1[cell1_id_indices[icell1]:cell1_id_indices[icell1+1]]
        y_icell1 = y1[cell1_id_indices[icell1]:cell1_id_indices[icell1+1]]
        z_icell1 = z1[cell1_id_indices[icell1]:cell1_id_indices[icell1+1]]

        ix1, iy1, iz1 = new_double_tree.cell_tuple_from_cell_id(icell1)

        num_x2_covering_steps = (
            new_double_tree.num_cells_to_cover_search_length(x1cell_size, rmax))
        num_y2_covering_steps = (
            new_double_tree.num_cells_to_cover_search_length(y1cell_size, rmax))
        num_z2_covering_steps = (
            new_double_tree.num_cells_to_cover_search_length(z1cell_size, rmax))

        x_icell2 = x2[cell2_id_indices[icell2]:cell2_id_indices[icell2+1]]
        y_icell2 = y2[cell2_id_indices[icell2]:cell2_id_indices[icell2+1]]
        z_icell2 = z2[cell2_id_indices[icell2]:cell2_id_indices[icell2+1]]





