import numpy as np
import multiprocessing
from functools import partial

from double_tree import FlatRectanguloidDoubleTree
from cpairs import npairs_no_pbc

__all__ = ('npairs', )

def npairs(x1, y1, z1, x2, y2, z2, rbins, Lbox, 
    num_threads = 1, approx_cell_size = None):

    xperiod, yperiod, zperiod = Lbox, Lbox, Lbox 
    rmax = np.max(rbins)
        
    ### Compute the estimates for the cell sizes
    if approx_cell_size is None:
        approx_x1cell_size, approx_y1cell_size, approx_z1cell_size = (
            rmax, rmax, rmax
            )
        approx_x2cell_size, approx_y2cell_size, approx_z2cell_size = (
            rmax, rmax, rmax
            )
    else:
        approx_x1cell_size = approx_cell_size
        approx_y1cell_size = approx_cell_size
        approx_z1cell_size = approx_cell_size
        approx_x2cell_size = approx_cell_size
        approx_y2cell_size = approx_cell_size
        approx_z2cell_size = approx_cell_size

    double_tree = FlatRectanguloidDoubleTree(
        x1, y1, z1, x2, y2, z2,  
        approx_x1cell_size, approx_y1cell_size, approx_z1cell_size, 
        approx_x2cell_size, approx_y2cell_size, approx_z2cell_size, 
        rmax, rmax, rmax, xperiod, yperiod, zperiod, PBCs=True)
        
    #number of cells
    Ncell1 = double_tree.num_x1divs*double_tree.num_y1divs*double_tree.num_z1divs
    Ncell2 = double_tree.num_x2divs*double_tree.num_y2divs*double_tree.num_z2divs
        
    #create a function to call with only one argument
    engine = partial(_npairs_engine, double_tree, rbins)
    
    #do the pair counting
    if num_threads > 1:
        pool = multiprocessing.Pool(num_threads)
        cell1_chunk_list = np.array_split(xrange(Ncell1), num_threads)
        result = pool.map(engine,cell1_chunk_list)
        pool.close()
        counts = np.sum(np.array(result), axis=0)
    if num_threads == 1:
        counts = engine(xrange(Ncell1))
    
    return counts


def _npairs_engine(double_tree, rbins, cell1_list):
    """
    pair counting engine for npairs function.  
    This code calls a cython function.
    """    
    counts = np.zeros(len(rbins))
    
    for icell1 in cell1_list:

        #extract the points in the cell
        s1 = double_tree.tree1.slice_array[icell1]
        x_icell1 = double_tree.tree1.x[s1]
        y_icell1 = double_tree.tree1.y[s1]
        z_icell1 = double_tree.tree1.z[s1]
            
        xsearch_length = rbins[-1]
        ysearch_length = rbins[-1]
        zsearch_length = rbins[-1]

        adj_cell_generator = double_tree.adjacent_cell_generator(
            icell1, xsearch_length, ysearch_length, zsearch_length)
                
        for icell2, xshift, yshift, zshift in adj_cell_generator:
            #extract the points in the cell
            s2 = double_tree.tree2.slice_array[icell2]
            x_icell2 = double_tree.tree2.x[s2] + xshift
            y_icell2 = double_tree.tree2.y[s2] + yshift 
            z_icell2 = double_tree.tree2.z[s2] + zshift


            #use cython functions to do pair counting
            counts += npairs_no_pbc(
                x_icell1, y_icell1, z_icell1,
                x_icell2, y_icell2, z_icell2,
                rbins)

            
    return counts


