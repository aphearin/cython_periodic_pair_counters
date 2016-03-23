import numpy as np
from math import ceil, floor

def sorted_points_and_ids(xin, yin, zin, xperiod, yperiod, zperiod, 
    approx_xcell_size, approx_ycell_size, approx_zcell_size):
    """ Determine the cell_id of every point, sort the points 
    according to cell_id, and return the sorted points as well as 
    the cell id indexing array. 

    Notes 
    -----
    The x-coordinates of points with cell_id = icell are given by 
    xout[cell_id_indices[icell]:cell_id_indices[icell+1]]. 
    """
    npts = len(xin)
    num_xdivs, xcell_size = determine_cell_size(xperiod, approx_xcell_size)
    num_ydivs, ycell_size = determine_cell_size(yperiod, approx_ycell_size)
    num_zdivs, zcell_size = determine_cell_size(zperiod, approx_zcell_size)
    ncells = num_xdivs*num_ydivs*num_zdivs

    ix = digitized_position(xin, xcell_size, num_xdivs)
    iy = digitized_position(yin, ycell_size, num_ydivs)
    iz = digitized_position(zin, zcell_size, num_zdivs)

    cell_ids = cell_id_from_cell_tuple(ix, iy, iz, num_ydivs, num_zdivs)
    cell_id_sorting_indices = np.argsort(cell_ids)

    cell_id_indices = np.searchsorted(cell_ids, np.arange(ncells), 
        sorter = cell_id_sorting_indices)
    cell_id_indices = np.append(cell_id_indices, npts)

    xout = np.ascontiguousarray(xin[cell_id_sorting_indices], dtype=np.float64)
    yout = np.ascontiguousarray(yin[cell_id_sorting_indices], dtype=np.float64)
    zout = np.ascontiguousarray(zin[cell_id_sorting_indices], dtype=np.float64)

    cell_id_indices = np.ascontiguousarray(cell_id_indices, dtype=np.int64)

    return xout, yout, zout, cell_id_indices

def determine_cell_size(period, approx_cell_size):
    """ Find closest cell size that evenly divides period.
    """
    num_divs = int(round(period/float(approx_cell_size)))
    cell_size = period/float(num_divs)
    return num_divs, cell_size

def cell_id_from_cell_tuple(ix, iy, iz, num_ydivs, num_zdivs):
    """ (ix, iy, iz) ==> icell
    """
    return ix*(num_ydivs*num_zdivs) + iy*num_zdivs + iz

def cell_tuple_from_cell_id(cell_id, num_ydivs, num_zdivs):
    """ icell ==> (ix, iy, iz)
    """
    nxny = num_ydivs*num_zdivs
    ix = cell_id / nxny
    iy = (cell_id - ix*nxny) / num_zdivs
    iz = cell_id - (ix*num_ydivs*num_zdivs) - (iy*num_zdivs)
    return ix, iy, iz

def digitized_position(p, cell_size, num_divs):
    ip = np.floor(p/cell_size).astype(int)
    return np.where(ip >= num_divs, num_divs-1, ip)

def set_sample1_cell_sizes(period, search_length, approx_cell_size, 
    max_cells_per_dimension = 25):
    if search_length > period/3.:
        msg = ("Input ``search_length`` cannot exceed period/3")
        raise ValueError(msg)

    ndivs = int(floor(period/float(approx_cell_size)))
    ndivs = max(ndivs, 1)
    ndivs = min(max_cells_per_dimension, ndivs)

    nsearch = int(floor(period/float(search_length)))
    nsearch = max(nsearch, 1)

    ndivs = min(ndivs, nsearch)
    ndivs = max(3, ndivs)
    cell_size = period/float(ndivs)

    return cell_size

def set_sample2_cell_sizes(sample1_cell_size, approx_cell_size, period, 
    max_cells_per_dimension = 25):
    ndivs_sample1_cells = int(round(sample1_cell_size/float(approx_cell_size)))
    ndivs_sample1_cells = max(1, ndivs_sample1_cells)
    ndivs_sample1_cells = min(max_cells_per_dimension, ndivs_sample1_cells)
    cell_size = sample1_cell_size/ndivs_sample1_cells
    min_cell_size = period/float(max_cells_per_dimension)
    cell_size = max(cell_size, min_cell_size)
    return cell_size

def lower_bound_cell_index(ip, num_cell2_per_cell1):
    return ip*num_cell2_per_cell1

def upper_bound_cell_index(ip, num_cell2_per_cell1):
    return (ip+1)*num_cell2_per_cell1 - 1

def num_cells_to_cover_search_length(cell_size, search_length):
    return int(np.ceil(search_length/float(cell_size)))




