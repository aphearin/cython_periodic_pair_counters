import numpy as np
from math import floor

def sorted_points_and_ids(x1in, y1in, z1in, x2in, y2in, z2in, 
    approx_x1cell_size, approx_y1cell_size, approx_z1cell_size, 
    approx_x2cell_size, approx_y2cell_size, approx_z2cell_size,
    xsearch_length, ysearch_length, zsearch_length, xperiod, yperiod, zperiod):
    """
    """
    x1cell_size, num_x1divs = sample1_cell_sizes(
        xperiod, xsearch_length, approx_x1cell_size)
    y1cell_size, num_y1divs = sample1_cell_sizes(
        yperiod, ysearch_length, approx_y1cell_size)
    z1cell_size, num_z1divs = sample1_cell_sizes(
        zperiod, zsearch_length, approx_z1cell_size)
    num_cell1 = num_x1divs*num_y1divs*num_z1divs

    x2cell_size, num_x2divs = sample2_cell_sizes(
        x1cell_size, approx_x2cell_size, xperiod)
    y2cell_size, num_y2divs = sample2_cell_sizes(
        y1cell_size, approx_y2cell_size, yperiod)
    z2cell_size, num_z2divs = sample2_cell_sizes(
        z1cell_size, approx_z2cell_size, zperiod)
    num_cell2 = num_x2divs*num_y2divs*num_z2divs

    ix1 = digitized_position(x1in, x1cell_size, num_x1divs)
    iy1 = digitized_position(y1in, y1cell_size, num_y1divs)
    iz1 = digitized_position(z1in, z1cell_size, num_z1divs)

    ix2 = digitized_position(x2in, x2cell_size, num_x2divs)
    iy2 = digitized_position(y2in, y2cell_size, num_y2divs)
    iz2 = digitized_position(z2in, z2cell_size, num_z2divs)

    cell1_ids = cell_id_from_cell_tuple(ix1, iy1, iz1, num_y1divs, num_z1divs)
    cell2_ids = cell_id_from_cell_tuple(ix2, iy2, iz2, num_y2divs, num_z2divs)

    cell1_id_sorting_indices = np.argsort(cell1_ids)
    cell2_id_sorting_indices = np.argsort(cell2_ids)

    npts1 = len(x1in)
    npts2 = len(x2in)

    cell1_id_indices = np.searchsorted(cell1_ids, np.arange(num_cell1), 
        sorter = cell1_id_sorting_indices)
    cell1_id_indices = np.append(cell1_id_indices, npts1)

    cell2_id_indices = np.searchsorted(cell2_ids, np.arange(num_cell2), 
        sorter = cell2_id_sorting_indices)
    cell2_id_indices = np.append(cell2_id_indices, npts2)

    x1out = np.ascontiguousarray(x1in[cell1_id_sorting_indices], dtype=np.float64)
    y1out = np.ascontiguousarray(y1in[cell1_id_sorting_indices], dtype=np.float64)
    z1out = np.ascontiguousarray(z1in[cell1_id_sorting_indices], dtype=np.float64)

    x2out = np.ascontiguousarray(x2in[cell2_id_sorting_indices], dtype=np.float64)
    y2out = np.ascontiguousarray(y2in[cell2_id_sorting_indices], dtype=np.float64)
    z2out = np.ascontiguousarray(z2in[cell2_id_sorting_indices], dtype=np.float64)

    cell1_id_indices = np.ascontiguousarray(cell1_id_indices, dtype=np.int64)
    cell2_id_indices = np.ascontiguousarray(cell2_id_indices, dtype=np.int64)

    return (x1out, y1out, z1out, cell1_id_indices, 
        x2out, y2out, z2out, cell2_id_indices)

def sample1_cell_sizes(period, search_length, approx_cell_size, 
    max_cells_per_dimension = 25):
    """
    """
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

    return cell_size, ndivs

def sample2_cell_sizes(sample1_cell_size, approx_cell_size, period, 
    max_cells_per_dimension = 25):
    ndivs_sample1_cells = int(round(sample1_cell_size/float(approx_cell_size)))
    ndivs_sample1_cells = max(1, ndivs_sample1_cells)
    ndivs_sample1_cells = min(max_cells_per_dimension, ndivs_sample1_cells)
    cell_size = sample1_cell_size/ndivs_sample1_cells
    min_cell_size = period/float(max_cells_per_dimension)
    cell_size = max(cell_size, min_cell_size)
    return cell_size, int(cell_size/period)
    
def digitized_position(p, cell_size, num_divs):
    ip = np.floor(p/cell_size).astype(int)
    return np.where(ip >= num_divs, num_divs-1, ip)

def cell_id_from_cell_tuple(ix, iy, iz, num_ydivs, num_zdivs):
    """ (ix, iy, iz) ==> icell
    """
    return ix*(num_ydivs*num_zdivs) + iy*num_zdivs + iz











