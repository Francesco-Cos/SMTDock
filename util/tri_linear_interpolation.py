from z3 import *

def find_lattice_neighbors_1d(x, y, z, enc, L, index_factor, grid_size):
    # Compute integer lattice indices, shifting for a centered grid
    if enc == 'z3':
        i = ToInt((x + L) * index_factor)
        j = ToInt((y + L) * index_factor)
        k = ToInt((z + L) * index_factor)
    else:
        i = int((x + L) * index_factor)
        j = int((y + L) * index_factor)
        k = int((z + L) * index_factor)
    
    # Compute 1D indices for the neighboring points
    def index(ix, iy, iz):
        return ix + grid_size * (iy + grid_size * iz)

    # neighbors = [
    #     index(i, j, k), index(i+1, j, k), index(i, j+1, k), index(i, j, k+1),
    #     index(i+1, j+1, k), index(i+1, j, k+1), index(i, j+1, k+1), index(i+1, j+1, k+1)
    # ]
    return index(i, j, k)
   
def trilinear_interpolation(values, neighbor, enc):
    # Extract function values at neighboring points
    if enc == 'z3':
        c = Select(values, neighbor)
    else:
        c = values[neighbor]
    return c if enc=='z3' else c


def apply_interpolation(x, y, z, values, enc, L, index_factor, grid_size=41):
    n = find_lattice_neighbors_1d(x, y, z, enc, L, index_factor, grid_size)
    return trilinear_interpolation(values, n, enc)

