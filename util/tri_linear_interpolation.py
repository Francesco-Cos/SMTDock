from z3 import *

def find_lattice_neighbors_1d(x, y, z, enc, res, grid_size=41, grid_spacing=0.375):
    # Compute integer lattice indices, shifting for a centered grid
    if enc == 'z3':
        L = (grid_size * res) / 2  # Half the total length, to shift coordinates

        i = ToInt((x + L) / res)
        j = ToInt((y + L) / res)
        k = ToInt((z + L) / res)
    else:
        L = (grid_size * grid_spacing) / 2  # Half the total length, to shift coordinates

        i = int((x + L) // grid_spacing)
        j = int((y + L) // grid_spacing)
        k = int((z + L) // grid_spacing)
    
    # Compute 1D indices for the 8 neighboring points
    def index(ix, iy, iz):
        return ix + grid_size * (iy + grid_size * iz)

    neighbors = [
        index(i, j, k), index(i+1, j, k), index(i, j+1, k), index(i, j, k+1),
        index(i+1, j+1, k), index(i+1, j, k+1), index(i, j+1, k+1), index(i+1, j+1, k+1)
    ]
    return neighbors

def trilinear_interpolation(values, neighbors, enc, grid_spacing=0.375):
    # Extract function values at the 8 neighboring points
    if enc == 'z3':
        c000 = Select(values, neighbors[0])
        # c100 = Select(values, neighbors[1])
        # c010 = Select(values, neighbors[2])
        # c001 = Select(values, neighbors[3])
        # c110 = Select(values, neighbors[4])
        # c101 = Select(values, neighbors[5])
        # c011 = Select(values, neighbors[6])
        # c111 = Select(values, neighbors[7])
    else:
        c000 = values[neighbors[0]]
        c100 = values[neighbors[1]]
        c010 = values[neighbors[2]]
        c001 = values[neighbors[3]]
        c110 = values[neighbors[4]]
        c101 = values[neighbors[5]]
        c011 = values[neighbors[6]]
        c111 = values[neighbors[7]]

        # Perform trilinear interpolation
        c00 = c000 * c100 
        c01 = c001 * c101
        c10 = c010 * c110
        c11 = c011 * c111

        c0 = c00 * c10
        c1 = c01 * c11

        c = c0 * c1

    return c000 if enc=='z3' else c000 


def apply_interpolation(x, y, z, values, enc, res):
    n = find_lattice_neighbors_1d(x, y, z, enc, res)
    return trilinear_interpolation(values, n, enc)
