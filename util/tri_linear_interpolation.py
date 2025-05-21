from z3 import *

def find_lattice_neighbors_1d(x, y, z, enc, res, L, index_factor, grid_size):
    # Compute integer lattice indices, shifting for a centered grid
    if enc == 'z3':
        i = ToInt((x + L) * index_factor)
        j = ToInt((y + L) * index_factor)
        k = ToInt((z + L) * index_factor)
    else:
        i = int((x + L) * index_factor)
        j = int((y + L) * index_factor)
        k = int((z + L) * index_factor)
        xd = (x + L - i / index_factor) * index_factor
        yd = (y + L - j / index_factor) * index_factor
        zd = (z + L - k / index_factor) * index_factor
    
    # Compute 1D indices for the 8 neighboring points
    def index(ix, iy, iz):
        return ix + grid_size * (iy + grid_size * iz)

    neighbors = [
        index(i, j, k), index(i+1, j, k), index(i, j+1, k), index(i, j, k+1),
        index(i+1, j+1, k), index(i+1, j, k+1), index(i, j+1, k+1), index(i+1, j+1, k+1)
    ]
    # if enc=='z3':
    return index(i, j, k)
    # else:
    #     return neighbors, [xd,yd,zd]
 
def trilinear_interpolation(values, neighbors, enc):
    # Extract function values at the 8 neighboring points
    if enc == 'z3':
        c000 = Select(values, neighbors)
        # c100 = Select(values, neighbors[1])
        # c010 = Selecat(values, neighbors[2])
        # c001 = Select(values, neighbors[3])
        # c110 = Select(values, neighbors[4])
        # c101 = Select(values, neighbors[5])
        # c011 = Select(values, neighbors[6])
        # c111 = Select(values, neighbors[7])
    else:
        # xd = neighbors[1][0]
        # yd = neighbors[1][1]
        # zd = neighbors[1][2]

        c000 = values[neighbors]
        # c100 = values[neighbors[0][1]]
        # c010 = values[neighbors[0][2]]
        # c001 = values[neighbors[0][3]]
        # c110 = values[neighbors[0][4]]
        # c101 = values[neighbors[0][5]]
        # c011 = values[neighbors[0][6]]
        # c111 = values[neighbors[0][7]]
        
        # # Perform trilinear interpolation
        # c00 = c000 * (1 - xd) + c100 * xd
        # c01 = c001 * (1 - xd) + c101 * xd
        # c10 = c010 * (1 - xd) + c110 * xd
        # c11 = c011 * (1 - xd) + c111 * xd

        # c0 = c00 * (1 - yd) + c10 * yd
        # c1 = c01 * (1 - yd) + c11 * yd

        # c = c0 * (1 - zd) + c1 * zd

    return c000 if enc=='z3' else c000


def apply_interpolation(x, y, z, values, enc, res, L, index_factor, grid_size=41):
    n = find_lattice_neighbors_1d(x, y, z, enc, res, L, index_factor, grid_size)
    return trilinear_interpolation(values, n, enc)
