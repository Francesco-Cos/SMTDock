from z3 import *

def find_lattice_neighbors_1d(x, y, z, enc, grid_size=41, grid_spacing=0.375):
    L = (grid_size * grid_spacing) / 2  # Half the total length, to shift coordinates

    # Compute integer lattice indices, shifting for a centered grid
    if enc == 'z3':
        i = ToInt((x + L) / grid_spacing)
        j = ToInt((y + L) / grid_spacing)
        k = ToInt((z + L) / grid_spacing)
    else:
        i = int((x + L) // grid_spacing)
        j = int((y + L) // grid_spacing)
        k = int((z + L) // grid_spacing)

    # Compute 1D indices for the 8 neighboring points
    def index(ix, iy, iz):
        return ix + grid_size * (iy + grid_size * iz)

    neighbors = [
        index(i, j, k)
        # , index(i+1, j, k), index(i, j+1, k), index(i, j, k+1),
        # index(i+1, j+1, k), index(i+1, j, k+1), index(i, j+1, k+1), index(i+1, j+1, k+1)
    ]
    # if enc == 'z3':
    #     A = Array('A', RealSort(), RealSort())
    #     ii = 0
    #     for elem in neighbors:
    #         A = Store(A, ii, elem)
    #         ii = ii + 1
    #     neighbors = A

    return neighbors

def trilinear_interpolation(x, y, z, values, neighbors, enc, grid_spacing=0.375):
    # Compute integer lattice indices
    if enc == 'z3':
        i = ToInt((x / grid_spacing))
        j = ToInt((y / grid_spacing))
        k = ToInt((z / grid_spacing))
    else:
        i = int((x // grid_spacing))
        j = int((y // grid_spacing))
        k = int((z // grid_spacing))
    
    # Compute interpolation weights
    xd = (x / grid_spacing) - i
    yd = (y / grid_spacing) - j
    zd = (z / grid_spacing) - k

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
        # c100 = values[neighbors[1]]
        # c010 = values[neighbors[2]]
        # c001 = values[neighbors[3]]
        # c110 = values[neighbors[4]]
        # c101 = values[neighbors[5]]
        # c011 = values[neighbors[6]]
        # c111 = values[neighbors[7]]
    
    # Perform trilinear interpolation
    # c00 = c000 * (1 - xd) + c100 * xd
    # c01 = c001 * (1 - xd) + c101 * xd
    # c10 = c010 * (1 - xd) + c110 * xd
    # c11 = c011 * (1 - xd) + c111 * xd
    
    # c0 = c00 * (1 - yd) + c10 * yd
    # c1 = c01 * (1 - yd) + c11 * yd
    
    # c = c0 * (1 - zd) + c1 * zd
    
    return ToInt(c000) if enc == 'z3' else c000

def apply_interpolation(x, y, z, values, enc):
    n = find_lattice_neighbors_1d(x, y, z, enc)
    return trilinear_interpolation(x, y, z, values, n, enc)
