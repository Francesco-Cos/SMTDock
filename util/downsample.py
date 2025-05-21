import numpy as np

def downsample(original_lattice, orig_size, new_size):
    orig_shape = (orig_size, orig_size, orig_size)
    new_shape = (new_size, new_size, new_size)

    new_lattice = np.zeros(np.prod(new_shape))

    def ravel_index(x, y, z, shape):
        return x * shape[1] * shape[2] + y * shape[2] + z

    bin_edges = np.linspace(0, orig_shape[0], new_shape[0] + 1, dtype=int)

    for i in range(new_shape[0]):
        for j in range(new_shape[1]):
            for k in range(new_shape[2]):
                x0, x1 = bin_edges[i], bin_edges[i+1]
                y0, y1 = bin_edges[j], bin_edges[j+1]
                z0, z1 = bin_edges[k], bin_edges[k+1]

                # Compute the center index of the current bin
                center_x = (x0 + x1 - 1) // 2
                center_y = (y0 + y1 - 1) // 2
                center_z = (z0 + z1 - 1) // 2

                center_idx = ravel_index(center_x, center_y, center_z, orig_shape)
                new_idx = ravel_index(i, j, k, new_shape)

                new_lattice[new_idx] = int(float(original_lattice[center_idx].strip()))
    
    return new_lattice
