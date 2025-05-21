import numpy as np

def downsample_average(original_lattice, orig_size, new_size):
    
    orig_shape = (orig_size, orig_size, orig_size)
    new_shape = (new_size, new_size, new_size)

    # Output array (flattened)
    new_lattice = np.zeros(np.prod(new_shape))

    # Helper: convert 3D -> 1D index
    def ravel_index(x, y, z, shape):
        return x * shape[1] * shape[2] + y * shape[2] + z

    # Compute boundaries of bins in the original grid
    # Use linspace to divide the 41 bins into 11 equal regions
    bin_edges = np.linspace(0, orig_shape[0], new_shape[0]+1, dtype=int)

    for i in range(new_shape[0]):
        for j in range(new_shape[1]):
            for k in range(new_shape[2]):
                # Define boundaries of the box in the original grid
                x0, x1 = bin_edges[i], bin_edges[i+1]
                y0, y1 = bin_edges[j], bin_edges[j+1]
                z0, z1 = bin_edges[k], bin_edges[k+1]

                # Collect all the values in that box
                values = []
                for xi in range(x0, x1):
                    for yj in range(y0, y1):
                        for zk in range(z0, z1):
                            flat_idx = ravel_index(xi, yj, zk, orig_shape)
                            values.append(int(float(original_lattice[flat_idx].strip())))

                # Average and store in new lattice
                new_idx = ravel_index(i, j, k, new_shape)
                new_lattice[new_idx] = np.mean(values)
    return new_lattice