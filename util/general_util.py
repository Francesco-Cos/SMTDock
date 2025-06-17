def generate_distance(p1, p2):
        diff = [e1 - e2 for e1, e2 in zip(p1, p2)]
        return sum([d**2 for d in diff])

def int_coord_distance(p1, p2, spacing):
        diff = [(e1 - e2)//spacing for e1, e2 in zip(p1, p2)]
        return diff

def center_of_mass(points):
    com = [0,0,0]
    for point in points:
        for i, coords in enumerate(point):
            com[i] += coords
    com = [el/len(points) for el in com]
    return com

def get_cube_bounds(iteration, num_iterations, bounds):
    n = round(num_iterations ** (1/3))
    assert n ** 3 == num_iterations, "num_iterations must be a perfect cube"

    x_min = y_min = z_min = -bounds
    x_max = y_max = z_max = bounds

    dx = (x_max - x_min) / n
    dy = (y_max - y_min) / n
    dz = (z_max - z_min) / n

    ix = iteration % n
    iy = (iteration // n) % n
    iz = (iteration // (n * n)) % n

    x0 = x_min + ix * dx
    x1 = x0 + dx
    y0 = y_min + iy * dy
    y1 = y0 + dy
    z0 = z_min + iz * dz
    z1 = z0 + dz

    return (x0, x1), (y0, y1), (z0, z1)