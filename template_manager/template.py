from z3 import *
import numpy as np
import matplotlib.pyplot as plt
from util.tri_linear_interpolation import apply_interpolation
from util.downsample import downsample
from util.general_util import center_of_mass
import time, json

set_option('smt.array.extensional', False)

def get_cube_bounds(iteration, num_iterations, bounds):
    import math

    # Ensure perfect cube root
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

grid_size = {{{{grid_size}}}}
new_size = {{{{new_size}}}}
grid_spacing = {{{{grid_spacing}}}}
res = {{{{res}}}}
new_spacing = grid_spacing * res * (grid_size / new_size)
center = [0, 0, 0]
space_bounds = {{{{space_bounds}}}}
index_factor = new_size/(space_bounds*2)
atomst = {{{{atomst}}}}
maps = {}
for atomt in atomst:
    if atomt not in maps.keys():
        B = Array('B', IntSort(), IntSort())
        with open(f'a2a_adenosine/input/a2a_h.{atomt}.map') as m:
            original_lattice = m.readlines()[6:]
        new_lattice = downsample(original_lattice, grid_size+1, new_size)
        with open(f'new_map_{atomt}.map', 'w') as f:
            f.write('\n'.join([f'{num}' for num in new_lattice]))
        i = 0
        for value in new_lattice:
            B = Store(B, i, int(float(value)))
            i = i + 1
        maps[atomt] = B

{{{{energy_variables}}}}
energy1 = Int('energy1')

{{{{trans_variables}}}}

trans_energy = np.inf
final_energy = np.inf
final_points = []
previous_models = []
number_models = 27
coms = []
for i in range(number_models):
    red_x, red_y, red_q = get_cube_bounds(i, number_models, space_bounds)
    com = center_of_mass([(x0,y0,q0), (x1,y1,q1), (x2,y2,q2)])

    solver = Solver()

    solver.add(
        {{{{space_constraints}}}}

        {{{{anchors_distance_constraints}}}}

        {{{{full_points_constraints}}}}

        {{{{energy}}}}
    )
    # for model in previous_models:
    #     solver.add(Or(model))
    result = solver.check()

    if result == sat:
        model = solver.model()
        selected_points = [
            {{{{returned_points}}}}
        ]
        energy = model.evaluate(energy1).as_long()
        com_dup = [Sum([selected_points[i][j] for i in range(3)])//3 for j in range(3)]
        coms.append(com_dup)
        previous_models.append([{{{{remove_duplicate}}}}])
    else:
        break
        print("No possible ligand position")
        selected_points = []

    com = [0,0,0]

    for point in selected_points:
        com = [cord + point[i] for i,cord in enumerate(com)]

    com = [cord//len(selected_points) for cord in com]

    {{{{rot_variables}}}}
    energy2 = Int('energy2')


    conds = [
        {{{{rotation_constraints}}}}
    ]

    saved_model = None
    best_energy = np.inf

    for i,cond in enumerate(conds):
        solver2 = Solver()
        solver2.add(
            {{{{space_constraints2}}}}

            {{{{space_constraints_all2}}}}

            cond,
            {{{{energy2}}}}
        )

        start = time.time()
        result = solver2.check()
        end = time.time()

        if result == sat:
            model = solver2.model()
            if model.evaluate(energy2).as_long() <= best_energy:
                saved_model = model
                best_energy = model.evaluate(energy2).as_long()
        else:
            continue
            selected_points = []
    if saved_model:
        selected_points = [
            {{{{returned_points2}}}}
        ]
        if best_energy < final_energy:
            final_energy = best_energy
            final_points = selected_points

print(coms)
print(final_points)