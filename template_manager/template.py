from z3 import *
import numpy as np
from util.tri_linear_interpolation import apply_interpolation
from util.downsample import downsample
import time

set_option('smt.array.extensional', False)

grid_size = {{{{grid_size}}}}
new_size = {{{{new_size}}}}
grid_spacing = {{{{grid_spacing}}}}
res = {{{{res}}}}
new_spacing = grid_spacing * res * (grid_size / new_size)
center = {{{{center}}}}
bounds = {{{{bounds}}}}
space_bounds = {{{{space_bounds}}}}
index_factor = new_size/(space_bounds*2)
atomst = {{{{atomst}}}}
maps = {}
for atomt in atomst:
    if atomt not in maps.keys():
        B = Array('B', IntSort(), IntSort())
        with open(f'a2a_adenosine/input/a2a_h.{atomt}.map') as m:
            original_lattice = m.readlines()[6:]
        new_lattice = downsample(original_lattice, grid_size+1, new_size+1)
        with open(f'new_map_{atomt}.map', 'w') as f:
            f.write('\n'.join([f'{num}' for num in new_lattice]))
        i = 0
        for value in new_lattice:
            B = Store(B, i, int(float(value)))
            i = i + 1
            # for line in m:
            #     try:
            #         B = Store(B, i, int(line.strip()))
            #         i = i + 1
            #     except ValueError:
            #         continue
        maps[atomt] = B

# Variables def
{{{{energy_variables}}}}
energy1 = Int('energy1')

{{{{trans_variables}}}}

# print('solving 1')
trans_energy = np.inf
# solver.minimize(energy1)
final_energy = np.inf
final_points = []
previous_models = []
number_models = 10
for i in range(number_models):

    solver = Solver()

    # Constraints def
    solver.add(
        {{{{space_constraints}}}}

        {{{{space_constraints_all}}}}

        {{{{anchors_distance_constraints}}}}

        {{{{full_points_constraints}}}}

        {{{{energy_constraints}}}}

        {{{{energy}}}}
        
        # energy < 0
    )
    for model in previous_models:
        solver.add(And(model))
    result = solver.check()
    # print(end - start)

    if result == sat:
        model = solver.model()
        selected_points = [
            {{{{returned_points}}}}
        ]
        # {{{{print_energy}}}}
        energy = model.evaluate(energy1).as_long()
        com_dup = [Sum([selected_points[i][j] for j in range(3)])//3 for i in range(3)]
        previous_models.append([{{{{remove_duplicate}}}}])
    else:
        break
        print("No possible ligand position")
        selected_points = []

    com = [0,0,0]

    for point in selected_points:
        com = [cord + point[i] for i,cord in enumerate(com)]

    com = [cord//len(selected_points) for cord in com]

    # bounding = {'min': 0, 'max': 0}
    maximal_bounding = 0

    for point in selected_points:
        for i,cord in enumerate(point):
            if abs(com[i] - cord) > maximal_bounding:
                maximal_bounding = abs(com[i] - cord)

    # {{{{trig_variables}}}}
    {{{{rot_variables}}}}
    energy2 = Int('energy2')


    conds = [
        {{{{rotation_constraints}}}}
    ]

    saved_model = None
    best_energy = energy

    for i,cond in enumerate(conds):
        solver2 = Solver()
        # Constraints def
        solver2.add(
            {{{{energy_constraints}}}}
            cond,
            {{{{energy2}}}}

            # energy < 0
        )

        # print('solving 2')
        start = time.time()
        # solver2.minimize(energy2)
        result = solver2.check()
        end = time.time()
        # print(end - start)

        if result == sat:
            model = solver2.model()
            if model.evaluate(energy2).as_long() <= best_energy:
                # print(f'better energy: {best_energy} -> {model.evaluate(energy2).as_long()}')
                saved_model = model
                best_energy = model.evaluate(energy2).as_long()
            # {{{{print_energy}}}}
            # print(model.evaluate(energy2))
        else:
            print("No possible ligand position")
            selected_points = []

    selected_points = [
        {{{{returned_points2}}}}
    ]
    if best_energy < final_energy:
        final_energy = best_energy
        final_points = selected_points
    # print(selected_points)
print(final_points)