from z3 import *
import numpy as np
from util.tri_linear_interpolation import apply_interpolation

set_option('smt.array.extensional', False)

res = {{{{res}}}}
center = {{{{center}}}}
bounds = {{{{bounds}}}}
atomst = {{{{atomst}}}}
maps = {}
for atomt in atomst:
    if atomt not in maps.keys():
        B = Array('B', IntSort(), IntSort())
        with open(f'a2a_adenosine/input/a2a_h.{atomt}.map') as m:
            i = 0
            # for line in m.readlines()[6:]:
            #     B = Store(B, i, int(float(line.strip())))
            #     i = i + 1
            for line in m:
                try:
                    B = Store(B, i, int(line.strip()))
                    i = i + 1
                except ValueError:
                    continue
        maps[atomt] = B

# Variables def
{{{{trig_variables}}}}
energy = Int('energy')

{{{{variables}}}}

solver = Solver()

# Constraints def
solver.add(
    {{{{trig_constraints}}}}

    {{{{space_constraints}}}}

    {{{{anchors_distance_constraints}}}}

    {{{{full_points_constraints}}}}

    {{{{rotation_constraints}}}}
    energy == {{{{energy}}}}
    # energy < 10
)

# print('solving')
# solver.minimize(energy)
result = solver.check()

if result == sat:
    model = solver.model()
    selected_points = [
        {{{{returned_points}}}}
    ]
    print(selected_points)
    # print(model.evaluate(energy))
else:
    print("No possible ligand position")
    selected_points = []
