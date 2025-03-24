from z3 import *
import numpy as np
from util.tri_linear_interpolation import apply_interpolation

center = {{{{center}}}}
spacing = {{{{spacing}}}}
atomst = {{{{atomst}}}}
maps = {}
for atomt in atomst:
    if atomt not in maps.keys():
        B = Array ('B', RealSort(), RealSort())
        with open(f'a2a_adenosine/input/a2a_h.{atomt}.map') as m:
            i = 0
            for line in m:
                try:
                    B = Store(B, i, float(line.strip()))
                    i = i + 1
                except ValueError:
                    continue
        maps[atomt] = B

# Variables def
sina, cosa = Reals('sina cosa')
sinb, cosb = Reals('sinb cosb')
sing, cosg = Reals('sing cosg')
energy = Int('energy')
x0, xx0, y0, yy0, q0, qq0 = Reals('x0 xx0 y0 yy0 q0 qq0')
x1, xx1, y1, yy1, q1, qq1 = Reals('x1 xx1 y1 yy1 q1 qq1')
x2, xx2, y2, yy2, q2, qq2 = Reals('x2 xx2 y2 yy2 q2 qq2')

{{{{variables}}}}

solver = Solver()

# Constraints def
solver.add(
    # -1 <= cost, cost <= 1,
    -1 <= sina, sina <= 1,
    cosa**2 == 1 - sina**2,
    -1 <= sinb, sinb <= 1,
    cosb**2 == 1 - sinb**2,
    -1 <= sing, sing <= 1,
    cosg**2 == 1 - sing**2,

    {{{{space_constraints}}}}
    
    {{{{anchors_distance_constraints}}}}

    {{{{full_points_constraints}}}}

    {{{{rotation_constraints}}}}
    energy == {{{{energy}}}}
    # energy < 100000
)

# print('solving')
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
