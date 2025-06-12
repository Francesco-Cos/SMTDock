import re
from typing import Dict
import util.general_util as utl

class IntSolver:
    def __init__(self, space_bounds, size, new_size, spacing, res, anchor_points, ligand_points, atomst):
        self._space_bounds = space_bounds
        self._size = size
        self._new_size = new_size
        self._spacing = spacing/res
        self._orig_spacing = spacing
        self._res = res
        self._anchor_points = anchor_points
        self._ligand_points = ligand_points
        self._atomst = atomst

    def _update_builder(self, builder):
        builder.update(grid_size=self._size)
        builder.update(new_size=self._new_size)
        builder.update(grid_spacing=self._spacing)
        builder.update(space_bounds=self._space_bounds)
        builder.update(res=self._res)
        builder.update(atomst=self._atomst)
        builder.update(trig_variables='')

        builder.update(trans_variables='\n'.join(
            f"x{i}, y{i}, q{i}= Ints('x{i} y{i} q{i}')"
            for i in range(0, len(self._ligand_points) + len(self._anchor_points))
        ))

        builder.update(rot_variables='\n'.join(
            f"xx{i}, yy{i}, qq{i} = Ints('xx{i} yy{i} qq{i}')"
            for i in range(0, len(self._ligand_points) + len(self._anchor_points))
        ))

        builder.update(energy_variables='')

        builder.update(space_constraints=',\n'.join(
            # ', '.join(
                f"center[{j}]+red_{d}[0] <= com[{j}], com[{j}] < center[{j}]+red_{d}[1]"
                for j, d in enumerate(['x', 'y', 'q'])
            # ) for i in range(len(self._anchor_points))
        ) + ','
        )

        builder.update(space_constraints2=',\n'.join(
            ', '.join(
                f"center[{j}]-space_bounds <= {d}{i}, {d}{i} < center[{j}]+space_bounds"
                for j, d in enumerate(['xx', 'yy', 'qq'])
            ) for i in range(len(self._anchor_points))
        ) + ','
        )

        # builder.update(space_constraints_all=',\n'.join(
        #     ', '.join(
        #         f"center[{j}]-space_bounds <= {d}{i}, {d}{i} < center[{j}]+space_bounds"
        #         for j, d in enumerate(['x', 'y', 'q'])
        #     ) for i in range(len(self._anchor_points), len(self._anchor_points) + len(self._ligand_points))
        # ) + ','
        # )

        builder.update(space_constraints_all2=',\n'.join(
            ', '.join(
                f"center[{j}]-space_bounds <= {d}{i}, {d}{i} < center[{j}]+space_bounds"
                for j, d in enumerate(['xx', 'yy', 'qq'])
            ) for i in range(len(self._anchor_points), len(self._anchor_points) + len(self._ligand_points))
        ) + ','
        )

        anchor_distances = []
        for i, pa in enumerate(self._anchor_points):
            anchor_distances.append(utl.int_coord_distance(pa, self._anchor_points[(i+1) % len(self._anchor_points)], self._spacing))

        builder.update(anchors_distance_constraints=',\n'.join(
            ', '.join(
                f"({d}{i} - {d}{(i+1) % len(self._anchor_points)}) >= {int(dist[j]) - 1}, ({d}{i} - {d}{(i+1) % len(self._anchor_points)}) <= {int(dist[j]) + 1}"
                for j, d in enumerate(['x', 'y', 'q'])
            )
            for i, dist in enumerate(anchor_distances)
        ) + ','
        )

        points_distances = []
        for j, pl in enumerate(self._ligand_points):
            points_distances.append([])
            for i, pa in enumerate(self._anchor_points):
                points_distances[j].append(utl.int_coord_distance(pa, pl, self._spacing))

        builder.update(full_points_constraints=',\n'.join(
            ',\n'.join(
                ', '.join(
                    f"({d}{g} - {d}{i + len(self._anchor_points)}) >= {int(dist[g][j])-1}, ({d}{g} - {d}{i + len(self._anchor_points)}) <= {int(dist[g][j])+1}"
                    for j, d in enumerate(['x', 'y', 'q'])
                )
                for g, _ in enumerate(self._anchor_points)
            )
            for i, dist in enumerate(points_distances)
        ) + ','
        )

        def get_rotations():
            return [
                ( 'selected_points[*][0]-com[0]',  'selected_points[*][1]-com[1]', ' selected_points[*][2]-com[2]'),   ( 'selected_points[*][0]-com[0]', '-selected_points[*][1]-com[1]', '-selected_points[*][2]-com[2]'),  ('-selected_points[*][0]-com[0]',  'selected_points[*][1]-com[1]', '-selected_points[*][2]-com[2]'),  ('-selected_points[*][0]-com[0]', '-selected_points[*][1]-com[1]',  'selected_points[*][2]-com[2]'),
                ( 'selected_points[*][1]-com[1]',  'selected_points[*][2]-com[2]', ' selected_points[*][0]-com[0]'),   ( 'selected_points[*][1]-com[1]', '-selected_points[*][2]-com[2]', '-selected_points[*][0]-com[0]'),  ('-selected_points[*][1]-com[1]',  'selected_points[*][2]-com[2]', '-selected_points[*][0]-com[0]'),  ('-selected_points[*][1]-com[1]', '-selected_points[*][2]-com[2]',  'selected_points[*][0]-com[0]'),
                ( 'selected_points[*][2]-com[2]',  'selected_points[*][0]-com[0]', ' selected_points[*][1]-com[1]'),   ( 'selected_points[*][2]-com[2]', '-selected_points[*][0]-com[0]', '-selected_points[*][1]-com[1]'),  ('-selected_points[*][2]-com[2]',  'selected_points[*][0]-com[0]', '-selected_points[*][1]-com[1]'),  ('-selected_points[*][2]-com[2]', '-selected_points[*][0]-com[0]',  'selected_points[*][1]-com[1]'),
                ( 'selected_points[*][0]-com[0]',  'selected_points[*][2]-com[2]', '-selected_points[*][1]-com[1]'),   ( 'selected_points[*][0]-com[0]', '-selected_points[*][2]-com[2]',  'selected_points[*][1]-com[1]'),  ('-selected_points[*][0]-com[0]',  'selected_points[*][2]-com[2]',  'selected_points[*][1]-com[1]'),  ('-selected_points[*][0]-com[0]', '-selected_points[*][2]-com[2]', '-selected_points[*][1]-com[1]'),
                ( 'selected_points[*][1]-com[1]',  'selected_points[*][0]-com[0]', '-selected_points[*][2]-com[2]'),   ( 'selected_points[*][1]-com[1]', '-selected_points[*][0]-com[0]',  'selected_points[*][2]-com[2]'),  ('-selected_points[*][1]-com[1]',  'selected_points[*][0]-com[0]',  'selected_points[*][2]-com[2]'),  ('-selected_points[*][1]-com[1]', '-selected_points[*][0]-com[0]', '-selected_points[*][2]-com[2]'),
                ( 'selected_points[*][2]-com[2]',  'selected_points[*][1]-com[1]', '-selected_points[*][0]-com[0]'),   ( 'selected_points[*][2]-com[2]', '-selected_points[*][1]-com[1]',  'selected_points[*][0]-com[0]'),  ('-selected_points[*][2]-com[2]',  'selected_points[*][1]-com[1]',  'selected_points[*][0]-com[0]'),  ('-selected_points[*][2]-com[2]', '-selected_points[*][1]-com[1]', '-selected_points[*][0]-com[0]'),
            ]

        rotations = get_rotations()

        builder.update(rotation_constraints='),\n'.join(
            'And(' + ', '.join(
                f'xx{i} == {rot[0].replace('*', f'{i}')}+com[0], yy{i} == {rot[1].replace('*', f'{i}')}+com[1], qq{i} == {rot[2].replace('*', f'{i}')}+com[2]'
                for i in range(len(self._anchor_points) + len(self._ligand_points))
            )
            for rot in rotations
        ) + ')'
        )

        builder.update(energy='energy1 == ' + '+'.join(
            f'apply_interpolation(x{i}, y{i}, q{i}, maps[atomst[{i}]], "z3", space_bounds, index_factor, new_size)'
            for i in range(len(self._anchor_points) + len(self._ligand_points))
        ) + ','
        )

        builder.update(energy2='energy2 == ' + '+'.join(
            f'apply_interpolation(xx{i}, yy{i}, qq{i}, maps[atomst[{i}]], "z3", space_bounds, index_factor, new_size)'
            for i in range(len(self._anchor_points) + len(self._ligand_points))
        ) + ','
        )

        builder.update(remove_duplicate=', '.join(
            f'x{i} != selected_points[{i}][0], y{i} != selected_points[{i}][1], q{i} != selected_points[{i}][2]'
            for i in range(len(self._anchor_points))
        )
        )

        # builder.update(remove_duplicate=', \n\t\t\t'.join(
        #     f'(x{i} - com_dup[0])**2 + (y{i} - com_dup[1])**2 + (q{i} - com_dup[2])**2 > (selected_points[0][0] - com_dup[0])**2 + (selected_points[0][1] - com_dup[1])**2 + (selected_points[0][2] - com_dup[2])**2' 
        #     for i in range(len(self._anchor_points))
        # )
        # )

        builder.update(returned_points=',\n'.join('(' +
                                                  ','.join(
                                                      f"(model.evaluate({d}{i}).as_long())"
                                                      for d in ['x', 'y', 'q'])+')'
                                                  for i in range(len(self._ligand_points) + len(self._anchor_points))
        )
        )

        builder.update(returned_points2=',\n'.join('(' +
                                                  ','.join(
                                                      f"(saved_model.evaluate({d}{i}))"
                                                      for d in ['xx', 'yy', 'qq'])+')'
                                                  for i in range(len(self._ligand_points) + len(self._anchor_points))
        )
        )

        builder.update(print_energy=','.join(
            f'print(model.evaluate(energy{i}), atomst[{i}])'
            for i in range(len(self._anchor_points) + len(self._ligand_points))
        )
        )

    def generate_script(self):
        builder = Builder.from_file('template_manager/template.py')

        self._update_builder(builder)

        with open('solver_script.py', 'w') as ofile:
            ofile.write(builder.finalize())


class Builder:
    LEFT_DELIMITER = '{{{{'
    RIGHT_DELIMITER = '}}}}'
    MAGIC_STRING_1 = 'zqwezrtyzuiozpaszdfgzhjkzlzxzcvbznm1z234z567z890z'
    MAGIC_STRING_2 = 'zmnbzvcxzzlkzjhgzfdszapoziuyztrezwq0z987z654z321z'

    def __init__(self, string: str) -> 'Builder':
        self._string: str = string
        self._kwargs: Dict[str, str] = dict()

    @classmethod
    def from_file(cls, filename: str):
        with open(filename, 'r') as ifile:
            return cls(ifile.read())

    def update(self, **kwargs: Dict[str, str]) -> None:
        self._kwargs.update(kwargs)

    def finalize(self) -> str:
        # get normalized template string
        normalized_string = (
            self._string
            .replace(self.LEFT_DELIMITER, self.MAGIC_STRING_1)
            .replace(self.RIGHT_DELIMITER, self.MAGIC_STRING_2)
            .replace('{', '{{')
            .replace('}', '}}')
            .replace(self.MAGIC_STRING_1, '{')
            .replace(self.MAGIC_STRING_2, '}')
        )

        # update kwargs with correctly tabulated values
        for key, value in self._kwargs.items():
            if m := re.search(rf'(?:\r\n|\r|\n)([\t ]+)\{{{key}\}}', normalized_string):
                tabulation = m.group(1)
                self._kwargs[key] = tabulation.join(value.splitlines(True))

        # apply kwargs to the tamplate
        return normalized_string.format(**self._kwargs)

    def __str__(self) -> str:
        return f'{self._string!r} <- {self._kwargs}'
