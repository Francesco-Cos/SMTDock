import re
import itertools
from typing import Dict


class RealSolver:
    def __init__(self, center, bounds, res, anchor_points, ligand_points, atomst):
        self._center = center
        self._bounds = bounds
        self._res = res
        self._anchor_points = anchor_points
        self._ligand_points = ligand_points
        self._atomst = atomst

    def _generate_distance(self, p1, p2):
        diff = [e1 - e2 for e1, e2 in zip(p1, p2)]
        return sum([d**2 for d in diff])

    def _coord_distance(self, p1, p2):
        diff = [e1 - e2 for e1, e2 in zip(p1, p2)]
        return diff

    def _update_builder(self, builder):
        builder.update(center=self._center)
        builder.update(bounds=self._bounds)
        builder.update(res=self._res)
        builder.update(atomst=self._atomst)

        builder.update(trig_variables='\n'.join(
            f"sin{ang}, cos{ang} = Reals('sin{ang} cos{ang}')"
            for ang in ['a', 'b', 'g']
        ))

        builder.update(variables='\n'.join(
            f"x{i}, xx{i}, y{i}, yy{i}, q{i}, qq{i} = Reals('x{i} xx{i} y{i} yy{i} q{i} qq{i}')"
            for i in range(0, len(self._ligand_points) + len(self._anchor_points))
        ))

        builder.update(trig_constraints='\n'.join(
            f"-1 <= sin{ang}, sin{ang} <= 1,\ncos{ang}**2 == 1 - sin{ang}**2,"
            for ang in ['a', 'b', 'g']
        ))

        builder.update(space_constraints=',\n'.join(
            ', '.join(
                f"center[{j}]-bounds[{j}] <= {d}{i}, {d}{i} < center[{j}]+bounds[{j}]"
                for j, d in enumerate(['x', 'y', 'q'])
            ) for i in range(len(self._anchor_points))
        ) + ','
        )

        anchor_distances = []
        for i, pa in enumerate(self._anchor_points):
            anchor_distances.append(self._coord_distance(pa, self._anchor_points[(i+1) % len(self._anchor_points)]))

        builder.update(anchors_distance_constraints=',\n'.join(
            ', '.join(
                f"({d}{i} - {d}{(i+1) % len(self._anchor_points)}) == {dist[j]:.3f}"
                for j, d in enumerate(['x', 'y', 'q'])
            )
            for i, dist in enumerate(anchor_distances)
        ) + ','
        )

        points_distances = []
        for j, pl in enumerate(self._ligand_points):
            points_distances.append([])
            for i, pa in enumerate(self._anchor_points):
                points_distances[j].append(self._coord_distance(pa, pl))
        builder.update(full_points_constraints=',\n'.join(
            ',\n'.join(
                ', '.join(
                    f"({d}{g} - {d}{i + len(self._anchor_points)}) == {dist[g][j]:.3f}"
                    for j, d in enumerate(['x', 'y', 'q'])
                )
                for g, pa in enumerate(self._anchor_points)
            )
            for i, dist in enumerate(points_distances)
        ) + ','
        )

        builder.update(rotation_constraints=',\n'.join(
            f"xx{i} == x{i}*cosb*cosg + y{i}*(sina*sinb*cosg - cosa*sing) + q{i}*(cosa*sinb*cosg + sina*sing), " +
            f"yy{i} == x{i}*cosb*sing + y{i}*(sina*sinb*sing + cosa*cosg) + q{i}*(cosa*sinb*sing - sina*cosg), " +
            f"qq{i} == x{i}*(-sinb) + y{i}*sina*cosb + q{i}*cosa*cosb"
            for i in range(len(self._anchor_points) + len(self._ligand_points))
        ) + ','
        )

        builder.update(energy=' + '.join(
            f'apply_interpolation(xx{i}, yy{i}, qq{i}, maps[atomst[{i}]], "z3", res)'
            for i in range(len(self._anchor_points) + len(self._ligand_points))
        ) + ','
        )

        builder.update(returned_points=',\n'.join('(' +
                                                  ','.join(
                                                      f"float(model.evaluate({d}{d}{i}).as_decimal(10).strip('?'))"
                                                      for d in ['x', 'y', 'q'])+')'
                                                  for i in range(len(self._ligand_points) + len(self._anchor_points)))
                       )

    def generate_script(self):
        builder = Builder.from_file('template_manager/template.py')

        self._update_builder(builder)

        with open('solver_script.py', 'w') as ofile:
            ofile.write(builder.finalize())

class IntSolver:
    def __init__(self, center, bounds, spacing, res, anchor_points, ligand_points, atomst):
        self._center = center
        self._bounds = bounds
        self._spacing = spacing
        self._res = res
        self._anchor_points = anchor_points
        self._ligand_points = ligand_points
        self._atomst = atomst

    def _generate_distance(self, p1, p2):
        diff = [e1 - e2 for e1, e2 in zip(p1, p2)]
        return sum([d**2 for d in diff])

    def _coord_distance(self, p1, p2):
        diff = [e1 - e2 for e1, e2 in zip(p1, p2)]
        return diff
    
    def _int_coord_distance(self, p1, p2):
        diff = [(e1 - e2)//self._spacing for e1, e2 in zip(p1, p2)]
        # diff = diff // self._dspacing
        return diff

    def _update_builder(self, builder):
        builder.update(center=self._center)
        builder.update(bounds=self._bounds)
        builder.update(res=self._res)
        builder.update(atomst=self._atomst)

        builder.update(trig_variables='')

        builder.update(variables='\n'.join(
            f"x{i}, y{i}, q{i}, = Ints('x{i} y{i} q{i}')"
            for i in range(0, len(self._ligand_points) + len(self._anchor_points))
        ))

        builder.update(trig_constraints='')

        builder.update(space_constraints=',\n'.join(
            ', '.join(
                f"center[{j}]-bounds[{j}] <= {d}{i}, {d}{i} < center[{j}]+bounds[{j}]"
                for j, d in enumerate(['x', 'y', 'q'])
            ) for i in range(len(self._anchor_points))
        ) + ','
        )

        anchor_distances = []
        for i, pa in enumerate(self._anchor_points):
            anchor_distances.append(self._int_coord_distance(pa, self._anchor_points[(i+1) % len(self._anchor_points)]))

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
                points_distances[j].append(self._int_coord_distance(pa, pl))
        builder.update(full_points_constraints=',\n'.join(
            ',\n'.join(
                ', '.join(
                    f"({d}{g} - {d}{i + len(self._anchor_points)}) >= {int(dist[g][j])-1}, ({d}{g} - {d}{i + len(self._anchor_points)}) <= {int(dist[g][j])+1}"
                    for j, d in enumerate(['x', 'y', 'q'])
                )
                for g, pa in enumerate(self._anchor_points)
            )
            for i, dist in enumerate(points_distances)
        ) + ','
        )

        builder.update(rotation_constraints='')

        builder.update(energy=' + '.join(
            f'apply_interpolation(x{i}, y{i}, q{i}, maps[atomst[{i}]], "z3", res)'
            for i in range(len(self._anchor_points) + len(self._ligand_points))
        ) + ','
        )

        builder.update(returned_points=',\n'.join('(' +
                                                  ','.join(
                                                      f"(model.evaluate({d}{i}))"
                                                      for d in ['x', 'y', 'q'])+')'
                                                  for i in range(len(self._ligand_points) + len(self._anchor_points)))
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
