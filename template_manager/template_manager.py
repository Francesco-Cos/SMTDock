import re
import itertools
from typing import Dict

class RealSolver:
    def __init__(self, center, spacing, anchor_points, ligand_points, atomst):
        self._center = center
        self._spacing = spacing
        self._anchor_points = anchor_points
        self._ligand_points = ligand_points
        self._atomst = atomst

    def _generate_distance(self, p1, p2):
        diff = [e1 - e2 for e1,e2 in zip(p1,p2)]
        return sum([d**2 for d in diff])
    
    def _coord_distance(self, p1, p2):
        diff = [e1 - e2 for e1,e2 in zip(p1,p2)]
        return diff

    def _update_builder(self, builder):
        builder.update(center=self._center)
        builder.update(spacing=self._spacing)
        builder.update(atomst=self._atomst)

        builder.update(variables='\n'.join(
            f"x{i}, xx{i}, y{i}, yy{i}, q{i}, qq{i} = Reals('x{i} xx{i} y{i} yy{i} q{i} qq{i}')"
            for i in range(len(self._anchor_points), len(self._ligand_points) + len(self._anchor_points))
        ))

        # solver.add(center[0]-spacing[0] <= x1, x1 < center[0]+spacing[0])
        # solver.add(center[1]-spacing[1] <= y1, y1 < center[1]+spacing[1])
        # solver.add(center[2]-spacing[2] <= q1, q1 < center[2]+spacing[2])
        builder.update(space_constraints=',\n'.join(
            ', '.join(
                f"center[{j}]-spacing[{j}] <= {d}{i}, {d}{i} < center[{j}]+spacing[{j}]"
                for j,d in enumerate(['x','y','q'])
            ) for i in range(len(self._anchor_points))
        ) + ','
        )

        # (x1 - x2)**2 + (y1 - y2)**2 + (q1 - q2)**2 == d12**2
        # builder.update(anchors_distance_constraints=', '.join(
        #     ' + '.join(
        #         f"({d}{i} - {d}{(i+1)%len(self._anchor_points)})**2"
        #         for d in ['x', 'y', 'q']
        #     ) + f" == {self._generate_distance(pa,self._anchor_points[(i+1)%len(self._anchor_points)]):.3f}"
        #     for i,pa in enumerate(self._anchor_points)
        # ) + ','
        # )
        anchor_distances = []
        for i,pa in enumerate(self._anchor_points):
            anchor_distances.append(self._coord_distance(pa,self._anchor_points[(i+1)%len(self._anchor_points)]))

        builder.update(anchors_distance_constraints=',\n'.join(
            ', '.join(
                f"({d}{i} - {d}{(i+1)%len(self._anchor_points)}) == {dist[j]:.3f}"
                for j, d in enumerate(['x', 'y', 'q'])
            )
        for i,dist in enumerate(anchor_distances)
        ) + ','
        )

        # (x1 - x2)**2 + (y1 - y2)**2 + (q1 - q2)**2 == d12**2
        # builder.update(full_points_constraints=',\n'.join(
        #     ', '.join(
        #         ' + '.join(
        #             f"({d}{i} - {d}{j + len(self._anchor_points)})**2"
        #             for d in ['x', 'y', 'q']
        #         ) + f" == {self._generate_distance(pa,pl):.3f}"
        #         for i,pa in enumerate(self._anchor_points)
        #     )
        #     for j,pl in enumerate(self._ligand_points)
        # ))
        points_distances = []
        for j, pl in enumerate(self._ligand_points):
            points_distances.append([])
            for i, pa in enumerate(self._anchor_points):
                points_distances[j].append(self._coord_distance(pa,pl))
        builder.update(full_points_constraints=',\n'.join(
            ',\n'.join(
                ', '.join(
                    f"({d}{g} - {d}{i + len(self._anchor_points)}) == {dist[g][j]:.3f}"
                    for j, d in enumerate(['x', 'y', 'q'])
                )
            for g, pa in enumerate(self._anchor_points)
            )
        for i,dist in enumerate(points_distances)
        ) + ','
        )

        builder.update(rotation_constraints=',\n'.join(
            # f"xx{i} == x{i} - y{i}*r2 + q{i}*r1, yy{i} == x{i}*r2 + y{i} - q{i}*r0, qq{i} == -x{i}*r1 + y{i}*r0 + q{i}"
            f"xx{i} == x{i}*cosb*cosg + y{i}*(sina*sinb*cosg - cosa*sing) + q{i}*(cosa*sinb*cosg + sina*sing), "+ 
            f"yy{i} == x{i}*cosb*sing + y{i}*(sina*sinb*sing + cosa*cosg) + q{i}*(cosa*sinb*sing - sina*cosg), "+ 
            f"qq{i} == x{i}*(-sinb) + y{i}*sina*cosb + q{i}*cosa*cosb"
        for i in range(len(self._anchor_points) + len(self._ligand_points))
        ) + ','
        )

        builder.update(energy=' + '.join(
            f'apply_interpolation(xx{i}, yy{i}, qq{i}, maps[atomst[{i}]], "z3")'
        for i in range(len(self._anchor_points) + len(self._ligand_points))
        ) + ','
        )

        # (float(model.evaluate(x1).as_decimal(10).strip('?')), 
        #  float(model.evaluate(y1).as_decimal(10).strip('?')), 
        #  float(model.evaluate(q1).as_decimal(10).strip('?'))),
        builder.update(returned_points=',\n'.join('('+
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

class Builder:
    LEFT_DELIMITER = '{{{{'
    RIGHT_DELIMITER = '}}}}'
    MAGIC_STRING_1 = 'zqwezrtyzuiozpaszdfgzhjkzlzxzcvbznm1z234z567z890z'
    MAGIC_STRING_2 = 'zmnbzvcxzzlkzjhgzfdszapoziuyztrezwq0z987z654z321z'

    def __init__(self, string: str) -> 'Builder':
        self._string: str = string
        self._kwargs: Dict[str, str] = dict()

    @ classmethod
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