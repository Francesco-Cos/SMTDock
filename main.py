import sys, math, time, subprocess, ast
from template_manager.template_manager import RealSolver, IntSolver
from util.ligand import select_ligand_anchors, get_ligand_points, get_soa_ligand
from util.create_pdb import create_pdb, center_protein
from util.tri_linear_interpolation import apply_interpolation
import numpy as np

def generate_distance(p1, p2):
        diff = [e1 - e2 for e1, e2 in zip(p1, p2)]
        return sum([d**2 for d in diff])

def center_of_mass(points):
    com = [0,0,0]
    for point in points:
        for i, coords in enumerate(point):
            com[i] += coords
    com = [el/len(points) for el in com]
    return com

def main():
    typing = sys.argv[1]
    with open('a2a_adenosine/input/a2a_h.gpf') as f:
        lines = f.readlines()
    ligand_path = 'a2a_adenosine/input/adenosine.pdbqt'
    # ligand_anchors = select_ligand_anchors(ligand_path)
    ligand_anchors, ligand_points, atoms, atomst = get_ligand_points(ligand_path)
    com = center_of_mass(ligand_points)
    distances_to_com = [generate_distance(p1,com) for p1 in ligand_points]
    pidx = [i for i in range(len(distances_to_com))]
    sorted_d, sorted_id = zip(*sorted(zip(distances_to_com, pidx)))
    remove_ids = [sorted_id[i] for i in range(3)]
    ligand_anchors = np.array([ligand_points[sorted_id[i]] for i in range(3)])
    ligand_points = [point for i,point in enumerate(ligand_points) if i not in remove_ids]
    print(f'Ligand anchors: \n{ligand_anchors}')
    points, center = [], []
    spacing = 0
    res = 4.0
    
    for line in lines:
        if 'npts' in line:
            size = [int(x) for x in line.split(' ')[1:4]]
        if 'spacing' in line:
            spacing = float(line.split(' ')[1])
        if 'gridcenter' in line:
            center = [float(x) for x in line.split(' ')[1:4]]
    print(f'Center: {center}')
    points = [10,10,10]
    soa_points = get_soa_ligand('a2a_adenosine/results/adenosine_a2a_out01.pdb', center)
    energy_soa = 0
    for i,p in enumerate(soa_points):
        values = []
        with open(f'a2a_adenosine/input/a2a_h.{atomst[i]}.map') as m:
            for line in m:
                try:
                    values.append(float(line.strip()))
                except ValueError:
                    continue
        ipol = apply_interpolation(p[0], p[1], p[2], values, 'notz3', res, spacing*size[0]/2, 1/spacing)
        energy_soa += ipol
    print(f'energy soa: {energy_soa}')
    center_protein(center)
    new_center = [0,0,0]

    bounds = [int(p*res) for p in points] if typing=='int' else [p*spacing for p in points]
    space_bounds = int((size[0]/2) * res)
    # print(ligand_points)
    solver = IntSolver(new_center, bounds, space_bounds, size[0], 20, spacing, int(res), ligand_anchors, ligand_points, atomst) if typing=='int' else RealSolver(center, bounds, res, ligand_anchors, ligand_points, atomst)
    
    solver.generate_script()
    start = time.time()
    process = subprocess.run(['python3', 'solver_script.py'], stdout=subprocess.PIPE)
    end = time.time()
    final_pose = ast.literal_eval(process.stdout.decode())
    mult = spacing/res
    if typing=='int':
        final_pose = [(x * mult, y * mult, z * mult) for x, y, z in final_pose]
    print(f'Final pose of ligand anchors returned from solver: {final_pose}\nTime elapsed: {end - start}')
    # transformed_points = rigid_transform(ligand_points, ligand_anchors[0], ligand_anchors[1], ligand_anchors[2], final_pose[0], final_pose[1], final_pose[2])
    create_pdb(atoms, final_pose, typing, center)
    energy = 0
    for i,p in enumerate(final_pose):
        values = []
        with open(f'a2a_adenosine/input/a2a_h.{atomst[i]}.map') as m:
            for line in m:
                try:
                    values.append(float(line.strip()))
                except ValueError:
                    continue
        ipol = apply_interpolation(p[0], p[1], p[2], values, "notz3", res, spacing*size[0]/2, 1/spacing)
        print(atomst[i], ipol)
        energy += ipol
    print(f'energy sat: {energy}')
    return final_pose

main()