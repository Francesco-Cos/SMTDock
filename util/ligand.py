import re
import numpy as np
import util.general_util as utl

def separate_anchors_points(ligand):
    com = utl.center_of_mass(ligand)
    distances_to_com = [utl.generate_distance(p1,com) for p1 in ligand]
    pidx = [i for i in range(len(distances_to_com))]
    _ , sorted_id = zip(*sorted(zip(distances_to_com, pidx)))
    remove_ids = [sorted_id[i] for i in range(3)]
    ligand_anchors = np.array([ligand[sorted_id[i]] for i in range(3)])
    ligand_points = [point for i,point in enumerate(ligand) if i not in remove_ids]
    return ligand_anchors, ligand_points

def get_ligand_points(ligand):
    with open(ligand) as f:
        data = f.read()
    points = []
    atoms = []
    atomst = []
    branches = re.split('ROOT|ENDROOT|BRANCH|ENDBRANCH', data)
    branches = [b.split('\n') for b in branches[1:]]
    branches = [sublist for sublist in [[item for item in sublist if ('ATOM' in item or 'HETATM' in item)] for sublist in branches] if sublist]
    i = 0
    for branch in branches:
        for j,a in enumerate(branch):
            atom = a.split()
            if len(atom) > 12:
                atom.pop(4)
            points.append([float(atom[5]), float(atom[6]), float(atom[7])])
            atoms.append(atom[2])
            atomst.append(atom[11])
    return np.array(points), atoms, atomst

def get_soa_ligand(ligand,center):
    with open(ligand) as f:
        data = f.readlines()
    points = []
    for line in data:
        if 'HETATM' in line:
            atom = line.split()
            if len(atom) > 11:
                atom.pop(4)
            points.append([float(atom[5]) - center[0], float(atom[6]) - center[1], float(atom[7]) - center[2]])
    return np.array(points)