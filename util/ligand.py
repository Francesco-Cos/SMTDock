import re
import numpy as np

def select_ligand_anchors(ligand):
    with open(ligand) as f:
        data = f.read()
    points = []
    branches = re.split('ROOT|ENDROOT|BRANCH|ENDBRANCH', data)
    branches = [b.split('\n') for b in branches[1:]]
    branches = [sublist for sublist in [[item for item in sublist if 'ATOM' in item] for sublist in branches] if sublist]
    for i in range(3):
        atom = branches[i][0].split()
        points.append([float(atom[5]), float(atom[6]), float(atom[7])])
    return np.array(points)

def get_ligand_points(ligand):
    with open(ligand) as f:
        data = f.read()
    anchors = []
    points = []
    atoms = []
    atomst = []
    branches = re.split('ROOT|ENDROOT|BRANCH|ENDBRANCH', data)
    branches = [b.split('\n') for b in branches[1:]]
    branches = [sublist for sublist in [[item for item in sublist if 'ATOM' in item] for sublist in branches] if sublist]
    i = 0
    for branch in branches:
        for j,a in enumerate(branch):
            atom = a.split()
            if j == 0 and i < 3:
                anchors.append([float(atom[5]), float(atom[6]), float(atom[7])])
                i += 1
            else:
                points.append([float(atom[5]), float(atom[6]), float(atom[7])])
            atoms.append(atom[2])
            atomst.append(atom[11])
    return np.array(anchors), np.array(points), atoms, atomst

def get_soa_ligand(ligand,center):
    with open(ligand) as f:
        data = f.readlines()
    points = []
    for line in data:
        if 'HETATM' in line:
            atom = line.split()
            points.append([float(atom[5]) - center[0], float(atom[6]) - center[1], float(atom[7]) - center[2]])
    return np.array(points)