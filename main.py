import sys, time, subprocess, ast, os
from template_manager.template_manager import IntSolver
import util.create_pdb as pdb
import util.ligand as lig
import util.tri_linear_interpolation as inter
import util.general_util as utl
import numpy as np
import matplotlib.pyplot as plt

def main():
    protein = sys.argv[1]
    ligand = sys.argv[2]
    res = int(sys.argv[3]) if len(sys.argv) > 3 else 4.0
    atom_contr = ''

    with open(f'protein_ligand/{protein}_{ligand}/input/{protein}_h.maps.fld') as f:
        lines = f.readlines()
    if not os.path.isdir(f'protein_ligand/{protein}_{ligand}/smtdock'):
        os.makedirs(f'protein_ligand/{protein}_{ligand}/smtdock')
    ligand_path = f'protein_ligand/{protein}_{ligand}/input/{ligand}.pdbqt'
    ligand_points, atoms, atomst = lig.get_ligand_points(ligand_path)
    ligand_anchors, ligand_points = lig.separate_anchors_points(ligand_points)
    print(f'Ligand anchors: \n{ligand_anchors}')
    
    for line in lines:
        # print(line)
        if '#NELEMENTS' in line:
            size = [int(x) for x in line.split(' ')[1:4]]
        if '#SPACING' in line:
            spacing = float(line.split(' ')[1])
        if '#CENTER' in line:
            center = [float(x) for x in line.split(' ')[1:4]]
    print(f'Center: {center}')
    soa_points = lig.get_soa_ligand(f'protein_ligand/{protein}_{ligand}/results/{ligand}_{protein}_out01.pdb', center)
    energy_soa = 0
    for i,p in enumerate(soa_points):
        values = []
        with open(f'protein_ligand/{protein}_{ligand}/input/{protein}_h.{atomst[i]}.map') as m:
            for line in m:
                try:
                    values.append(float(line.strip()))
                except ValueError:
                    continue
        ipol = inter.apply_interpolation(p[0], p[1], p[2], values, 'notz3', spacing*size[0]/2, 1/spacing, size[0]+1)
        print(atomst[i], ipol)
        atom_contr += f'{atomst[i]} {ipol}\n'
        energy_soa += ipol
    print(f'energy soa: {energy_soa}')
    atom_contr += f'energy soa: {energy_soa}\n\n'
    pdb.center_protein(center, protein, ligand)

    space_bounds = int((size[0]/2) * res)
    solver = IntSolver(protein, ligand, space_bounds, size[0], spacing, res, ligand_anchors, ligand_points, atomst)

    solver.generate_script()
    start = time.time()
    process = subprocess.run(['python3', 'solver_script.py'], stdout=subprocess.PIPE)
    end = time.time()
    # lines = process.stdout.strip().split("\n")
    lines = process.stdout.decode().split("\n")
    final_pose = ast.literal_eval(lines[1])
    coms = ast.literal_eval(lines[0])
    print(coms)
    
    mult = spacing/res
    
    final_pose = [(x * mult, y * mult, z * mult) for x, y, z in final_pose]
    print(f'Final pose of ligand anchors returned from solver: {final_pose}\nTime elapsed: {end - start}')
    # transformed_points = rigid_transform(ligand_points, ligand_anchors[0], ligand_anchors[1], ligand_anchors[2], final_pose[0], final_pose[1], final_pose[2])
    pdb.create_pdb(atoms, final_pose, center, protein, ligand)
    energy = 0
    for i,p in enumerate(final_pose):
        values = []
        with open(f'protein_ligand/{protein}_{ligand}/input/{protein}_h.{atomst[i]}.map') as m:
            for line in m:
                try:
                    values.append(float(line.strip()))
                except ValueError:
                    continue
        ipol = inter.apply_interpolation(p[0], p[1], p[2], values, "notz3", spacing*size[0]/2, 1/spacing, size[0]+1)
        print(atomst[i], ipol)
        atom_contr += f'{atomst[i]} {ipol}\n'
        energy += ipol
    print(f'energy sat: {energy}')
    atom_contr += f'energy sat: {energy}\n'
    with open(f'protein_ligand/{protein}_{ligand}/smtdock/atom_contribution.txt', 'w') as f:
        f.write(atom_contr)
    with open(f'protein_ligand/{protein}_{ligand}/smtdock/energy.csv', 'w') as f:
        f.write(f'energy: {energy}\ntime: {end-start}')
    return final_pose

main()