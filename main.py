import sys, time, subprocess, ast
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

    with open(f'{protein}_{ligand}/input/{protein}_h.gpf') as f:
        lines = f.readlines()
    ligand_path = f'{protein}_{ligand}/input/{ligand}.pdbqt'
    ligand_points, atoms, atomst = lig.get_ligand_points(ligand_path)
    ligand_anchors, ligand_points = lig.separate_anchors_points(ligand_points)
    print(f'Ligand anchors: \n{ligand_anchors}')
    
    for line in lines:
        if 'npts' in line:
            size = [int(x) for x in line.split(' ')[1:4]]
        if 'spacing' in line:
            spacing = float(line.split(' ')[1])
        if 'gridcenter' in line:
            center = [float(x) for x in line.split(' ')[1:4]]
    print(f'Center: {center}')

    soa_points = lig.get_soa_ligand(f'{protein}_{ligand}/results/{ligand}_{protein}_out01.pdb', center)
    energy_soa = 0
    for i,p in enumerate(soa_points):
        values = []
        with open(f'{protein}_{ligand}/input/{protein}_h.{atomst[i]}.map') as m:
            for line in m:
                try:
                    values.append(float(line.strip()))
                except ValueError:
                    continue
        ipol = inter.apply_interpolation(p[0], p[1], p[2], values, 'notz3', spacing*size[0]/2, 1/spacing)
        print(atomst[i], ipol)
        energy_soa += ipol
    print(f'energy soa: {energy_soa}')
    pdb.center_protein(center, protein, ligand)

    space_bounds = int((size[0]/2) * res)
    solver = IntSolver(space_bounds, size[0], 21, spacing, res, ligand_anchors, ligand_points, atomst)

    solver.generate_script()
    start = time.time()
    process = subprocess.run(['python3', 'solver_script.py'], stdout=subprocess.PIPE)
    end = time.time()
    # lines = process.stdout.strip().split("\n")
    lines = process.stdout.decode().split("\n")
    final_pose = ast.literal_eval(lines[1])
    coms = ast.literal_eval(lines[0])
    print(coms)
    
    # final_pose = ast.literal_eval(process.stdout.decode())

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot with color mapping
    sc = ax.scatter3D([point[0] for point in coms], [point[1] for point in coms], [point[2] for point in coms])
    plt.colorbar(sc, ax=ax, label='Energy values')

    # Labels
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    ax.set_title('3D Scatter Plot with Color Mapping')

    plot_name = f'centers_of_mass.pdf'
    plt.savefig(plot_name)
    plt.show()
    mult = spacing/res
    
    final_pose = [(x * mult, y * mult, z * mult) for x, y, z in final_pose]
    print(f'Final pose of ligand anchors returned from solver: {final_pose}\nTime elapsed: {end - start}')
    # transformed_points = rigid_transform(ligand_points, ligand_anchors[0], ligand_anchors[1], ligand_anchors[2], final_pose[0], final_pose[1], final_pose[2])
    pdb.create_pdb(atoms, final_pose, center, protein, ligand)
    energy = 0
    for i,p in enumerate(final_pose):
        values = []
        with open(f'{protein}_{ligand}/input/{protein}_h.{atomst[i]}.map') as m:
            for line in m:
                try:
                    values.append(float(line.strip()))
                except ValueError:
                    continue
        ipol = inter.apply_interpolation(p[0], p[1], p[2], values, "notz3", spacing*size[0]/2, 1/spacing)
        print(atomst[i], ipol)
        energy += ipol
    print(f'energy sat: {energy}')
    return final_pose

main()