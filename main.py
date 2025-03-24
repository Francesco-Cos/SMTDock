import sys, math, time, subprocess, ast
from template_manager.template_manager import RealSolver
from solver_real import run_solver_real
from solver_int import run_solver_int
from util.ligand import select_ligand_anchors, get_ligand_points, get_soa_ligand
from util.create_pdb import create_pdb, center_protein
from util.tri_linear_interpolation import apply_interpolation

def ceildiv(a, b):
    return -(a // -b)

def main():
    typing = sys.argv[1]
    with open('a2a_adenosine/input/a2a_h.gpf') as f:
        lines = f.readlines()
    ligand_path = 'a2a_adenosine/input/adenosine.pdbqt'
    # ligand_anchors = select_ligand_anchors(ligand_path)
    ligand_anchors, ligand_points, atoms, atomst = get_ligand_points(ligand_path)

    print(f'Ligand anchors: \n{ligand_anchors}')
    d12 = [ligand_anchors[0][i] - ligand_anchors[1][i] for i in range(3)]
    d13 = [ligand_anchors[0][i] - ligand_anchors[2][i] for i in range(3)]
    d23 = [ligand_anchors[1][i] - ligand_anchors[2][i] for i in range(3)]
    points, center = [], []
    spacing = 0
    
    for line in lines:
        if 'npts' in line:
            points = [int(x) for x in line.split(' ')[1:4]]
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
        ipol = apply_interpolation(p[0], p[1], p[2], values, 'notz3')
        energy_soa += ipol
    print(f'energy soa: {energy_soa}')
    center_protein(center)
    center = [0,0,0]

    if typing == 'real':
        # Real version of the solver
        print('Real solver')
        d12 = float(f'{math.sqrt(sum(x**2 for x in d12)):.3f}')
        d13 = float(f'{math.sqrt(sum(x**2 for x in d13)):.3f}')
        d23 = float(f'{math.sqrt(sum(x**2 for x in d23)):.3f}')

        spac = [p*spacing for p in points]
        # print(ligand_points)
        rs = RealSolver(center, spac, ligand_anchors, ligand_points, atomst)
        
        rs.generate_script()
        start = time.time()
        process = subprocess.run(['python3', 'solver_script.py'], stdout=subprocess.PIPE)
        end = time.time()
        final_pose = ast.literal_eval(process.stdout.decode())
        # print(final_pose)
        energy = 0
        for i,p in enumerate(final_pose):
            values = []
            with open(f'a2a_adenosine/input/a2a_h.{atomst[i]}.map') as m:
                for line in m:
                    try:
                        values.append(float(line.strip()))
                    except ValueError:
                        continue
            ipol = apply_interpolation(p[0], p[1], p[2], values, "notz3")
            energy += ipol
        print(f'energy sat: {energy}')
        print(f'Final pose of ligand anchors returned from solver: {final_pose}\nTime elapsed: {end - start}')
        # transformed_points = rigid_transform(ligand_points, ligand_anchors[0], ligand_anchors[1], ligand_anchors[2], final_pose[0], final_pose[1], final_pose[2])
        create_pdb(atoms, final_pose, typing)
    else:
        # Integer version of the solver (lattice)
        print('Integer solver')
        d12 = [ceildiv(el, spacing) for el in d12]
        d13 = [ceildiv(el, spacing) for el in d13]

        d12 = float(f'{math.sqrt(sum(x**2 for x in d12)/spacing):.3f}')
        d13 = float(f'{math.sqrt(sum(x**2 for x in d13)/spacing):.3f}')
        d23 = float(f'{math.sqrt(sum(x**2 for x in d23)/spacing):.3f}')
        start = time.time()
        final_points = run_solver_int(points,d12,d13,d23)
        final_points = [(x - points[0]/2, y - points[1]/2, z - points[2]/2) for x, y, z in final_points]
        print(f'Final lattice points: {final_points}')
        final_pose = [(center[0] + final_point[0]*spacing, center[1] + final_point[1]*spacing, center[2] + final_point[2]*spacing) for final_point in final_points]
        end = time.time()
        print(f'Final pose of ligand anchors returned from solver: {final_pose}\nTime elapsed: {end - start}')
    return final_pose

main()