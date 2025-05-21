from Bio.PDB import PDBParser, PDBIO

def mod_pdb_coords(structure, atoms, new_points):
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() in atoms:
                        atom.coord = new_points[atoms.index(atom.get_name())]

def translate_protein(structure, center):
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.coord = [atom.coord[i] - center[i] for i in range(3)]

def shift_ligand(structure, center):
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.coord = [atom.coord[i] + center[i] for i in range(3)]

def center_protein(center):
    pdb_file = "a2a_adenosine/a2a.pdb"
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure("original_pdb1", pdb_file)

    translate_protein(structure, center)

    output_pdb = f'centered_a2a.pdb'
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(output_pdb , write_end = False , preserve_atom_numbering = False)

    return output_pdb

def create_pdb(atoms, new_points, typ, center):
    pdb_file = "a2a_adenosine/input/adenosine.pdbqt"
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure("original_pdb1", pdb_file)

    mod_pdb_coords(structure, atoms, new_points)
    shift_ligand(structure, center)

    output_pdb = f'new_adenosine_{typ}.pdb'
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(output_pdb , write_end = False , preserve_atom_numbering = False)

    return output_pdb