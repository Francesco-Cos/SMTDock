from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
import numpy as np

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

def center_protein(center, protein, ligand):
    pdb_file = f'protein_ligand/{protein}_{ligand}/input/{protein}_h.pdb'
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure("original_pdb1", pdb_file)

    translate_protein(structure, center)

    output_pdb = f'protein_ligand/{protein}_{ligand}/smtdock/centered_{protein}.pdb'
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(output_pdb , write_end = False , preserve_atom_numbering = False)

    return output_pdb

def create_pdb(atoms, new_points, center, protein, ligand):
    # Initialize structure
    structure = Structure("MyStructure")
    model = Model(0)
    chain = Chain("A")
    residue = Residue((" ", 1, " "), "RES", " ")  # Generic residue
    numerate = True

    # Count atoms of each element type for unique naming
    element_counts = {}

    if any(len(atom)>1 for atom in atoms):
        numerate = False
    # Add atoms to residue
    for i, (coord, element) in enumerate(zip(new_points, atoms), start=1):
        # Increment element-specific counter
        count = element_counts.get(element, 0) + 1
        element_counts[element] = count
        if numerate:
            atom_name = f"{element}{count}"
        else:
            atom_name = f"{element}"
        padded_name = f"{atom_name:<4}"  # PDB requires 4-char names

        atom = Atom(
            name=atom_name,
            coord=np.array(coord),
            bfactor=0.0,
            occupancy=1.0,
            altloc=" ",
            fullname=padded_name,
            serial_number=i,
            element=element[0]
        )
        residue.add(atom)

    # Build the structure hierarchy
    chain.add(residue)
    model.add(chain)
    structure.add(model)
    shift_ligand(structure, center)

    # Save to file
    io = PDBIO()
    io.set_structure(structure)
    io.save(f'protein_ligand/{protein}_{ligand}/smtdock/new_{ligand}.pdb')