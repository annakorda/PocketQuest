from utils.classes import *
import utils.check_files as cf
import subprocess
import os
import glob
import sys
import numpy as np
from scipy.spatial import cKDTree

def read_pbs(protein):
    print(f"Processing binding sites of {protein.get_name()}.pdb")
    protein_dir = protein.get_dir().replace(".pdb", "")
    name = protein_dir + "_bs_*.pdb"
    files = glob.glob(name)
    atoms = set()
    for file in files:
        with open(file, "r") as open_file:
            for line in open_file:
                if line[0:4] != "ATOM":
                    pass
                else:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coordinates = (x, y, z)
                    atoms.add(coordinates)
        open_file.close()
    print(f"Finished processing binding sites of {protein.get_name()}.pdb")
    return atoms

def read_pdb(protein, atom_list = set()):
    protein_dir = protein.get_dir()
    name = protein.get_name() + ".pdb"
    print(f"Reading {name}")
    with open(protein_dir, "r") as file:
        for line in file:
            if line[0:4] != "ATOM":
                pass
            else:
                number = line[6:11].strip()
                type = line[12:16].strip()
                residue = line[17:20].strip()
                residue_number = line[22:26].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coordinates = (x, y, z)
                bfactor = float(line[60:66])
                element = line[76:78].strip()
                new_atom = Atom(type, number, residue, residue_number, coordinates, bfactor, element)
                if coordinates in atom_list:
                    new_atom.set_in_pbs()
                protein.append_atom(new_atom)      
    file.close()
    print(f"Finished reading {name}")

def generate_xyzr(protein):
    print("Generating XYZR coordinates.")
    output_dir = "files"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, protein.get_name() + ".xyzr")
    pdb_file = protein.get_dir()
    script_dir = os.path.dirname(__file__)
    pdb_to_xyzr_path = os.path.join(script_dir, "..", "MSMS", "pdb_to_xyzr")
    try:
        with open(output_file, "w") as outfile:
            subprocess.run([pdb_to_xyzr_path, pdb_file], stdout=outfile, stderr=subprocess.DEVNULL, check=True)
        outfile.close()
        print("XYZR coordinates have been generated.")
    except subprocess.CalledProcessError:
        print("XYZR coordinates couldn't be generated")
        sys.exit(1) 

def generate_points(protein):
    print("Generating Connolly Points")
    output_dir = "files"
    os.makedirs(output_dir, exist_ok=True)
    xyzr_file = os.path.join(output_dir, protein.get_name() + ".xyzr")
    output_file = os.path.join(output_dir, protein.get_name())
    script_dir = os.path.dirname(__file__)
    msms_path = os.path.join(script_dir, "..", "MSMS", "MSMS")
    try:
        subprocess.run([msms_path, "-if", xyzr_file, "-of", output_file, "-probe_radius", str(1.6), "-no_header", "-density", str(0.5)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        print("Finished generating Connolly points.")
    except subprocess.CalledProcessError as e:
        print(f"Skipping {protein.get_name()} — MSMS error: {e}")
        return

def map_points_atoms(protein):
    output_dir = "files"
    os.makedirs(output_dir, exist_ok=True)
    points_file = os.path.join(output_dir, protein.get_name() + ".vert")
    atom_list = protein.get_atoms()
    atom_coord = np.array([atom.get_coord() for atom in atom_list])
    tree = cKDTree(atom_coord)
    try: 
        with open(points_file, "r") as file:
            for line in file:
                if line.startswith("#") or len(line.strip()) == 0:
                    continue
                x, y, z = map(float, line.strip().split()[:3])
                coord =(x, y, z)
                near_atoms = tree.query_ball_point(coord, r = 6)
                protrusion = len(tree.query_ball_point(coord, r = 10))
                point = Point(coord, protrusion)
                for i in near_atoms:
                    atom = atom_list[i]
                    point.append_atom(atom)
                protein.append_point(point)
        file.close()
    except Exception as e:
        print(f"Error mapping atoms to points: {e}")

def generate_features(protein, ML = False, freq = 0.3):
    print("Started generating feature vectors.")
    feature_vector = []

    for point in protein.get_points():
        point_vector = np.zeros(24)

        total_atoms = 0
        atom_density = 0
        num_carbon = 0
        num_oxigen = 0
        num_nitrogen = 0
        num_donors = 0
        num_acceptors = 0
        in_pbs = 0
        protrusion = point.get_protrusion()

        for atom in point.get_atoms(): 
            atom_obj, distance = atom[0], atom[1]
            total_atoms += 1

            residue = atom_obj.get_residue()
            atom_type = atom_obj.get_type()
            element = atom_obj.get_element()
            bfactor = atom_obj.get_bfactor()
            atom_in_pbs = atom_obj.get_in_pbs()
            scale = 1 - distance / 6

            if residue in residues:
                residue_vector = np.array(residues[residue][0])
                atom_dict = residues[residue][1]

                if atom_type in atom_dict:
                    atom_data = atom_dict[atom_type]
                    atom_vector = np.array(atom_data + [bfactor])
                    combined_vector = np.concatenate([residue_vector, atom_vector])
                    combined_vector = scale * combined_vector
                    point_vector += combined_vector

                    atom_density += scale
                    if element == "C":
                        num_carbon += 1
                    elif element == "O":
                        num_oxigen += 1
                    elif element == "N":
                        num_nitrogen += 1
                    donor = atom_data[5]
                    acceptor = atom_data[4]
                    num_donors += int(donor != 0)
                    num_acceptors += int(acceptor != 0)
                    if ML:
                        in_pbs += int(atom_in_pbs == True)
                # else:
                    # print(f"ATOM {atom_type} is not defined for residue {residue} and will be skipped.")
            # else:
                # print(f"Residue {residue} in ATOM number {atom[0].get_number()} is not defined and will be skipped.") 

        if total_atoms == 0:
            continue  # skip this point to avoid division by zero

        point_in_pbs = int((in_pbs / total_atoms) >= freq)
        if (ML == True):
             scalar_vector = np.array([total_atoms, atom_density, num_carbon, num_oxigen, num_nitrogen, num_donors, num_acceptors, protrusion, point_in_pbs])
        else:
             scalar_vector = np.array([total_atoms, atom_density, num_carbon, num_oxigen, num_nitrogen, num_donors, num_acceptors, protrusion])        
        final_vector = np.concatenate([point_vector, scalar_vector])
        feature_vector.append(final_vector)

    print("Finished generating feature vectors.")
    return feature_vector


def generate_vector(protein, ML = False):
    if (ML == True):
        atoms_in_pbs = read_pbs(protein)
        read_pdb(protein, atoms_in_pbs)
    else:
        read_pdb(protein)
    generate_xyzr(protein)
    generate_points(protein)
    map_points_atoms(protein)
    vector = generate_features(protein, ML)
    return vector