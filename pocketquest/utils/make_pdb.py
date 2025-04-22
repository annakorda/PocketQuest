from pocketquest.utils.classes import *
import os

def make_pdb(protein, clusters, results_dir="results"):
    """
    Writes out separate PDB files for each predicted cluster.

    Parameters:
    - protein: Protein object with original .pdb file path and name
    - clusters: List of (label, cluster) tuples
    - results_dir: Directory where output .pdb files will be written
    """
    print("Starting writing PDB")

    os.makedirs(results_dir, exist_ok=True)

    for label, cluster in clusters:
        pdb_path = protein.get_dir()
        residues = cluster.get_residues()
        output = ""

        with open(pdb_path, "r") as pdb_file:
            for line in pdb_file:
                if not line.startswith("ATOM"):
                    continue
                residue = line[17:20].strip()
                residue_number = line[22:26].strip()
                residue_tuple = (residue, residue_number)
                if residue_tuple in residues:
                    output += line

        filename = f"{protein.get_name()}cluster_{label}.pdb"
        output_file = os.path.join(results_dir, filename)
        with open(output_file, "w") as file:
            file.write(output)

    print("Finished writing PDB")
