from utils.classes import *
import os

def make_pdb(protein, clusters):
    print("Starting writing PDB")
    for label, cluster in clusters:
        pdb_dir = protein.get_dir()
        residues = cluster.get_residues()
        output = ""
        with open(pdb_dir, "r") as pdb_file:
            for line in pdb_file:
                if line[0:4] != "ATOM":
                    pass
                else:
                    residue = line[17:20].strip()
                    residue_number = line[22:26].strip()
                    residue_tuple = (residue, residue_number)
                    if residue_tuple in residues:
                        output += line
        pdb_file.close()

        filename = f"{protein.get_name()}cluster{label}.pdb"
        output_dir = "results"
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, filename)
        with open(output_file, "w") as file:
            file.write(output)
        file.close()
        
    print("Finished writing PDB")