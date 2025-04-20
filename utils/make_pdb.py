from utils.classes import *
import os

def make_pdb(clusters):
    print("Starting writing PDB")
    for cluster in clusters:
        filename = "cluster_" + str(cluster[0]) + ".pdb"
        output_dir = "results"
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, filename)

        with open(output_file, "w") as file:
            num = 1
            for atom_tuple in cluster[1].get_atoms():
                atom = atom_tuple[0]
                x, y, z = atom.get_coord()
                atom_name = atom.get_type()
                residue_name = atom.get_residue() 
                element = atom.get_element().strip()
                bfactor = atom.get_bfactor()
                pdb_line = (
                    f"ATOM  {num:5d} {atom_name:<4} {residue_name:<3} A{1:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{bfactor:6.2f}{element:>2}\n"
                )
                file.write(pdb_line)
                num += 1
        file.close()
    print("Finished writing PDB")