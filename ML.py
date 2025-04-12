from utils.classes import *
import utils.check_files as cf 
import utils.feature_extraction as fe 
import sys
import os

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python ML.py <file_with_pdb_names.txt>")
        sys.exit(1)
    
    pdb_list_file = sys.argv[1]

    if not cf.file_exists(pdb_list_file):
        sys.exit(1)

    if not cf.can_read_file:
        sys.exit(1)

    proteins = []
    with open(pdb_list_file, "r") as file:
        for line in file:
            pdb_dir = line.strip()
            if(len(pdb_dir) == 0):
                continue
            pdb_base_name = os.path.basename(pdb_dir)
            pdb_name = os.path.splitext(pdb_base_name)[0]

            if not cf.file_exists(pdb_dir) or not cf.can_read_file(pdb_dir):
                print(f"{pdb_dir} will be skipped.")
            else:
                proteins.append(Protein(pdb_dir, pdb_name))
    file.close()

    for protein in proteins:
        print(fe.generate_vector(protein))
        
        
        