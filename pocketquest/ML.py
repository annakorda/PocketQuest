from utils.classes import *
import utils.check_files as cf 
import utils.feature_extraction as fe 
import numpy as np
import sys
import os
from pathlib import Path

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python ML.py <file_with_pdb_names.txt> <output_file_name.npy>")
        sys.exit(1)
    
    pdb_list_file = sys.argv[1]
    output_name = sys.argv[2]  # just the name, not full path

    if not cf.file_exists(pdb_list_file):
        sys.exit(1)

    if not cf.can_read_file(pdb_list_file):
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

    all_vectors = []
    processed_count = 0
    for protein in proteins:
        vectors = fe.generate_vector(protein, ML=True)
        if vectors is None:
          continue
        all_vectors.extend(vectors)
        processed_count += 1

    # Create input_numpy directory in current folder if not exists
    output_dir = Path.cwd() / "input_numpy"
    output_dir.mkdir(exist_ok=True)

    # Save final .npy file into that directory
    output_path = output_dir / output_name
    np.save(output_path, np.array(all_vectors))
    print(f"Saved {len(all_vectors)} feature vectors to: {output_path}")

    print(f"\nSuccessfully processed {processed_count} out of {len(proteins)} PDBs.")


