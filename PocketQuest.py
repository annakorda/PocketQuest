from utils.classes import *
import utils.check_files as cf 
import utils.feature_extraction as fe 
import utils.clustering as cl
import utils.make_pdb as pdb
# import numpy as np
import sys
import os
# from pathlib import Path


if __name__ == "__main__":
    # temp measure
    if len(sys.argv) != 2:
        print("Usage: python PocketQuest.py <protein.pdb> ")
        sys.exit(1)
    
    pdb_dir = sys.argv[1]
    pdb_base_name = os.path.basename(pdb_dir)
    pdb_name = os.path.splitext(pdb_base_name)[0]

    if not cf.file_exists(pdb_dir):
        sys.exit(1)

    if not cf.can_read_file(pdb_dir):
        sys.exit(1)

    protein = Protein(pdb_dir, pdb_name)
    feature_vector = fe.generate_vector(protein)
    clusters = cl.cluster_points(protein, feature_vector)
    pdb.make_pdb(protein, clusters)