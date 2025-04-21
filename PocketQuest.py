from utils.classes import *
import utils.check_files as cf 
import utils.feature_extraction as fe 
import utils.clustering as cl
import utils.make_pdb as pdb
import utils.visualization as vis
import sys
import os

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python PocketQuest.py <protein.pdb> [--vis N] [--rotate]")
        sys.exit(1)

    pdb_dir = sys.argv[1]
    pdb_base_name = os.path.basename(pdb_dir)
    pdb_name = os.path.splitext(pdb_base_name)[0]

    vis_n = None
    rotate = False
    if "--vis" in sys.argv:
        vis_index = sys.argv.index("--vis")
        if vis_index + 1 < len(sys.argv):
            vis_n = int(sys.argv[vis_index + 1])

    if "--rotate" in sys.argv:
        rotate = True

    if not cf.file_exists(pdb_dir):
        sys.exit(1)

    if not cf.can_read_file(pdb_dir):
        sys.exit(1)

    protein = Protein(pdb_dir, pdb_name)
    feature_vector = fe.generate_vector(protein)
    clusters = cl.cluster_points(protein, feature_vector)
    pdb.make_pdb(protein, clusters)

    if vis_n is not None:
        vis.visualize_clusters(pdb_name, top_n=vis_n, rotate=rotate)
