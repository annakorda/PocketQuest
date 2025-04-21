from utils.classes import *
import utils.check_files as cf 
import utils.feature_extraction as fe 
import utils.clustering as cl
import utils.make_pdb as pdb
import utils.visualization as vis
import argparse
import sys
import os

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-in", dest = "pdb", required=True, help = "PDB file to check")
    parser.add_argument("-distance", type = cf.positive_float, help = "Distance cutoff while performing the clustering. Positive float.")
    parser.add_argument("-size", type = cf.positive_int, help = "Minimum number of points in a valid cluster. Positive integer.")
    parser.add_argument("-max_residues", type = cf.positive_int, help = "Maximum number of residues in a valid cluster. Positive integer.")
    parser.add_argument("-prob", type = cf.probability_float, help = "Threshold probability to consider a point to be part of a binding site. Float between 0 and 1.")
    parser.add_argument("-vis", type= cf.positive_int, help = "Number of clusters to be visualized. Positive integer.")
    # FINISH HELP
    parser.add_argument("-rotate", action = "store_true", help="")

    args = parser.parse_args()
    print(args)
    pdb_dir = args.pdb
    pdb_base_name = os.path.basename(pdb_dir)
    pdb_name = os.path.splitext(pdb_base_name)[0]

    if not cf.file_exists(pdb_dir):
        sys.exit(1)

    if not cf.can_read_file(pdb_dir):
        sys.exit(1)

    # vis_n = None
    # rotate = False
    # if "--vis" in sys.argv:
    #     vis_index = sys.argv.index("--vis")
    #     if vis_index + 1 < len(sys.argv):
    #         vis_n = int(sys.argv[vis_index + 1])

    # if "--rotate" in sys.argv:
    #     rotate = True

    protein = Protein(pdb_dir, pdb_name)
    feature_vector = fe.generate_vector(protein)
    clusters = cl.cluster_points(protein, feature_vector, args)
    pdb.make_pdb(protein, clusters)

    if args.vis is not None:
        vis.visualize_clusters(pdb_name, top_n=args.vis, rotate=args.rotate)
