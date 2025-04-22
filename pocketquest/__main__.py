from pocketquest.utils.classes import *
import pocketquest.utils.check_files as cf
import pocketquest.utils.feature_extraction as fe
import pocketquest.utils.clustering as cl
import pocketquest.utils.make_pdb as pdb
import pocketquest.utils.visualization as vis
import argparse
import sys
import os
import textwrap
import shutil

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            PocketQuest: Predict and visualize ligand binding sites from protein PDB files.

            Example:
              pocketquest -i /path/to/2brd.pdb -distance 2.0 -size 10 -max_residues 40 -prob 0.8 -vis 3 --rotate
        """)
    )

    parser.add_argument("-i", "--input", dest="pdb", required=True, help="Path to input PDB file.")
    parser.add_argument("-distance", type=cf.positive_float, help="Distance cutoff for clustering (positive float).")
    parser.add_argument("-size", type=cf.positive_int, help="Minimum number of points in a cluster (positive int).")
    parser.add_argument("-max_residues", type=cf.positive_int, help="Maximum number of residues in a cluster (positive int).")
    parser.add_argument("-prob", type=cf.probability_float, help="Minimum prediction probability (0–1).")
    parser.add_argument("-vis", type=cf.positive_int, help="Number of clusters to visualize (positive int).")
    parser.add_argument("-rotate", action="store_true", help="Rotate structure for better Chimera orientation.")

    args = parser.parse_args()

    #  Capture the user's working directory to save results there
    user_run_dir = os.getcwd()

    #  Absolute path to the input PDB file
    input_pdb_path = os.path.abspath(args.pdb)

    if not os.path.isfile(input_pdb_path):
        print(f"Error: File '{args.pdb}' does not exist.")
        sys.exit(1)

    #  Extract PDB name from the filename
    pdb_base_name = os.path.basename(input_pdb_path)
    pdb_name = os.path.splitext(pdb_base_name)[0]

    #  Set repo root so the code works regardless of run location
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(repo_root)

    #  Save results in the user’s run directory
    results_dir = os.path.join(user_run_dir, "results")
    os.makedirs(results_dir, exist_ok=True)

    #  Copy the input PDB to the results directory
    copied_pdb_path = os.path.join(results_dir, f"{pdb_name}.pdb")
    shutil.copy(input_pdb_path, copied_pdb_path)
    print(f"Using PDB: {copied_pdb_path}")

    if not cf.file_exists(copied_pdb_path) or not cf.can_read_file(copied_pdb_path):
        print(f"Error: Could not read copied PDB at {copied_pdb_path}")
        sys.exit(1)

    #  Run PocketQuest core pipeline
    protein = Protein(copied_pdb_path, pdb_name)
    feature_vector = fe.generate_vector(protein)
    clusters = cl.cluster_points(protein, feature_vector, args, results_dir=results_dir)


    pdb.make_pdb(protein, clusters, results_dir=results_dir)

    if args.vis is not None:
        vis.visualize_clusters(pdb_name, top_n=args.vis, rotate=args.rotate, results_dir=results_dir)

if __name__ == "__main__":
    main()
