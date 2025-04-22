from utils.classes import *
import utils.check_files as cf 
import utils.feature_extraction as fe 
import utils.clustering as cl
import utils.make_pdb as pdb
import os

### evaluation.py ###
# In this file we are using 10 batches that were created for Active Learning that was not implemented in the end, in order to assess our results. These 
# batches will be available in the data drive folder (Assessment_Batches.tar.gz) and the batch_evaluation.log contains detailed info for the percentage of alignment of each predicted cluster
# with the real binding site per PDB (Only top 10 Clusters used). The class_distribution_topscoring.csv file contains the protein class distribution of the proteins 
# that were recovered more than 70% from top 10 clusters produced from PocketQuest.
#
# We are basically using PocketQuest (the exact same components) to evaluate the batches and additionally we define functions for metrics, but we exclude visualization for computational reasons. 


def evaluate_top_clusters(predicted_clusters, binding_residues):
    top_10_clusters = predicted_clusters[:10]
    hits_per_cluster = []
    for label, cluster in top_10_clusters:
        cluster_residues = set(cluster.get_residues())
        intersection = cluster_residues.intersection(binding_residues)
        percent_overlap = 100 * len(intersection) / len(binding_residues) if binding_residues else 0
        hits_per_cluster.append((label, percent_overlap))
    return hits_per_cluster

def get_binding_residues(bs_file):
    residues = set()
    with open(bs_file, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                resname = line[17:20].strip()
                resnum = line[22:26].strip()
                residues.add((resname, resnum))
    return residues

def run_evaluation_on_batch(batch_dir, f, global_stats):
    for filename in os.listdir(batch_dir):
        if filename.endswith(".pdb") and "_bs_" not in filename:
            global_stats['total_pdbs'] += 1

            pdb_path = os.path.join(batch_dir, filename)
            pdb_base = os.path.splitext(filename)[0]

            bs_files = [bs for bs in os.listdir(batch_dir) if bs.startswith(pdb_base + "_bs_")]
            if not bs_files or not cf.file_exists(pdb_path) or not cf.can_read_file(pdb_path):
                global_stats['skipped_pdbs'] += 1
                continue

            protein = Protein(pdb_path, pdb_base)
            feature_vector = fe.generate_vector(protein)
            clusters = cl.cluster_points(protein, feature_vector)

            if not clusters:
                global_stats['skipped_pdbs'] += 1
                f.write(f"\n{pdb_base} - Skipped: No clusters found after filtering.\n")
                continue

            pdb.make_pdb(protein, clusters)
            pdb_hit = False
            #Calculate a couple of metrics for the predicted clusters
            for bs_file in bs_files:
                global_stats['total_sites'] += 1
                bs_path = os.path.join(batch_dir, bs_file)
                binding_residues = get_binding_residues(bs_path)
                hits = evaluate_top_clusters(clusters, binding_residues)

                best_cluster, best_overlap = max(hits, key=lambda x: x[1]) if hits else (None, 0)
                global_stats['best_overlaps'][f"{pdb_base} vs {bs_file}"] = best_overlap

                if best_overlap >= 50:
                    global_stats['site_hit_count'] += 1
                    pdb_hit = True

                if best_overlap >= 70:
                    global_stats['overlap_distribution']["High (≥70%)"] += 1
                elif best_overlap >= 50:
                    global_stats['overlap_distribution']["Moderate (50–69%)"] += 1
                elif best_overlap >= 30:
                    global_stats['overlap_distribution']["Low (30–49%)"] += 1
                else:
                    global_stats['overlap_distribution']["Very Low (<30%)"] += 1

                f.write(f"\n{pdb_base} vs {bs_file}:\n")
                for label, score in hits:
                    f.write(f"  Cluster {label}: {score:.2f}% of binding site residues recovered\n")
                f.flush()

            if pdb_hit:
                global_stats['pdb_hit_count'] += 1

def summarize(f, stats):
    f.write("\n=== Evaluation Summary ===\n")
    f.write(f"Total PDBs processed: {stats['total_pdbs']}\n")
    f.write(f"PDBs skipped due to missing clusters or validation: {stats['skipped_pdbs']}\n")
    f.write(f"Total binding sites evaluated: {stats['total_sites']}\n")
    f.write(f"Binding sites with ≥ 50% recovery: {stats['site_hit_count']} ({stats['site_hit_count']/stats['total_sites']:.2%})\n")
    f.write(f"PDBs with at least one good cluster: {stats['pdb_hit_count']}\n")

    mean_best = sum(stats['best_overlaps'].values()) / len(stats['best_overlaps']) if stats['best_overlaps'] else 0
    f.write(f"Mean best overlap per BS: {mean_best:.2f}%\n")

    f.write("Overlap distribution (based on best cluster per binding site):\n")
    for label, count in stats['overlap_distribution'].items():
        f.write(f"  {label}: {count} binding site(s)\n")

    top5 = sorted(stats['best_overlaps'].items(), key=lambda x: x[1], reverse=True)[:5]
    f.write("\nTop 5 best binding site recoveries:\n")
    for entry, score in top5:
        f.write(f"  - {entry}: {score:.2f}%\n")

if __name__ == "__main__": 
    parent_dir = "Final_Batches"
    #These final batches are new and were used for the active learning exploration so they are in the format we want, but we didnt use active learning in the end
    #So I use them for evaluation of the whole PocketQuest on unknown data
    output_file = "batch_evaluation.log" #This file has an important summary in the end showing what we find and what not in top 10 clusters.

    global_stats = {
        "total_pdbs": 0,
        "skipped_pdbs": 0,
        "site_hit_count": 0,
        "total_sites": 0,
        "pdb_hit_count": 0,
        "overlap_distribution": {
            "High (≥70%)": 0,
            "Moderate (50–69%)": 0,
            "Low (30–49%)": 0,
            "Very Low (<30%)": 0,
        },
        "best_overlaps": {}
    }

    with open(output_file, "w") as f:
        for batch_folder in sorted(os.listdir(parent_dir)):
            batch_path = os.path.join(parent_dir, batch_folder)
            if os.path.isdir(batch_path):
                f.write(f"\n\n=== Batch: {batch_folder} ===\n")
                run_evaluation_on_batch(batch_path, f, global_stats)
        summarize(f, global_stats)
