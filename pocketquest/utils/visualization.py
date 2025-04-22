import os
import subprocess
import colorsys
import shutil

def visualize_clusters(pdb_id, top_n=1, rotate=False, results_dir="results"):
    """
    Visualize predicted binding site clusters using Chimera.
    Saves Chimera .cmd script and screenshot in the given results_dir.
    """

    os.makedirs(results_dir, exist_ok=True)

    # Read top N clusters from results.log
    log_path = os.path.join(results_dir, "results.log")
    all_clusters = []
    with open(log_path, "r") as f:
        for line in f:
            if line.startswith("Cluster"):
                parts = line.split(":")
                label = int(parts[0].split()[1])
                score = float(parts[1].split("=")[1])
                all_clusters.append((label, score))

    # Sort clusters by score and get top N
    sorted_clusters = sorted(all_clusters, key=lambda x: x[1], reverse=True)
    if top_n is not None:
        sorted_clusters = sorted_clusters[:int(top_n)]

    # Find the base PDB file (original structure)
    input_dirs = [results_dir, "."]
    original_input_path = None

    for directory in input_dirs:
        test_path = os.path.join(directory, f"{pdb_id}.pdb")
        if os.path.exists(test_path):
            original_input_path = test_path
            break

    if original_input_path is None:
        raise FileNotFoundError(f"{pdb_id}.pdb not found in any of the expected directories: {input_dirs}")

    base_pdb = os.path.join(results_dir, f"{pdb_id}.pdb")
    if not os.path.exists(base_pdb):
        shutil.copyfile(original_input_path, base_pdb)

    cluster_files = [os.path.join(results_dir, f"{pdb_id}cluster_{label}.pdb") for label, _ in sorted_clusters]
    chimera_script_path = os.path.join(results_dir, f"visualize_{pdb_id}.cmd")

    with open(chimera_script_path, "w") as f:
        f.write(f"open {base_pdb}\n")
        f.write("preset apply publication 1\n")

        for i, cluster_file in enumerate(cluster_files):
            if not os.path.exists(cluster_file):
                print(f"[warning] Skipping missing cluster file: {cluster_file}")
                continue

            hue = (i * 0.07) % 1.0
            r, g, b = colorsys.hsv_to_rgb(hue, 0.8, 0.9)
            color_name = f"cluster_color_{i+1}"
            f.write(f"colordef {color_name} {r:.2f} {g:.2f} {b:.2f}\n")

            residues_seen = set()
            with open(cluster_file, "r") as cf:
                for line in cf:
                    if line.startswith("ATOM"):
                        resnum = line[22:26].strip()
                        chain = line[21].strip()
                        res_id = f":{resnum}.{chain}"
                        if res_id not in residues_seen:
                            f.write(f"surface {res_id}\n")
                            f.write(f"color {color_name} {res_id}\n")
                            f.write(f"surftransparency 50 {res_id}\n")
                            residues_seen.add(res_id)

        if rotate:
            f.write("cofr\n")
            f.write("reset\n")
            f.write("turn x 90\n")
            f.write("turn y 180\n")

        f.write("focus\n")
        screenshot_path = os.path.abspath(os.path.join(results_dir, f"{pdb_id}_clusters.png"))
        f.write(f"copy file {screenshot_path}\n")

    subprocess.run(["chimera", chimera_script_path])
