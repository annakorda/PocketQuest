import os
import subprocess
import colorsys
import shutil


def visualize_clusters(pdb_id, top_n=1, rotate=False):
    """
    Visualize predicted binding site clusters using Chimera.
    Loads the original structure and highlights top N scoring clusters (by residue)
    using Chimera-compatible :resid.chain format and distinct RGB colors.
    Optionally applies rotation for transmembrane proteins if rotate=True.
    Keeps Chimera open and saves a snapshot image.
    """
    # Read top N clusters from results.log
    log_path = os.path.join("results", "results.log")
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
        top_n = int(top_n)
        sorted_clusters = sorted_clusters[:top_n]

    # Dynamically determine input PDB path
    input_dirs = ["input_pdbs", "."]
    original_input_path = None
    for directory in input_dirs:
        test_path = os.path.join(directory, f"{pdb_id}.pdb")
        if os.path.exists(test_path):
            original_input_path = test_path
            break

    if original_input_path is None:
        raise FileNotFoundError(f"{pdb_id}.pdb not found in any of the expected directories: {input_dirs}")

    base_pdb = os.path.abspath(os.path.join("results", f"{pdb_id}.pdb"))
    if not os.path.exists(base_pdb):
        shutil.copyfile(original_input_path, base_pdb)

    cluster_files = [os.path.join("results", f"cluster_{label}.pdb") for label, _ in sorted_clusters]

    chimera_script_path = os.path.join("results", f"visualize_{pdb_id}.cmd")

    with open(chimera_script_path, "w") as f:
        f.write(f"open {base_pdb}\n")
        f.write("preset apply publication 1\n")

        for i, cluster_file in enumerate(cluster_files):
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
        f.write(f"copy file {os.path.abspath(os.path.join('results', f'{pdb_id}_clusters.png'))}\n")

    subprocess.run(["chimera", chimera_script_path])
