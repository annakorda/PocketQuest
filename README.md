# PocketQuest  
**Ligand Binding Site Prediction using Surface Geometry and Machine Learning**

PocketQuest is a Python-based machine learning pipeline developed as the final project for the Structural Bioinformatics and Python Programming courses. It predicts protein-ligand binding sites by analyzing 3D molecular surface features at Connolly points and classifies them using an optimized XGBoost model trained on tens of millions of surface points.

---

##  Project Overview

We developed a tool that:
- Extracts **Connolly points** from protein surfaces using MSMS.
- Computes **geometric, chemical, and residue-based features** around each point.
- Labels surface points as **binding vs non-binding** using ligand proximity.
- Trains an XGBoost model with **Optuna hyperparameter optimization**.
- Classifies and **clusters high-probability binding points**.
- Outputs **predicted binding sites as PDB files and Chimera visualizations**.

Our training dataset is based on curated protein-ligand complexes from **BindingDB**, covering a wide diversity of functional classes.

---

## 📂 Repository Contents

### Core Executable

- **`PocketQuest.py`**  
  Runs the full binding site prediction pipeline.
  ```bash
  python PocketQuest.py -in input.pdb -prob 0.9 -distance 3 -size 4 -max_residues 30 -vis 3 -rotate
  ```
---

### `utils/`
Core scripts for feature engineering and internal processing:
- `feature_extraction.py`: Extracts 32 structural, chemical, and residue features per Connolly point.
- `classes.py`: Defines `Protein`, `Point`, `Atom`, and `Cluster` objects used throughout the pipeline.
- `check_files.py`: Ensures command-line inputs are valid (type-safe, accessible).
- `make_pdb.py`: Generates cluster-specific PDBs from predicted residues.
- `visualization.py`: Prepares `.cmd` Chimera scripts to visualize clusters.

---

### `MLTraining/`
Model training and evaluation artifacts:
- `Training.ipynb`: Jupyter notebook for Optuna tuning and model training.
- `xgboost_binding_site_model_*.pkl/json`: Final trained XGBoost model.
- `xgboost_binding_site_model_*.csv`: Feature importance scores.

---

### `Preprocess/`
Resources and scripts for dataset preparation:
- `Preprocessing.ipynb`: Pipeline for preparing labeled binding site datasets.
- `binding_sites.csv`: Residue-level labels used during training.
- `batch_content/`: Balanced `.csv` files listing PDBs per batch.
- `cmd_scripts/`: PyMOL `.cmd` files to view known binding sites.
- `screenshots/`: Chimera images of annotated sites.
- `final_panther_class_counts.csv`: Functional class distribution.
- `other_files/`: UniProt ↔ Panther functional mappings.

---

### `MSMS/`
Utilities for Connolly surface generation:
- `MSMS`: Binary executable for Connolly point generation.
- `pdb_to_xyzr`: Format converter for MSMS input.
- `atmtypenumbers`: Atom type lookup table used during preprocessing.

---

### `ML.py`  
Processes a list of PDBs into a NumPy array of Connolly point features.
```bash
python ML.py pdb_list.txt output.npy
```

### `Assessment/`
PocketQuest evaluation:
- `evaluation.py`: Custom metrics to evaluate predicted clusters.
- `batch_evaluation.log`: Score logs across test batches.
- `class_distribution_topscoring.csv`: Top-scoring cluster distribution by function class.

---

## Input & Labeling

- Connolly points generated with a **1.6 Å probe**.
- Points labeled as binding if within **4.0 Å of a ligand non-H atom**.
- Over **10 million points** used for training (positive class ratio: ~5.5%).
Labels:
- `0 = non-binding`  
- `1 = binding-site point`

---

## Output Files
If you clone the repository and run PocketQuest.py, two new directories will be created, files
(contains temporary files) and **results**:
- `results/*.pdb`: Predicted binding-site residue files.
- `results/results.log`: Cluster scores (based on summed probabilities).
- `results/*.cmd`: Chimera visualization scripts.
- `results/*.png`: Snapshot images (if `-vis` option is enabled).

---

## Dataset Availability

Due to GitHub storage limits:
- Only scripts, logs, and trained model are provided
- Full datasets (raw PDBs, Connolly surfaces) can be shared upon request via Google Drive.




