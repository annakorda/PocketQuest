# PocketQuest
Ligand Binding Site Prediction using Machine Learning

This repository contains the full machine learning pipeline developed for the final project of the Structural Bioinformatics & Python courses. The aim is to predict protein-ligand binding sites from structure using surface geometry and chemical descriptors around Connolly points, trained with a GPU-accelerated XGBoost model.

Due to size constraints, only the code, processed features, and trained model are provided here. Full raw .pdb datasets and molecular surface outputs can be shared via a Drive folder upon request.

## Project Overview
We built a ligand binding site predictor using structural features derived from surface-accessible points of proteins (Connolly points).

Our dataset includes protein-ligand complexes from BindingDB (article-curated).

Workflow Summary:
1. Connolly point generation via MSMS
2. Feature extraction for each surface point
3. Labeling of binding vs non-binding points
4. XGBoost training with extensive Optuna hyperparameter tuning
5. Final model used for classification of Connolly points.
6. Connolly points predicted to belong in the binding site are clustered.
7. Evaluation of clusters and extraction of final binding site(s).

## Repository Contents
utils/

Core scripts for processing and feature extraction:

feature_extraction.py: Generates 32 features per Connolly point

classes.py: Data objects for PDB parsing

check_files.py: Validates input folder contents before batch runs

MLTraining/

Model training, evaluation, and outputs:

Training.ipynb: Main notebook for merging batches, tuning, and training

xgboost_binding_site_model_*.pkl/json: Final trained model in two formats

xgboost_binding_site_model_*.csv: Feature importance rankings

Preprocess/

Preprocessing, mapping, and annotations:

Preprocessing.ipynb: From initial PDB complex files to curated PDB files with their defined binding sites

binding_sites.csv: Residue-level labels for training

batch_content/: Balanced CSV files listing PDBs per batch

other_files/: UniProt ↔ Panther mappings and functional annotations

cmd_scripts/: PyMOL .cmd scripts for viewing previously defined binding sites

screenshots/: Chimera PNG snapshots of binding sites 

final_panther_class_counts.csv: Breakdown of functional classes

low_quality.txt: Structures removed due to preprocessing issues

MSMS/

Tools for Connolly surface generation:

MSMS: Compiled binary for Connolly point generation

pdb_to_xyzr: Converter script to prepare MSMS input

atmtypenumbers: Atom-type lookup table

ML.py: Main script that converts PDB batch folders into numpy arrays ready for training.

## Input & Labeling Notes
1. Connolly points are calculated using a 1.4 Å probe and mapped to residues.
2. Points within 4.0 Å of any non-H atom of a ligand are labeled as binding.
3. Model trained on >10 million points with a 5.5% positive class ratio.
4. Label in numpy arrays: 0=not in binding-site, 1=in binding-site.




