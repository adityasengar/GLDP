# Scripts Directory

This directory contains all the Python scripts and shell wrappers necessary to run the GLDP pipeline.

## Main Execution Scripts

These are the primary entry points for users.

-   `run_pipeline.sh`: The master script that automates the entire workflow for the Koopman and Neural Network methods, from data preprocessing to final analysis.
-   `install_dependencies.sh`: Sets up a Python virtual environment (`venv`) and installs all required packages with tested versions.

## Pipeline Stages

The pipeline is composed of scripts that handle distinct stages of the workflow.

### 1. Preprocessing

-   `extract_res.py`: Processes a raw PDB and trajectory (XTC/DCD) to extract heavy-atom coordinates and dihedral angle information. It generates a `heavy_chain.pdb` and detailed JSON files.
-   `extract_heavy_xtc.py`: Creates a new trajectory file containing only the heavy atoms from the original trajectory, which is used as the ground truth for analysis.

### 2. Autoencoder Training

-   `train_autoencoder.py`: Trains the graph neural network (GNN) based autoencoder. This script learns to encode all-atom structures into a low-dimensional latent space and decode them back. It produces the core `pooled_embedding.h5` file used by the propagators.

### 3. Latent Space Propagation

These scripts simulate dynamics within the latent space.

-   **Koopman Propagator:**
    -   `fit_koopman_model.py`: Fits a linear Koopman operator to the latent space trajectory and generates a new trajectory ("rollout").
-   **Neural Network Propagator:**
    -   `train_neural_propagator.py`: Trains an MLP to learn the non-linear dynamics of the latent space and generates a new trajectory.
-   **Langevin Propagator (Manual Workflow):**
    -   `train_score_model.py`: Trains a time-conditional diffusion model on the latent embeddings to learn the score function (effective force). This is a prerequisite for the Langevin simulation.
    -   `run_langevin_simulation.py`: Performs the Langevin dynamics simulation in the latent space using the pre-trained score model.

### 4. Decoding and Analysis

-   `decode_trajectory.py`: Takes a latent space trajectory (e.g., from a propagator rollout) and uses the trained decoder to reconstruct the full, all-atom 3D coordinates.
-   `convert_h5_to_xtc.py`: Converts the HDF5 coordinate files produced by the decoder into the standard `.xtc` trajectory format for visualization and analysis.
-   `analyze_rmsd.py`: Compares the generated `.xtc` trajectories against the ground truth heavy-atom trajectory and calculates the average subsequent RMSD to find the best-performing model.

## Configuration

-   `param_template.yaml`: A template file containing all the hyperparameters for the autoencoder training and pipeline configuration. The `run_pipeline.sh` script uses this to generate a run-specific `param.yaml`.
-   `requirements.txt` / `requirements_tested.txt`: Lists of Python dependencies.
