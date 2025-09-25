# Graph Latent Dynamics Propagator (GLDP)

This repository contains the official implementation for the paper "Beyond Ensembles: Simulating All-Atom Protein Dynamics in a Learned Latent Space". We introduce the Graph Latent Dynamics Propagator (GLDP), a modular framework for simulating long-timescale protein dynamics.

The core idea is to use a pre-trained Graph Neural Network to encode high-dimensional, all-atom protein structures into a simplified latent space. Within this space, we can efficiently simulate the system's temporal evolution using one of several "propagator" methods. The resulting latent trajectory is then decoded back into all-atom structures.

## Visual Overview

The GLDP framework follows a modular encoder-propagator-decoder pipeline:

![GLDP Framework Overview](https://i.imgur.com/your_figure_1_url_here.png)
*(Note: You can replace the URL above with a link to an image of Figure 1 from your paper.)*

1.  **Encoder**: A pre-trained ChebNet GNN maps all-atom coordinates to a low-dimensional latent vector `z(t)`.
2.  **Propagator**: Advances the state in the latent space (`z(t) -> z(t+1)`) using one of three methods:
    *   **Score-Guided Langevin Dynamics**: A physics-informed stochastic simulation.
    *   **Koopman Operator**: A linear, data-driven model.
    *   **Autoregressive Neural Network**: A non-linear, expressive neural network.
3.  **Decoder**: A pre-trained network maps the new latent vector `z(t+1)` back to all-atom coordinates.

## Installation

A virtual environment is strongly recommended. The `install_dependencies.sh` script automates the setup of a `venv` and installs all required packages at their tested versions.

```bash
# Navigate to the scripts directory
cd gldp_repository/scripts

# Run the installation script
./install_dependencies.sh

# Activate the environment to run the pipeline
source venv/bin/activate
```

## Usage

The entire pipeline is orchestrated by the `run_pipeline.sh` script. It handles data preprocessing, model training, simulation, and analysis in a single command.

**Command:**
```bash
./run_pipeline.sh <path_to_pdb> <path_to_xtc> <method>
```

**Arguments:**
*   `<path_to_pdb>`: Path to your input PDB file.
*   `<path_to_xtc>`: Path to your input XTC trajectory file.
*   `<method>`: The propagation method to use. Choose from:
    *   `koopman`: For the Koopman operator.
    *   `neural`: For the autoregressive neural network.
    *   `langevin`: For score-guided Langevin dynamics.

**Example:**
```bash
# Assuming PDB/XTC files are in the current directory
./run_pipeline.sh 7jfl_C.pdb 7jfl_C_prod_R1_fit.xtc koopman
```

Each run creates a unique, timestamped output directory (e.g., `run_20250915_210144/`) containing all logs, models, and generated data. The final analysis, which identifies the best-generated trajectory, will be printed to `pipeline_run.log` inside that directory.

## Repository Structure

```
gldp_repository/
├── data/
│   ├── alanine_dipeptide/
│   │   └── README.md  (placeholder for latent reps)
│   ├── 7jfl_C/
│   │   └── README.md  (placeholder for latent reps)
│   ├── A1AR/
│   │   └── README.md  (placeholder for latent reps)
│   └── A2AR/
│       └── README.md  (placeholder for latent reps)
├── scripts/
│   ├── extract_res.py             # Preprocessing: Extracts heavy atoms and dihedrals
│   ├── chebnet_blind.py           # Trains the GNN Encoder-Decoder
│   ├── fit_koopman_model.py       # Koopman propagator
│   ├── train_neural_propagator.py # Neural Network propagator
│   ├── simulate_dynamics.py       # Langevin propagator (score-based)
│   ├── new_diff.py                # Trains the diffusion/score model for Langevin
│   ├── inference_old.py           # Decodes latent trajectories to 3D structures
│   ├── convert_h5_to_xtc.py       # Converts HDF5 coordinates to XTC format
│   ├── analyze_rmsd.py            # Performs final RMSD analysis
│   ├── run_pipeline.sh            # Main execution script
│   ├── install_dependencies.sh    # Sets up the environment
│   └── param_template.yaml        # Configuration template
└── README.md                      # This file
```

## Workflow Steps

1.  **Preprocessing (`extract_res.py`)**: The pipeline begins by processing the input PDB and XTC files to extract heavy-atom coordinates and dihedral angle information, creating a `heavy_chain.pdb` and standardized JSON files.
2.  **Encoder-Decoder Training (`chebnet_blind.py`)**: A Graph Neural Network autoencoder is trained to learn a mapping between the 3D structure and a low-dimensional latent space.
3.  **Score Model Training (`new_diff.py`)**: (For Langevin method only) A diffusion model is trained on the latent space embeddings to learn the score function, which approximates the energy gradient.
4.  **Latent Space Propagation**: A new trajectory is generated in the latent space using one of the chosen methods:
    *   `simulate_dynamics.py` (Langevin)
    *   `fit_koopman_model.py` (Koopman)
    *   `train_neural_propagator.py` (Neural Network)
5.  **Decoding & Analysis**: The new latent trajectory is decoded back into all-atom 3D coordinates (`inference_old.py`), converted to a standard format (`convert_h5_to_xtc.py`), and evaluated against the native trajectory (`analyze_rmsd.py`).

## Citation

If you use this work, please cite the original paper:

```bibtex
@article{sengar2025gldp,
  title={Beyond Ensembles: Simulating All-Atom Protein Dynamics in a Learned Latent Space},
  author={Antiquus S. Hippocampus and Natalia Cerebro and Amelie P. Amygdale and Aditya Sengar},
  journal={bioRxiv},
  year={2025}
}
```
