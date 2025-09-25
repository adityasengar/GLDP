#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -euo pipefail

log() {
    printf "[%s] %s\n" "$(date +'%Y-%m-%d %H%M%S')" "$1"
}

log "Starting dependency installation..."

# Create a virtual environment
python3 -m venv venv
log "Virtual environment 'venv' created."

# Activate the virtual environment
source venv/bin/activate
log "Virtual environment activated."

# Install dependencies from requirements_tested.txt
pip install --upgrade pip

# Install non-PyTorch dependencies
pip install numpy MDAnalysis==2.7.0 tqdm==4.66.5 pyyaml==6.0.2 h5py==3.12.1 scikit-learn==1.5.1 mdtraj==1.10.1

# Install PyTorch with CUDA 12.1 support
pip install torch==2.4.0 --index-url https://download.pytorch.org/whl/cu121

# Install torch-geometric and its dependencies (it will find the installed PyTorch)
# These are installed with a specific index URL for CUDA compatibility.
pip install torch-scatter -f https://data.pyg.org/whl/torch-2.4.0+cu121.html
pip install torch-sparse -f https://data.pyg.org/whl/torch-2.4.0+cu121.html
pip install torch-cluster -f https://data.pyg.org/whl/torch-2.4.0+cu121.html
pip install torch-geometric==2.5.3

log "Dependencies installed successfully."
log "To activate the environment: source venv/bin/activate"
log "To deactivate the environment: deactivate"