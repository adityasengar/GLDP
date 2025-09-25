import h5py
import numpy as np
import mdtraj as md
import argparse
import os

def convert_h5_to_xtc(h5_file: str, h5_key: str, ref_pdb: str, output_xtc: str):
    """Converts coordinates from an HDF5 file to an XTC trajectory using a reference PDB.

    Args:
        h5_file (str): Path to the input HDF5 file containing coordinates.
        h5_key (str): Dataset key within the HDF5 file for the coordinates.
        ref_pdb (str): Path to the reference PDB file for topology.
        output_xtc (str): Path for the output XTC trajectory file.
    """
    print(f"Loading coordinates from {h5_file} (key: {h5_key})...")
    try:
        with h5py.File(h5_file, 'r') as f:
            coords_np = f[h5_key][:]
        print(f"Loaded coordinates with shape: {coords_np.shape}")
    except Exception as e:
        print(f"Error loading HDF5 file {h5_file}: {e}")
        return

    print(f"Loading topology from reference PDB: {ref_pdb}...")
    try:
        # Load only topology from PDB, coordinates will come from H5
        # mdtraj expects coordinates in nanometers, so we assume H5 is in Angstroms and convert
        # If your H5 is already in nm, remove the /10.0
        dummy_traj = md.load_pdb(ref_pdb)
        topology = dummy_traj.topology
        print(f"Loaded topology with {topology.n_atoms} atoms and {topology.n_residues} residues.")
    except Exception as e:
        print(f"Error loading PDB file {ref_pdb}: {e}")
        return

    # Ensure coordinates are float32 and in nanometers
    # Assuming H5 coordinates are in Angstroms, convert to nanometers for mdtraj
    coords_nm = coords_np.astype(np.float32) / 10.0

    print(f"Creating mdtraj Trajectory object...")
    try:
        # Create a trajectory object with the loaded coordinates and topology
        # The unitcell_lengths and unitcell_angles are optional, set to None for now
        traj = md.Trajectory(xyz=coords_nm, topology=topology, unitcell_lengths=None, unitcell_angles=None)
        print(f"Trajectory created with {traj.n_frames} frames.")
    except Exception as e:
        print(f"Error creating mdtraj Trajectory: {e}")
        return

    print(f"Saving trajectory to XTC: {output_xtc}...")
    try:
        output_dir = os.path.dirname(output_xtc)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        traj.save_xtc(output_xtc)
        print(f"Successfully converted and saved to {output_xtc}")
    except Exception as e:
        print(f"Error saving XTC file {output_xtc}: {e}")
        return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert HDF5 coordinates to XTC trajectory.")
    parser.add_argument("--h5_file", type=str, required=True, help="Path to the input HDF5 file.")
    parser.add_argument("--h5_key", type=str, default="full_coords_diff", help="Dataset key within the HDF5 file.")
    parser.add_argument("--ref_pdb", type=str, required=True, help="Path to the reference PDB file for topology.")
    parser.add_argument("--output_xtc", type=str, required=True, help="Path for the output XTC trajectory file.")
    
    args = parser.parse_args()

    convert_h5_to_xtc(args.h5_file, args.h5_key, args.ref_pdb, args.output_xtc)





