import argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from pathlib import Path
import sys

def calculate_avg_rmsd_subsequent_frames(xtc_file: Path, ref_pdb_file: Path, num_frames: int = 100) -> float:
    """
    Calculates the average RMSD between subsequent frames for the first `num_frames`
    of a trajectory, using the provided PDB as topology.
    """
    try:
        universe = mda.Universe(str(ref_pdb_file), str(xtc_file))
    except Exception as e:
        print(f"ERROR: Could not load trajectory {xtc_file} with topology {ref_pdb_file}: {e}", file=sys.stderr)
        return float('inf') # Return infinity to indicate failure

    if len(universe.trajectory) < 2:
        print(f"WARNING: Trajectory {xtc_file} has fewer than 2 frames. Cannot calculate subsequent RMSD.", file=sys.stderr)
        return float('inf')

    rmsd_values = []
    # Ensure we don't go beyond the available frames or the requested num_frames
    max_frames_to_analyze = min(num_frames, len(universe.trajectory))

    # Set the first frame as the initial reference for the loop
    universe.trajectory[0]
    ref_coords = universe.select_atoms("all").positions.copy()

    for i in range(1, max_frames_to_analyze):
        universe.trajectory[i]
        current_coords = universe.select_atoms("all").positions
        
        # Calculate RMSD directly
        if ref_coords.shape != current_coords.shape:
            print(f"WARNING: Atom count mismatch in {xtc_file} at frame {i}. Skipping RMSD calculation for this frame.", file=sys.stderr)
            continue

        current_rmsd = np.sqrt(np.mean(np.sum((current_coords - ref_coords)**2, axis=1)))
        rmsd_values.append(current_rmsd)
        
        # Update ref_coords for the next iteration to be the current frame's coordinates
        ref_coords = current_coords.copy()

    if not rmsd_values:
        return float('inf')

    return np.mean(rmsd_values)

def main():
    parser = argparse.ArgumentParser(description="Analyze RMSD between subsequent frames for XTC files.")
    parser.add_argument("--ref_pdb", required=True, type=Path, help="Path to the heavy_chain.pdb file (topology).")
    parser.add_argument("--native_xtc", required=True, type=Path, help="Path to the native XTC file.")
    parser.add_argument("--output_dir", required=True, type=Path, help="Path to the pipeline's output directory containing generated XTCs.")
    parser.add_argument("--num_frames", type=int, default=100, help="Number of initial frames to analyze for RMSD.")
    args = parser.parse_args()

    if not args.ref_pdb.exists():
        print(f"Error: Reference PDB file not found at {args.ref_pdb}", file=sys.stderr)
        sys.exit(1)
    if not args.native_xtc.exists():
        print(f"Error: Native XTC file not found at {args.native_xtc}", file=sys.stderr)
        sys.exit(1)
    if not args.output_dir.is_dir():
        print(f"Error: Output directory not found at {args.output_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Analyzing RMSD for the first {args.num_frames} frames...")

    # 1. Calculate RMSD for the native XTC
    print(f"Calculating RMSD for native XTC: {args.native_xtc.name}")
    native_avg_rmsd = calculate_avg_rmsd_subsequent_frames(args.native_xtc, args.ref_pdb, args.num_frames)
    if native_avg_rmsd == float('inf'):
        print("Error: Could not calculate RMSD for native XTC. Exiting.", file=sys.stderr)
        sys.exit(1)
    print(f"Native XTC average subsequent RMSD: {native_avg_rmsd:.4f}")

    # 2. Calculate RMSD for generated XTC files
    generated_xtc_files = list(args.output_dir.glob("*.xtc"))
    if not generated_xtc_files:
        print(f"No generated XTC files found in {args.output_dir}", file=sys.stderr)
        sys.exit(0)

    generated_rmsds = {}
    for xtc_file in generated_xtc_files:
        print(f"Calculating RMSD for generated XTC: {xtc_file.name}")
        avg_rmsd = calculate_avg_rmsd_subsequent_frames(xtc_file, args.ref_pdb, args.num_frames)
        if avg_rmsd != float('inf'):
            generated_rmsds[xtc_file.name] = avg_rmsd

    if not generated_rmsds:
        print("No valid RMSD values could be calculated for generated XTC files.", file=sys.stderr)
        sys.exit(0)

    # 3. Find the closest generated XTC
    closest_xtc_name = None
    min_diff = float('inf')
    closest_xtc_path = None

    for name, rmsd in generated_rmsds.items():
        diff = abs(rmsd - native_avg_rmsd)
        if diff < min_diff:
            min_diff = diff
            closest_xtc_name = name
            closest_xtc_path = args.output_dir / name

    print("\n--- RMSD Analysis Summary ---")
    print(f"Native XTC ({args.native_xtc.name}) average subsequent RMSD: {native_avg_rmsd:.4f}")
    print("\nGenerated XTC files and their average subsequent RMSD:")
    for name, rmsd in generated_rmsds.items():
        print(f"  - {name}: {rmsd:.4f}")

    if closest_xtc_name:
        print(f"\nThe generated XTC file with average subsequent RMSD closest to the native is: {closest_xtc_name} (Difference: {min_diff:.4f})")
        print(f"BEST_XTC_PATH: {closest_xtc_path.resolve()}") # Print absolute path for easy parsing
    else:
        print("Could not determine the closest generated XTC file.")

if __name__ == "__main__":
    main()