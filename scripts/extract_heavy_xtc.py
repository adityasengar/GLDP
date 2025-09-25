import argparse
import MDAnalysis as mda
from pathlib import Path
import sys

def extract_heavy_atoms_xtc(original_pdb: Path, original_xtc: Path, heavy_pdb: Path, output_xtc: Path):
    """
    Extracts only heavy atoms from an original trajectory and saves them to a new XTC file.
    The heavy atom selection is based on the heavy_pdb file.
    """
    try:
        # Load the original universe with full PDB and XTC
        full_universe = mda.Universe(str(original_pdb), str(original_xtc))
        
        # Load the heavy_chain PDB to get the heavy atom selection count
        heavy_atom_universe_for_count = mda.Universe(str(heavy_pdb))
        
        # Select all non-hydrogen atoms from the full universe
        heavy_atoms_selection = full_universe.select_atoms("not element H")
        
        # Get the expected number of heavy atoms from the heavy_chain.pdb
        expected_heavy_atom_count = heavy_atom_universe_for_count.select_atoms("all").n_atoms

        if heavy_atoms_selection.n_atoms != expected_heavy_atom_count:
            print(f"ERROR: Mismatch in heavy atom count. Original XTC selection has {heavy_atoms_selection.n_atoms} heavy atoms, but heavy_chain.pdb has {expected_heavy_atom_count}.", file=sys.stderr)
            sys.exit(1)

        # Write the selected heavy atoms to a new XTC file
        with mda.Writer(str(output_xtc), heavy_atoms_selection.n_atoms) as W:
            for ts in full_universe.trajectory:
                W.write(heavy_atoms_selection)
        
        print(f"Successfully extracted {heavy_atoms_selection.n_atoms} heavy atoms to {output_xtc}")

    except Exception as e:
        print(f"ERROR: Failed to extract heavy atoms from XTC: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Extracts heavy atoms from an XTC trajectory.")
    parser.add_argument("--original_pdb", required=True, type=Path, help="Path to the original full PDB file.")
    parser.add_argument("--original_xtc", required=True, type=Path, help="Path to the original full XTC file.")
    parser.add_argument("--heavy_pdb", required=True, type=Path, help="Path to the heavy_chain.pdb file (for atom count verification).")
    parser.add_argument("--output_xtc", required=True, type=Path, help="Path for the output heavy-atom-only XTC file.")
    args = parser.parse_args()

    if not args.original_pdb.exists():
        print(f"Error: Original PDB file not found at {args.original_pdb}", file=sys.stderr)
        sys.exit(1)
    if not args.original_xtc.exists():
        print(f"Error: Original XTC file not found at {args.original_xtc}", file=sys.stderr)
        sys.exit(1)
    if not args.heavy_pdb.exists():
        print(f"Error: Heavy PDB file not found at {args.heavy_pdb}", file=sys.stderr)
        sys.exit(1)

    extract_heavy_atoms_xtc(args.original_pdb, args.original_xtc, args.heavy_pdb, args.output_xtc)

if __name__ == "__main__":
    main()