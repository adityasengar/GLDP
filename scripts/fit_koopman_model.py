
import numpy as np
import h5py
import argparse
import os

def fit_dmd(X, Y, rank=None):
    """
    Fits a Dynamic Mode Decomposition model.
    Finds the best-fit linear operator A such that Y ≈ A @ X.

    Args:
        X (np.ndarray): Matrix of snapshots at time t. Shape (state_dim, num_snapshots).
        Y (np.ndarray): Matrix of snapshots at time t+1. Shape (state_dim, num_snapshots).
        rank (int, optional): The rank for SVD truncation. If None, uses the full rank.
                              Truncating can help de-noise and find dominant modes.

    Returns:
        np.ndarray: The Koopman operator A.
        np.ndarray: Eigenvalues of A.
        np.ndarray: Eigenvectors (modes) of A.
    """
    print(f"Performing DMD on snapshot matrices of shape {X.shape}")

    # Step 1: SVD of the input matrix X
    # U is (state_dim, rank), s is (rank,), Vh is (rank, num_snapshots)
    U, s, Vh = np.linalg.svd(X, full_matrices=False)
    
    # Truncate SVD matrices if a rank is specified
    if rank is not None and rank > 0 and rank < len(s):
        print(f"Truncating SVD to rank {rank}")
        U = U[:, :rank]
        s = s[:rank]
        Vh = Vh[:rank, :]
        
    # Step 2: Compute the Koopman operator A
    # A = Y @ V @ Σ⁻¹ @ U.T
    # V is Vh.T
    # Σ⁻¹ is the inverse of the diagonal matrix of singular values
    s_inv = np.diag(1. / s)
    V = Vh.T
    
    # This is the full operator, which can be very large (state_dim x state_dim)
    A = Y @ V @ s_inv @ U.T
    
    # For analysis, it's often better to work with the low-rank representation
    # A_tilde = U.T @ Y @ V @ s_inv
    # This is a smaller matrix (rank x rank) whose eigenvalues approximate those of A
    A_tilde = U.T @ A @ U

    # Step 3: Eigendecomposition of the low-rank operator
    eigenvalues, eigenvectors = np.linalg.eig(A_tilde)
    
    # Convert eigenvectors of A_tilde back to modes of the full system
    # These are the dynamic modes (Koopman modes)
    dynamic_modes = U @ eigenvectors

    print("DMD computation complete.")
    return A, eigenvalues, dynamic_modes

def main():
    parser = argparse.ArgumentParser(description="Fit a linear Koopman operator using Dynamic Mode Decomposition (DMD).")
    parser.add_argument("--h5_file", type=str, default="latent_reps/pooled_embedding.h5",
                        help="Path to the HDF5 file containing the time-series latent embeddings.")
    parser.add_argument("--dataset_key", type=str, default="pooled_embedding",
                        help="The key within the HDF5 file for the dataset.")
    parser.add_argument("--output_file", type=str, default="koopman_operator.npy",
                        help="Path to save the learned Koopman operator matrix.")
    parser.add_argument("--svd_rank", type=int, default=None,
                        help="Rank to truncate the SVD to. Helps in de-noising. If None, uses full rank.")
    parser.add_argument("--train_fraction", type=float, default=1.0,
                        help="Fraction of data to use for training (e.g., 0.8 for 80%%). Default is 1.0 (all data).")
    parser.add_argument("--rollout_noise_std", type=float, default=0.0,
                        help="Standard deviation of Gaussian noise to add at each step of the rollout. Default is 0.0 (deterministic).")

    args = parser.parse_args()

    print(f"Loading data from {args.h5_file} with key {args.dataset_key}")
    if not os.path.exists(args.h5_file):
        print(f"Error: Input file not found at {args.h5_file}")
        return

    with h5py.File(args.h5_file, 'r') as f:
        data = f[args.dataset_key][:]

    if data.ndim > 2:
        print(f"Original data shape: {data.shape}. Flattening to 2D.")
        data = data.reshape(data.shape[0], -1)
    
    print(f"Loaded data with shape: {data.shape}")
    num_timesteps, state_dim = data.shape

    # --- Create Snapshot Matrices based on train_fraction ---
    if args.train_fraction < 0.0 or args.train_fraction > 1.0:
        raise ValueError("train_fraction must be between 0.0 and 1.0")
    
    num_train_samples = int(num_timesteps * args.train_fraction)
    
    if num_train_samples >= num_timesteps:
        train_data = data
        print("Training Koopman model on the full dataset.")
    else:
        train_data = data[:num_train_samples]
        print(f"Training Koopman model on the first {num_train_samples} samples ({args.train_fraction * 100:.1f}% of the data).")

    X = train_data[:-1, :].T
    Y = train_data[1:, :].T

    # --- Fit the DMD model ---
    A, eigenvalues, modes = fit_dmd(X, Y, rank=args.svd_rank)
    print(f"Successfully computed Koopman operator A of shape: {A.shape}")

    # --- Save Koopman Operator, Eigenvalues, and Modes ---
    koopman_output_dir = os.path.join("latent_reps", "koopman_npy_outputs")
    os.makedirs(koopman_output_dir, exist_ok=True)

    # Extract just the filename from args.output_file (e.g., "koopman_operator.npy")
    base_output_filename = os.path.basename(args.output_file)
    output_prefix = os.path.splitext(base_output_filename)[0]

    # Construct full paths for saving
    operator_filepath = os.path.join(koopman_output_dir, base_output_filename)
    eigenvalues_filepath = os.path.join(koopman_output_dir, f"{output_prefix}_eigenvalues.npy")
    modes_filepath = os.path.join(koopman_output_dir, f"{output_prefix}_modes.npy")

    np.save(operator_filepath, A)
    print(f"Koopman operator saved to {operator_filepath}")
    np.save(eigenvalues_filepath, eigenvalues)
    np.save(modes_filepath, modes)
    print(f"Eigenvalues and modes saved to {eigenvalues_filepath} and {modes_filepath}")

    # --- Perform Koopman Rollout (Prediction) ---
    print(f"\nPerforming Koopman rollout with noise std: {args.rollout_noise_std}...")
    predicted_trajectory = np.zeros_like(data)
    predicted_trajectory[0, :] = data[0, :]

    for i in range(1, num_timesteps):
        deterministic_step = A @ predicted_trajectory[i-1, :]
        noise = np.random.normal(0, args.rollout_noise_std, size=state_dim)
        predicted_trajectory[i, :] = deterministic_step + noise

    print(f"Rollout complete. Predicted trajectory shape: {predicted_trajectory.shape}")

    # --- Save the Rollout Trajectory ---
    rollout_dir = "latent_reps"
    os.makedirs(rollout_dir, exist_ok=True)
    
    train_str = f"train{args.train_fraction:.2f}".replace('.', 'p') if args.train_fraction < 1.0 else "full"
    rank_str = f"rank{args.svd_rank}" if args.svd_rank is not None else "rankfull"
    noise_str = f"noise{args.rollout_noise_std:.2e}".replace('.', 'p') if args.rollout_noise_std > 0 else "noise0"
    rollout_filename = f"koopman_rollout_{train_str}_{rank_str}_{noise_str}.h5"
    rollout_filepath = os.path.join(rollout_dir, rollout_filename)

    with h5py.File(rollout_filepath, 'w') as f:
        f.create_dataset('koopman_rollout', data=predicted_trajectory)
        f.attrs['source_h5_file'] = args.h5_file
        f.attrs['svd_rank'] = args.svd_rank if args.svd_rank is not None else -1
        f.attrs['train_fraction'] = args.train_fraction
        f.attrs['rollout_noise_std'] = args.rollout_noise_std

    print(f"Koopman rollout trajectory saved to {rollout_filepath}")


if __name__ == "__main__":
    main()


