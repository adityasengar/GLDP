
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader, random_split, Subset
import numpy as np
import h5py
import argparse
import os

# --- Define the Neural Network Model ---
class PropagatorMLP(nn.Module):
    def __init__(self, input_dim, hidden_dim=512, dropout_rate=0.1):
        """
        A simple MLP to learn the dynamics function z_{t+n} = f(z_t).
        """
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(hidden_dim, input_dim) # Output dimension must match input dimension
        )

    def forward(self, x):
        return self.net(x)

# --- Define a custom Dataset ---
class DynamicsDataset(Dataset):
    def __init__(self, data, frame_skip=1):
        """
        Args:
            data (np.ndarray): The full time-series data of shape (num_timesteps, latent_dim).
            frame_skip (int): The number of steps to predict into the future.
        """
        if frame_skip < 1:
            raise ValueError("frame_skip must be at least 1.")
        # Create (input, target) pairs: (z_t, z_{t + frame_skip})
        self.X = torch.tensor(data[:-frame_skip, :], dtype=torch.float32)
        self.Y = torch.tensor(data[frame_skip:, :], dtype=torch.float32)

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.Y[idx]

def main():
    parser = argparse.ArgumentParser(description="Train a neural network to propagate latent space dynamics.")
    parser.add_argument("--h5_file", type=str, default="latent_reps/pooled_embedding.h5",
                        help="Path to the HDF5 file containing the time-series latent embeddings.")
    parser.add_argument("--dataset_key", type=str, default="pooled_embedding",
                        help="The key within the HDF5 file for the dataset.")
    parser.add_argument("--output_model_path", type=str, default="neural_propagator.pth",
                        help="Path to save the trained neural network model.")
    parser.add_argument("--epochs", type=int, default=100, help="Number of training epochs.")
    parser.add_argument("--batch_size", type=int, default=128, help="Batch size for training.")
    parser.add_argument("--learning_rate", type=float, default=1e-4, help="Learning rate for the optimizer.")
    parser.add_argument("--hidden_dim", type=int, default=512, help="Hidden dimension size in the MLP.")
    parser.add_argument("--train_fraction", type=float, default=0.8,
                        help="Fraction of data to use for training (e.g., 0.8 for 80%%). The rest is used for validation.")
    parser.add_argument("--frame_skip", type=int, default=1,
                        help="Number of frames to skip for the prediction target (z_{t+n}). Default is 1 (predict next frame).")
    parser.add_argument("--rollout_noise_std", type=float, default=0.0,
                        help="Standard deviation of Gaussian noise to add at each step of the rollout. Default is 0.0 (deterministic).")

    args = parser.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # --- Load Data ---
    print(f"Loading data from {args.h5_file}")
    if not os.path.exists(args.h5_file):
        print(f"Error: Input file not found at {args.h5_file}")
        return

    with h5py.File(args.h5_file, 'r') as f:
        data = f[args.dataset_key][:]

    if data.ndim > 2:
        data = data.reshape(data.shape[0], -1)
    
    num_timesteps, latent_dim = data.shape
    print(f"Data loaded. Shape: {data.shape}, Latent Dimension: {latent_dim}")

    # --- Create Dataset and DataLoaders ---
    full_dataset = DynamicsDataset(data, frame_skip=args.frame_skip)
    
    if not 0.0 < args.train_fraction <= 1.0:
        raise ValueError("train_fraction must be between 0.0 and 1.0.")
        
    train_size = int(len(full_dataset) * args.train_fraction)
    val_size = len(full_dataset) - train_size
    
    train_indices = list(range(train_size))
    val_indices = list(range(train_size, len(full_dataset)))

    train_dataset = Subset(full_dataset, train_indices)
    val_dataset = Subset(full_dataset, val_indices)

    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True, num_workers=4, pin_memory=True)
    val_loader = DataLoader(val_dataset, batch_size=args.batch_size, shuffle=False, num_workers=4, pin_memory=True)
    
    print(f"Predicting {args.frame_skip} frame(s) ahead.")
    print(f"Using {len(train_dataset)} samples for training ({args.train_fraction * 100:.1f}%)")
    if val_dataset:
        print(f"Using {len(val_dataset)} samples for validation.")

    # --- Initialize Model, Loss, and Optimizer ---
    model = PropagatorMLP(input_dim=latent_dim, hidden_dim=args.hidden_dim).to(device)
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=args.learning_rate)
    
    print("Model initialized:")
    print(model)

    # --- Training Loop ---
    best_val_loss = float('inf')
    for epoch in range(args.epochs):
        model.train()
        train_loss = 0.0
        for inputs, targets in train_loader:
            inputs, targets = inputs.to(device), targets.to(device)
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, targets)
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * inputs.size(0)
        train_loss /= len(train_loader.dataset)

        val_loss = 0.0
        if val_loader:
            model.eval()
            with torch.no_grad():
                for inputs, targets in val_loader:
                    inputs, targets = inputs.to(device), targets.to(device)
                    outputs = model(inputs)
                    loss = criterion(outputs, targets)
                    val_loss += loss.item() * inputs.size(0)
            val_loss /= len(val_loader.dataset)
        
       	print(f"Epoch {epoch+1}/{args.epochs} | Train Loss: {train_loss:.6f} | Val Loss: {val_loss:.6f}")

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), args.output_model_path)
            print(f"  -> Validation loss improved. Model saved to {args.output_model_path}")

    print("Training complete.")
    print(f"Best validation loss: {best_val_loss:.6f}")

    # --- Perform Neural Network Rollout ---
    print(f"\nPerforming Neural Network rollout with noise std: {args.rollout_noise_std}...")
    model.load_state_dict(torch.load(args.output_model_path))
    model.to(device)
    model.eval()

    predicted_trajectory = np.zeros_like(data)
    current_state = torch.tensor(data[0, :], dtype=torch.float32).to(device)
    predicted_trajectory[0, :] = current_state.cpu().numpy()

    with torch.no_grad():
        for i in range(1, num_timesteps):
            deterministic_step = model(current_state)
            noise = torch.randn_like(deterministic_step) * args.rollout_noise_std
            current_state = deterministic_step + noise
            predicted_trajectory[i, :] = current_state.cpu().numpy()

    print(f"Rollout complete. Predicted trajectory shape: {predicted_trajectory.shape}")

    # --- Save the Rollout Trajectory ---
    rollout_dir = "latent_reps"
    os.makedirs(rollout_dir, exist_ok=True)
    
    train_str = f"train{args.train_fraction:.2f}".replace('.', 'p')
    skip_str = f"skip{args.frame_skip}"
    noise_str = f"noise{args.rollout_noise_std:.2e}".replace('.', 'p') if args.rollout_noise_std > 0 else "noise0"
    rollout_filename = f"nn_rollout_{train_str}_{skip_str}_{noise_str}.h5"
    rollout_filepath = os.path.join(rollout_dir, rollout_filename)

    with h5py.File(rollout_filepath, 'w') as f:
        f.create_dataset('nn_rollout', data=predicted_trajectory)
        f.attrs['source_h5_file'] = args.h5_file
        f.attrs['train_fraction'] = args.train_fraction
        f.attrs['frame_skip'] = args.frame_skip
        f.attrs['rollout_noise_std'] = args.rollout_noise_std
        f.attrs['model_path'] = args.output_model_path

    print(f"Neural network rollout trajectory saved to {rollout_filepath}")


if __name__ == "__main__":
    main()






