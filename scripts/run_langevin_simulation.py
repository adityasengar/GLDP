import argparse
import torch
import h5py
import numpy as np
import yaml
import os
import sys
import pathlib # For robust path handling

# --- Copied full model definitions from new_diff.py for self-containment ---
class DiffusionMLPBase(torch.nn.Module):
    """Base class for MLP models to handle time embedding and conservative forces."""
    def __init__(self, is_conservative=False):
        super().__init__()
        self.is_conservative = is_conservative

    @staticmethod
    def compute_forces(energy, positions, training):
        """
        Computes forces by taking the negative gradient of energy with respect to positions.
        Ensures gradients are computed only when training.
        """
        # Ensure positions require gradients
        if not positions.requires_grad:
            positions.requires_grad_(True)

        # Compute gradients
        gradient = torch.autograd.grad(
            outputs=energy,
            inputs=positions,
            grad_outputs=torch.ones_like(energy),
            create_graph=training,
            retain_graph=training,
            allow_unused=True
        )[0]

        # Forces are the negative gradient of potential energy
        return -1 * gradient

    def _prepare_input(self, x, t, current_diffusion_steps):
        # x shape: (B, FlatInputDim)
        # t shape: (B,)
        t_norm = (t.float().unsqueeze(1) / current_diffusion_steps) # Shape (B, 1)
        x_in = torch.cat([x, t_norm], dim=1) # Shape (B, FlatInputDim + 1)
        return x_in

class DiffusionMLP(DiffusionMLPBase):
    def __init__(self, input_dim, hidden_dim, is_conservative=False):
        super().__init__(is_conservative)
        output_dim = 1 if is_conservative else input_dim
        self.net = torch.nn.Sequential(
            torch.nn.Linear(input_dim + 1, hidden_dim), torch.nn.ReLU(),
            torch.nn.Linear(hidden_dim, hidden_dim), torch.nn.ReLU(),
            torch.nn.Linear(hidden_dim, output_dim)
        )
    def forward(self, x, t, current_diffusion_steps): # x is expected to be flattened (B, D_flat)
        x_in = self._prepare_input(x, t, current_diffusion_steps)
        model_output = self.net(x_in)
        if self.is_conservative:
            return self.compute_forces(model_output, x, self.training)
        else:
            return model_output

class DiffusionMLP_v2(DiffusionMLPBase):
    def __init__(self, input_dim, hidden_dim, is_conservative=False):
        super().__init__(is_conservative)
        output_dim = 1 if is_conservative else input_dim
        self.net = torch.nn.Sequential(
            torch.nn.Linear(input_dim + 1, hidden_dim), torch.nn.ReLU(),
            torch.nn.Linear(hidden_dim, hidden_dim), torch.nn.ReLU(),
            torch.nn.Linear(hidden_dim, hidden_dim), torch.nn.ReLU(),
            torch.nn.Linear(hidden_dim, output_dim)
        )
    def forward(self, x, t, current_diffusion_steps): # x is expected to be flattened (B, D_flat)
        x_in = self._prepare_input(x, t, current_diffusion_steps)
        model_output = self.net(x_in)
        if self.is_conservative:
            return self.compute_forces(model_output, x, self.training)
        else:
            return model_output

class DiffusionMLP_v3(DiffusionMLPBase):
    def __init__(self, input_dim, hidden_dim, is_conservative=False):
        super().__init__(is_conservative)
        output_dim = 1 if is_conservative else input_dim
        self.fc1 = torch.nn.Linear(input_dim + 1, hidden_dim * 2)
        self.relu = torch.nn.ReLU()
        self.dropout = torch.nn.Dropout(0.2)
        self.fc2 = torch.nn.Linear(hidden_dim * 2, hidden_dim * 2)
        self.fc3 = torch.nn.Linear(hidden_dim * 2, output_dim)
    def forward(self, x, t, current_diffusion_steps): # x is expected to be flattened (B, D_flat)
        x_in = self._prepare_input(x, t, current_diffusion_steps)
        out = self.relu(self.fc1(x_in))
        out = self.dropout(out)
        out = self.relu(self.fc2(out))
        out = self.dropout(out)
        model_output = self.fc3(out)
        if self.is_conservative:
            return self.compute_forces(model_output, x, self.training)
        else:
            return model_output

class DiffusionConv2D(DiffusionMLPBase):
    """ Conv2D denoiser. Input x expected shape (B, 1, H, W). """
    def __init__(self, input_channels=1, hidden_channels=64, is_conservative=False):
        super().__init__(is_conservative)
        self.input_channels = input_channels
        output_channels = 1 if is_conservative else input_channels

        self.conv1 = torch.nn.Conv2d(input_channels + 1, hidden_channels, kernel_size=3, padding=1)
        self.relu1 = torch.nn.ReLU()
        self.conv2 = torch.nn.Conv2d(hidden_channels, hidden_channels, kernel_size=3, padding=1)
        self.relu2 = torch.nn.ReLU()
        self.conv3 = torch.nn.Conv2d(hidden_channels, output_channels, kernel_size=3, padding=1)

    def forward(self, x, t, current_diffusion_steps):
        B, C, H_in, W_in = x.shape
        t_norm = (t.float().view(B, 1, 1, 1) / current_diffusion_steps)
        t_map = t_norm.expand(B, 1, H_in, W_in)
        x_in = torch.cat([x, t_map], dim=1)

        h = self.relu1(self.conv1(x_in))
        h = self.relu2(self.conv2(h))
        model_output = self.conv3(h)

        if self.is_conservative:
            energy = model_output.sum(dim=[1, 2, 3])
            return self.compute_forces(energy, x, self.training)
        else:
            return model_output

# --- End copied section ---

def linear_beta_schedule(timesteps, beta_start, beta_end, device):
    return torch.linspace(beta_start, beta_end, timesteps, device=device)

def main():
    parser = argparse.ArgumentParser(description="Simulate Langevin Dynamics in Latent Space")
    parser.add_argument("--diff_config", type=str, required=True, help="Path to the param_diff.yaml file used for training the diffusion model.")
    parser.add_argument("--diff_ckpt", type=str, required=True, help="Path to the trained diffusion model checkpoint.")
    
    # New arguments for specifying the starting embedding file
    parser.add_argument("--start_emb_file_name", type=str, help="Absolute path to the HDF5 file containing the starting latent embedding. If provided, --start_emb_dir, --start_exp_idx, and --start_conservative_mode are ignored for filename construction.")
    parser.add_argument("--start_emb_dir", type=str, default="latent_reps", help="(Used if --start_emb_file_name is not provided) Directory containing the starting generated embedding file.")
    parser.add_argument("--start_exp_idx", type=int, default=0, help="Index of the specific embedding within the HDF5 file to use as the starting point. Defaults to 0 (first embedding).")
    parser.add_argument("--start_conservative_mode", action="store_true", help="(Used if --start_emb_file_name is not provided) Set if the starting embedding is from a conservative diffusion run (e.g., generated_embeddings_expX_conservative.h5).")
    parser.add_argument("--start_emb_key", type=str, default="generated_embeddings", help="Dataset key for the starting embedding in HDF5 file.")

    parser.add_argument("--num_steps", type=int, required=True, help="The number of simulation steps to perform.")
    parser.add_argument("--temperature", type=float, required=True, help="A float representing the temperature for the Langevin dynamics diffusion term.")
    parser.add_argument("--output_file", type=str, required=True, help="The name of the output HDF5 file.")
    parser.add_argument("--time_step_size", type=float, default=1e-3, help="The time step size (dt) for Langevin dynamics.")
    parser.add_argument("--noise_level_t", type=int, default=0, help="The integer t for noise_level_sigma (sqrt_one_minus_alphas_cumprod[t]).")
    
    # New argument for save frequency
    parser.add_argument("--save_frequency", type=int, default=1, help="Frequency (in steps) at which to save states to the output HDF5 file. Default is 1 (save every step).")
    # New argument for score clipping
    parser.add_argument("--score_clip_value", type=float, default=None, help="Maximum norm for the score/force vector to prevent instability. If None, no clipping is applied.")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility. If None, no seed is set.")

    args = parser.parse_args()

    # Device setup
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # Set random seed for reproducibility
    if args.seed is not None:
        torch.manual_seed(args.seed)
        np.random.seed(args.seed)
        print(f"Set random seed to {args.seed} for torch and numpy.")

    # --- Determine path to Starting Embedding ---
    if args.start_emb_file_name:
        full_start_emb_path = pathlib.Path(args.start_emb_file_name)
    else:
        # Fallback to old logic if direct filename is not provided
        # The start_exp_idx refers to the index within the generated_embeddings_expX_...h5 file.
        # The filename itself should be based on the experiment index used during new_diff.py training.
        # Assuming the user's previous new_diff.py run used exp_idx 1.
        # This 'diffusion_gen_exp_idx' is for constructing the filename, not selecting the sample.
        # It's a bit ambiguous, but keeping it for backward compatibility if start_emb_file_name is not used.
        diffusion_gen_exp_idx = args.start_exp_idx # Use start_exp_idx for filename construction if no direct filename
        start_emb_filename = f"generated_embeddings_exp{diffusion_gen_exp_idx}"
        if args.start_conservative_mode:
            start_emb_filename += "_conservative.h5"
        else:
            start_emb_filename += "_non_conservative.h5"
        full_start_emb_path = pathlib.Path(args.start_emb_dir) / start_emb_filename

    # Load checkpoint
    if not pathlib.Path(args.diff_ckpt).is_file():
        print(f"Error: Checkpoint file not found at {args.diff_ckpt}.")
        sys.exit(1)
    print(f"Loading model checkpoint from: {args.diff_ckpt}")
    checkpoint = torch.load(args.diff_ckpt, map_location=device)
    
    # Extract parameters from checkpoint
    ckpt_params = checkpoint.get('params', {})
    if not ckpt_params:
        print("Error: 'params' key not found in checkpoint. Cannot infer model configuration.")
        sys.exit(1)

    model_type = ckpt_params.get('model_type', 'mlp_v2')
    hidden_dim = ckpt_params.get('hidden_dim', 1024)
    beta_start = ckpt_params.get('beta_start', 5e-6)
    beta_end = ckpt_params.get('beta_end', 0.03)
    diffusion_steps = ckpt_params.get('diffusion_steps', 1400)
    is_conservative_model = ckpt_params.get('conservative', False) # Assuming 'conservative' flag is saved

    print(f"Parameters loaded from checkpoint: model_type={model_type}, hidden_dim={hidden_dim}, "
          f"beta_start={beta_start}, beta_end={beta_end}, diffusion_steps={diffusion_steps}, "
          f"conservative={is_conservative_model}")

    # Load starting vector
    try:
        with h5py.File(full_start_emb_path, 'r') as f:
            start_vectors = f[args.start_emb_key][:]
        print(f"Loaded starting embeddings from: {full_start_emb_path}")
    except Exception as e:
        print(f"Error loading starting embedding from {full_start_emb_path}: {e}")
        sys.exit(1)
    
    # Validate start_exp_idx
    if args.start_exp_idx >= start_vectors.shape[0] or args.start_exp_idx < 0:
        print(f"Error: start_exp_idx {args.start_exp_idx} is out of bounds for file with {start_vectors.shape[0]} vectors.")
        sys.exit(1)

    # Determine latent_dim from the starting vector
    # The generated embeddings from new_diff.py are (N_gen, num_selected_residues, embedding_dim_per_residue)
    # We need the flattened dimension for MLP input_dim.
    # For Conv2D, the model_input_dim from new_diff.py was (1, H, W) where H*W = R*D_per_res.
    # We need to handle this carefully.
    
    # If the input is (N_gen, R, D_per_res), then latent_dim for MLP is R * D_per_res.
    # For Conv2D, the model_input_dim from new_diff.py was (1, H, W) where H*W = R*D_per_res.
    # We need to pass the correct `input_dim` to the model constructor.
    
    # Let's assume `latent_dim` is the flattened dimension for all models for now.
    # The model classes themselves will handle internal reshaping if they are Conv2D.
    # This is a simplification. A more robust solution for Conv2D would involve passing H,W explicitly.
    latent_dim = start_vectors.shape[1] * start_vectors.shape[2] # Assuming (N, R, D_per_res) -> R*D_per_res
    
    # Instantiate and load diffusion model based on inferred type and conservative flag
    model = None
    if model_type == "mlp":
        model = DiffusionMLP(input_dim=latent_dim, hidden_dim=hidden_dim, is_conservative=is_conservative_model)
    elif model_type == "mlp_v2":
        model = DiffusionMLP_v2(input_dim=latent_dim, hidden_dim=hidden_dim, is_conservative=is_conservative_model)
    elif model_type == "mlp_v3":
        model = DiffusionMLP_v3(input_dim=latent_dim, hidden_dim=hidden_dim, is_conservative=is_conservative_model)
    elif model_type == "conv2d":
        # For Conv2D, we need the original H, W from decoder2_settings in param_diff.yaml
        # Re-read `param_diff.yaml` for `decoder2_settings` to get H, W for Conv2D model instantiation.
        with open(args.diff_config, 'r') as f: # Re-open diff_config to get decoder2_settings
            full_diff_config = yaml.safe_load(f)
        d2_settings = full_diff_config.get('parameters', {}).get('decoder2_settings', {})
        conv2d_h = d2_settings.get('output_height')
        conv2d_w = d2_settings.get('output_width')
        if not conv2d_h or not conv2d_w:
            print("Error: Conv2D model type requires 'decoder2_settings' (output_height, output_width) in param_diff.yaml.")
            sys.exit(1)
        # The input_channels for Conv2D is typically 1 for latent embeddings.
        model = DiffusionConv2D(input_channels=1,
                                hidden_channels=hidden_dim, # Re-using hidden_dim for conv_hidden_channels
                                is_conservative=is_conservative_model)
    else:
        print(f"Error: Unknown model type {model_type} from HDF5 attributes.")
        sys.exit(1)

    # Load checkpoint
    if pathlib.Path(args.diff_ckpt).is_file():
        print(f"Loading model checkpoint from: {args.diff_ckpt}")
        checkpoint = torch.load(args.diff_ckpt, map_location=device)
        model.load_state_dict(checkpoint['model_state_dict'])
    else:
        print(f"Warning: Checkpoint file not found at {args.diff_ckpt}. Running with untrained model.")

    model = model.to(device)
    model.eval()

    # Compute diffusion schedule parameters
    betas = linear_beta_schedule(diffusion_steps, beta_start, beta_end, device)
    alphas = torch.clamp(1.0 - betas, min=1e-9)
    alphas_cumprod = torch.cumprod(alphas, dim=0)
    sqrt_one_minus_alphas_cumprod = torch.sqrt(1.0 - alphas_cumprod)

    # Get noise_level_sigma from the schedule
    if args.noise_level_t >= diffusion_steps or args.noise_level_t < 0:
        print(f"Error: noise_level_t {args.noise_level_t} is out of bounds for diffusion_steps {diffusion_steps}.")
        sys.exit(1)
    noise_level_sigma = sqrt_one_minus_alphas_cumprod[args.noise_level_t]

    # Prepare starting latent vector
    # Select the starting vector using start_exp_idx
    z_t_initial_np = start_vectors[args.start_exp_idx]
    z_t_initial = torch.tensor(z_t_initial_np, dtype=torch.float32).to(device)

    # Reshape z_t_initial for the model's expected input format
    if model_type.startswith("mlp"):
        z_t = z_t_initial.view(-1) # Flatten to (FlatDim,)
    elif model_type == "conv2d":
        # Reshape to (C, H, W) for Conv2D model input. C=1.
        # The original shape of z_t_initial_np is (num_selected_residues, embedding_dim_per_residue)
        # which corresponds to (H, W) if C=1.
        z_t = z_t_initial.unsqueeze(0) # Add channel dimension (1, H, W)
    else:
        z_t = z_t_initial # Should not happen

    trajectory = []
    # Always save the initial state
    trajectory.append(z_t.cpu().numpy())

    print(f"Starting Langevin dynamics simulation for {args.num_steps} steps...")
    print(f"Initial latent vector shape: {z_t.shape}")
    print(f"Noise level sigma (t={args.noise_level_t}): {noise_level_sigma.item():.6f}")
    print(f"Saving states every {args.save_frequency} steps.")
    if args.score_clip_value is not None:
        print(f"Using score clipping with max norm: {args.score_clip_value}")

    # Langevin Dynamics Simulation Loop
    with torch.no_grad(): # Operations within this block are not tracked for gradients unless explicitly enabled
        for step in range(args.num_steps):
            # For score calculation, we need to enable gradients if the model is conservative
            # and we are in eval mode. During training, it's handled by model.train().
            # Here, we are in eval mode, so we need torch.enable_grad() for compute_forces.
            
            # The model expects input_dim, so if z_t is (D,), make it (1, D) or (1, C, H, W)
            if model_type.startswith("mlp"):
                model_input_for_grad = z_t.unsqueeze(0).requires_grad_(True) # (1, FlatDim)
            elif model_type == "conv2d":
                model_input_for_grad = z_t.unsqueeze(0).requires_grad_(True) # (1, C, H, W)
            else:
                model_input_for_grad = z_t.unsqueeze(0).requires_grad_(True)

            t_for_model = torch.full((1,), args.noise_level_t, device=device, dtype=torch.long)
            
            # If the model is conservative, its forward pass will compute forces.
            # If not conservative, it predicts noise.
            if is_conservative_model:
                with torch.enable_grad(): # Ensure gradients are computed for force calculation
                    # The model's forward method will return forces directly if is_conservative=True
                    predicted_forces = model(model_input_for_grad, t_for_model, diffusion_steps).squeeze(0) # (FlatDim,) or (C, H, W)
                # In conservative mode, the model directly predicts forces (negative gradient of energy)
                # So, the 'score' is directly `predicted_forces`.
                score = predicted_forces
            else:
                # Non-conservative model predicts noise
                predicted_noise = model(model_input_for_grad, t_for_model, diffusion_steps).squeeze(0) # (FlatDim,) or (C, H, W)
                score = -predicted_noise / noise_level_sigma

            if step == 0:
                print(f"Initial score/force vector norm: {torch.linalg.norm(score).item():.6f}")

            # --- Score Clipping for Stability ---
            if args.score_clip_value is not None:
                score_norm = torch.linalg.norm(score)
                if score_norm > args.score_clip_value:
                    score = score * (args.score_clip_value / score_norm)

            # Langevin dynamics update
            drift_term = args.time_step_size * score
            diffusion_term = torch.sqrt(torch.tensor(2 * args.time_step_size * args.temperature, device=device)) * torch.randn_like(z_t)
            
            z_t = z_t + drift_term + diffusion_term
            
            # Detach z_t to prevent accumulating graph across steps if it was part of a grad computation
            z_t = z_t.detach()

            if (step + 1) % args.save_frequency == 0:
                trajectory.append(z_t.cpu().numpy())

            if (step + 1) % 500 == 0:
                print(f"  Step {step + 1}/{args.num_steps} completed.")

    # Save trajectory
    output_trajectory_np = np.array(trajectory, dtype=np.float32)
    print(f"Saving trajectory of shape {output_trajectory_np.shape} to {args.output_file}")
    
    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with h5py.File(args.output_file, 'w') as f:
        f.create_dataset('langevin_trajectory', data=output_trajectory_np, compression="gzip")

    print("Langevin dynamics simulation complete.")

if __name__ == "__main__":
    main()
