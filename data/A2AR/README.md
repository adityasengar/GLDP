# A2AR Data

This directory contains the pre-computed `pooled_embedding.h5` file for the A2AR system.

### Important: Latent Space Dimensions

To use this file for running the propagation and decoding steps, you **must** configure the `param_template.yaml` to match the dimensions of this specific embedding file.

The total latent space dimension is the product of `output_height` and `output_width`. You can find the dimension of the provided file by inspecting the shape of the `pooled_embedding` dataset within it. For example, if the dataset shape is `(num_frames, 100)`, the total dimension is 100.

Once you know the dimension, update the following parameters in `gldp_repository/scripts/param_template.yaml`:

```yaml
decoder2_settings:
  # ... other settings
  output_height: [VALUE]  # e.g., 50
  output_width: [VALUE]   # e.g., 2
  # Ensure output_height * output_width equals the latent dimension of this file.
```