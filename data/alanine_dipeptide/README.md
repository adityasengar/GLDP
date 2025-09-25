# Alanine Dipeptide Data

This directory contains the pre-computed latent embeddings for the Alanine Dipeptide system. The files have been compressed and split to meet GitHub's file size limit.

### 1. Reassemble the Files

Before using the data, you must reassemble the split files. For each set of parts, run the following commands from this directory:

```bash
# Reassemble 100_1.h5
cat 100_1.h5.gz.part_* > 100_1.h5.gz
gunzip 100_1.h5.gz

# Reassemble 50_1.h5
cat 50_1.h5.gz.part_* > 50_1.h5.gz
gunzip 50_1.h5.gz

# Reassemble 50_2.h5
cat 50_2.h5.gz.part_* > 50_2.h5.gz
gunzip 50_2.h5.gz
```
After running these commands, you will have the original `100_1.h5`, `50_1.h5`, and `50_2.h5` files.

### 2. Important: Latent Space Dimensions

To use these files for running the propagation and decoding steps, you **must** configure the `param_template.yaml` to match the dimensions of the specific embedding file you are using.

The total latent space dimension is the product of `output_height` and `output_width`. You can find the dimension of a given file by inspecting the shape of the `pooled_embedding` dataset within it. For example, if the dataset shape is `(num_frames, 100)`, the total dimension is 100.

Once you know the dimension, update the following parameters in `gldp_repository/scripts/param_template.yaml`:

```yaml
decoder2_settings:
  # ... other settings
  output_height: [VALUE]  # e.g., 50
  output_width: [VALUE]   # e.g., 2
  # Ensure output_height * output_width equals the latent dimension of the file.
```