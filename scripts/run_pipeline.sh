#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -euo pipefail

# --- Configuration ---
readonly PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Define standard deviation values for rollout noise (simplified for this pipeline)
# These are hardcoded for simplicity, but could be made configurable if needed.
readonly KOOPMAN_STDEVS=(0.0001 0.001 0.01)
readonly KOOPMAN_RANKS=(10 20)
readonly NEURAL_STDEVS=(0.0001 0.001 0.01)

# --- Helper Functions ---

log() {
    printf "[%s] %s\n" "$(date +'%Y-%m-%d %H%M%S')" "$1"
}

# --- Main Execution ---
main() {
    if [ "$#" -ne 3 ]; then
        log "Usage: $0 <pdb_filename_or_path> <xtc_filename_or_path> <method>"
        log "  <pdb_filename_or_path>: Filename or path to the input PDB file. If relative, it's resolved from the script's directory."
        log "  <xtc_filename_or_path>: Filename or path to the input XTC file. If relative, it's resolved from the script's directory."
        log "  <method>: 'koopman' or 'neural'."
        exit 1
    fi

    local USER_PDB_INPUT="$1"
    local USER_XTC_INPUT="$2"
    local METHOD="$3"

    # Resolve absolute paths for PDB and XTC inputs
    local USER_PDB_PATH="$(realpath "${PIPELINE_DIR}/${USER_PDB_INPUT}")"
    local USER_XTC_PATH="$(realpath "${PIPELINE_DIR}/${USER_XTC_INPUT}")"

    if [[ ! -f "${USER_PDB_PATH}" ]]; then
        log "Error: PDB file not found at resolved path: ${USER_PDB_PATH}" >&2
        exit 1
    fi
    if [[ ! -f "${USER_XTC_PATH}" ]]; then
        log "Error: XTC file not found at resolved path: ${USER_XTC_PATH}" >&2
        exit 1
    fi
    if [[ "${METHOD}" != "koopman" && "${METHOD}" != "neural" ]]; then
        log "Error: Invalid method '${METHOD}'. Choose 'koopman' or 'neural'." >&2
        exit 1
    fi

    # Create a unique working directory for this run
    local RUN_ID="run_$(date +%Y%m%d_%H%M%S)"
    local WORK_DIR="${PIPELINE_DIR}/${RUN_ID}"
    mkdir -p "${WORK_DIR}/checkpoints" "${WORK_DIR}/latent_reps" "${WORK_DIR}/structures" "${WORK_DIR}/output"

    # Redirect all subsequent output to a log file within the working directory
    local LOG_FILE="${WORK_DIR}/pipeline_run.log"
    exec > >(tee -a "$LOG_FILE") 2>&1
    log "Pipeline output redirected to ${LOG_FILE}"
    log "Starting protein dynamics pipeline with method: ${METHOD}"
    log "Created working directory: ${WORK_DIR}"
    log "Resolved PDB input path: ${USER_PDB_PATH}"
    log "Resolved XTC input path: ${USER_XTC_PATH}"

    # Change to the working directory for all subsequent operations
    cd "${WORK_DIR}"

    # --- Step 1: Extract Residues and Heavy Chain PDB ---
    log "Step 1/X: Extracting heavy atoms and residue data..."
    python "${PIPELINE_DIR}/extract_res.py" \
        --pdb "${USER_PDB_PATH}" \
        --traj "${USER_XTC_PATH}" \
        --pdb_out "heavy_chain.pdb" \
        --json_out "residues_data.json" \
        --condensed_out "condensed.json"
    log "Extraction complete. heavy_chain.pdb, residues_data.json, condensed.json generated."

    # --- Step 1.5: Extract Heavy Atoms from Native XTC ---
    local NATIVE_HEAVY_XTC="native_heavy_chain.xtc"
    log "Step 1.5/X: Extracting heavy atoms from native XTC to ${NATIVE_HEAVY_XTC}..."
    python "${PIPELINE_DIR}/extract_heavy_xtc.py" \
        --original_pdb "${USER_PDB_PATH}" \
        --original_xtc "${USER_XTC_PATH}" \
        --heavy_pdb "heavy_chain.pdb" \
        --output_xtc "${NATIVE_HEAVY_XTC}"
    log "Native heavy atom XTC created."

    # --- Step 2: Generate param.yaml from template ---
    log "Step 2/X: Generating param.yaml from template..."
    local PARAM_TEMPLATE_PATH="${PIPELINE_DIR}/param_template.yaml"
    local PARAM_FILE="param.yaml"

    # Use cat and pipe to sed for robust substitution
    cat "${PARAM_TEMPLATE_PATH}" | sed \
        -e "s|{{HEAVY_CHAIN_PDB_PATH}}|heavy_chain.pdb|g" \
        -e "s|{{RESIDUES_DATA_JSON_PATH}}|residues_data.json|g" \
        -e "s|{{CONDENSED_JSON_PATH}}|condensed.json|g" \
        -e "s|checkpoint_dir: checkpoints|checkpoint_dir: ./checkpoints|g" \
        -e "s|latent_dir: latent_reps|latent_dir: ./latent_reps|g" \
        -e "s|structure_dir: structures|structure_dir: ./structures|g" \
        > "${PARAM_FILE}"

    log "param.yaml generated and updated with correct paths."

    # --- Step 3: Run chebnet_blind.py (HNO Encoder + Decoder2 Training) ---
    log "Step 3/X: Running chebnet_blind.py to train HNO and Decoder2..."
    python "${PIPELINE_DIR}/chebnet_blind.py" --config "${PARAM_FILE}"
    log "chebnet_blind.py execution complete. Models trained and initial outputs generated."

    # --- Step 4: Run Koopman or Neural Propagator Pipeline ---
    if [[ "${METHOD}" == "koopman" ]]; then
        log "Step 4/X: Running Koopman pipeline..."
        for rank in "${KOOPMAN_RANKS[@]}"; do
            for stdev_float in "${KOOPMAN_STDEVS[@]}"; do
                log "--- Koopman Experiment: Rank=${rank}, Stdev=${stdev_float} ---"
                local stdev_str
                if (( $(echo "$stdev_float == 0.0" | bc -l) )); then
                    stdev_str="noise0"
                else
                    stdev_str=$(printf "%.2e" "$stdev_float" | sed 's/\./p/')
                    stdev_str="noise${stdev_str}"
                fi

                local output_prefix="koopman_operator_case_${RUN_ID}_rank${rank}_${stdev_str}"
                local diff_emb_file="latent_reps/koopman_rollout_full_rank${rank}_${stdev_str}.h5"
                local output_h5_file="structures/full_coords_koopman_rollout_rank${rank}_${stdev_str}.h5"
                local output_xtc_file="structures/full_coords_koopman_rollout_rank${rank}_${stdev_str}.xtc"
                local output_h5_key="full_coords_koopman_rollout_rank${rank}_${stdev_str}"

                python "${PIPELINE_DIR}/fit_koopman_model.py" \
                    --h5_file "latent_reps/pooled_embedding.h5" \
                    --dataset_key "pooled_embedding" \
                    --svd_rank "${rank}" \
                    --rollout_noise_std "${stdev_float}" \
                    --output_file "${output_prefix}.npy"

                python "${PIPELINE_DIR}/inference_old.py" \
                    --config "${PARAM_FILE}" \
                    --hno_ckpt "checkpoints/hno_checkpoint.pth" \
                    --decoder2_ckpt "checkpoints/decoder2_checkpoint.pth" \
                    --diff_emb_file "${diff_emb_file}" \
                    --diff_emb_key "koopman_rollout" \
                    --conditioner_x_ref_pt "structures/X_ref_coords.pt" \
                    --conditioner_z_ref_pt "latent_reps/z_ref_embedding.pt" \
                    --output_file "${output_h5_file}" \
                    --output_key "${output_h5_key}"

                python "${PIPELINE_DIR}/convert_h5_to_xtc.py" \
                    --h5_file "${output_h5_file}" \
                    --h5_key "${output_h5_key}" \
                    --ref_pdb "heavy_chain.pdb" \
                    --output_xtc "${output_xtc_file}"
                
                # Copy relevant outputs to the final output directory
                cp "${output_h5_file}" "${WORK_DIR}/output/"
                cp "${output_xtc_file}" "${WORK_DIR}/output/"
            done
        done
        log "Koopman pipeline complete."

    elif [[ "${METHOD}" == "neural" ]]; then
        log "Step 4/X: Running Neural Propagator pipeline..."
        for stdev_float in "${NEURAL_STDEVS[@]}"; do
            log "--- Neural Propagator Experiment: Stdev=${stdev_float} ---"
            local stdev_str
            if (( $(echo "$stdev_float == 0" | bc -l) )); then
                stdev_str="0"
            else
                stdev_str=$(printf "%.2e" "${stdev_float}" | sed 's/\./p/g')
            fi

            local OUTPUT_MODEL_PATH="checkpoints/nn_propagator_stdev${stdev_str}.pth"
            local TRAIN_NN_ROLLOUT_H5="latent_reps/nn_rollout_train0p80_skip1_noise${stdev_str}.h5"
            local INFERENCE_OUTPUT_H5="structures/full_coords_nn_rollout_stdev${stdev_str}.h5"
            local OUTPUT_XTC="structures/full_coords_nn_rollout_stdev${stdev_str}.xtc"

            python "${PIPELINE_DIR}/train_neural_propagator.py" \
                --h5_file "latent_reps/pooled_embedding.h5" \
                --dataset_key "pooled_embedding" \
                --output_model_path "${OUTPUT_MODEL_PATH}" \
                --epochs 200 \
                --batch_size 128 \
                --learning_rate 1e-4 \
                --hidden_dim 512 \
                --train_fraction 0.8 \
                --frame_skip 1 \
                --rollout_noise_std "${stdev_float}"

            python "${PIPELINE_DIR}/inference_old.py" \
                --config "${PARAM_FILE}" \
                --hno_ckpt "checkpoints/hno_checkpoint.pth" \
                --decoder2_ckpt "checkpoints/decoder2_checkpoint.pth" \
                --diff_emb_file "${TRAIN_NN_ROLLOUT_H5}" \
                --diff_emb_key "nn_rollout" \
                --conditioner_x_ref_pt "structures/X_ref_coords.pt" \
                --conditioner_z_ref_pt "latent_reps/z_ref_embedding.pt" \
                    --output_file "${INFERENCE_OUTPUT_H5}" \
                --output_key "full_coords_nn_rollout"

            python "${PIPELINE_DIR}/convert_h5_to_xtc.py" \
                --h5_file "${INFERENCE_OUTPUT_H5}" \
                --h5_key "full_coords_nn_rollout" \
                --ref_pdb "heavy_chain.pdb" \
                --output_xtc "${OUTPUT_XTC}"

            # Copy relevant outputs to the final output directory
            cp "${INFERENCE_OUTPUT_H5}" "${WORK_DIR}/output/"
            cp "${OUTPUT_XTC}" "${WORK_DIR}/output/"
        done
        log "Neural Propagator pipeline complete."
    fi

    # Copy heavy_chain.pdb to the final output directory
    cp "heavy_chain.pdb" "${WORK_DIR}/output/"

    log "Step X/X: Analyzing RMSD of generated trajectories..."
    # Capture the output of analyze_rmsd.py
    local RMSD_ANALYSIS_OUTPUT
    RMSD_ANALYSIS_OUTPUT=$(python "${PIPELINE_DIR}/analyze_rmsd.py" \
        --ref_pdb "heavy_chain.pdb" \
        --native_xtc "${NATIVE_HEAVY_XTC}" \
        --output_dir "${WORK_DIR}/output" \
        --num_frames 100)
    log "${RMSD_ANALYSIS_OUTPUT}"
    log "RMSD analysis complete."

    # Extract and log the best XTC path
    local BEST_XTC_PATH=$(echo "${RMSD_ANALYSIS_OUTPUT}" | grep "^BEST_XTC_PATH:" | cut -d' ' -f2-)
    if [[ -n "${BEST_XTC_PATH}" ]]; then
        log "The XTC file with average subsequent RMSD closest to the native is located at: ${BEST_XTC_PATH}"
    else
        log "Could not determine the best XTC file from RMSD analysis."
    fi

    log "Pipeline finished. All results, including the log, are in ${WORK_DIR}/"
}

main "$@"
