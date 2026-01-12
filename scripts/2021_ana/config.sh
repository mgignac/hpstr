#!/bin/bash
# Central configuration for 2021 analysis job submission
# Source this file in other scripts: source config.sh

# User configuration
USER_NAME="${USER:-$(whoami)}"

# SLURM configuration
SLURM_PARTITION="roma"
SLURM_ACCOUNT="hps"

# Year for processing
YEAR="2021"

# Input data paths
INPUT_DATA_DIR="/sdf/data/hps/physics2021/data/recon/pass5_v9/merged"

# Output configuration
OUTPUT_DIR="/sdf/data/hps/physics2021/data/recon/pass5_v9/testrun"

# File list configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_LISTS_DIR="${SCRIPT_DIR}/data_lists/pass5"
FILE_LIST="${DATA_LISTS_DIR}/all_files.list"

# Log directory
LOG_DIR="${SCRIPT_DIR}/log"

# hpstr configuration file
HPSTR_CONFIG="config_2021.py"
