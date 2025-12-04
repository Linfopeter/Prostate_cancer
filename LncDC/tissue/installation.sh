#!/bin/bash

# ===================================================================================
# Bioinformatics Tools Installation Script for LncDC Preprocessing (Conda Version)
#
# This script sets up a dedicated and reproducible Conda environment with all
# necessary tools. It avoids manual compilation and system-wide PATH conflicts.
#
# Run with: bash installation.sh
# ===================================================================================

# Exit immediately if a command exits with a non-zero status.
set -e

# --- 0. Pre-flight Check: Do not run as root! ---
if [ "$EUID" -eq 0 ]; then
    echo "❌ ERROR: This script should not be run with sudo or as the root user."
    echo "Please run it as your normal user: ./installation.sh"
    exit 1
fi

# --- 1. Define Environment and Tool list ---
ENV_NAME="lncdc_env"
TOOLS_TO_INSTALL=(
    "hisat2=2.2.1"
    "star=2.7.5a"
    "cufflinks=2.2.1"
    "stringtie=2.1.2"
    "trim-galore=0.6.5"
    "samtools=1.9"
    "gffread"
    "gffcompare"
    "fastqc"
    "cd-hit"
    "cutadapt"
    "bwa"
    "blast"          # <--- añade blastn (BLAST+ suite)
)

# --- 2. Check for Conda Installation and Install if Necessary ---
echo "--- Checking for Conda installation... ---"
if ! command -v conda &> /dev/null; then
    echo "Conda not found. Proceeding with Miniconda installation."
    MINICONDA_DIR="$HOME/miniconda3"
    MINICONDA_SCRIPT="Miniconda3-latest-Linux-x86_64.sh"

    wget "https://repo.anaconda.com/miniconda/$MINICONDA_SCRIPT" -O "$MINICONDA_SCRIPT"
    bash "$MINICONDA_SCRIPT" -b -u -p "$MINICONDA_DIR"
    rm "$MINICONDA_SCRIPT"

    echo "Initializing Conda for your shell..."
    source "$MINICONDA_DIR/bin/activate"
    conda init

    echo "Miniconda installation complete. Please restart your terminal for changes to take full effect and re-run this script."
    exit 0
else
    echo "✅ Conda is already installed."
fi

# --- 3. Configure Conda Channels and Accept ToS ---
echo "--- Configuring Conda channels and accepting Terms of Service... ---"
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# --- 4. Create or Update the Dedicated Environment ---
echo "--- Setting up the '$ENV_NAME' environment... ---"
if conda env list | grep -q "$ENV_NAME"; then
    echo "Environment '$ENV_NAME' already exists. Removing and re-creating to ensure a clean install."
    conda env remove -n "$ENV_NAME" -y
fi

echo "Creating new environment: '$ENV_NAME'."
conda create -n "$ENV_NAME" python=3.6 -y

echo "--- Installing tools: ${TOOLS_TO_INSTALL[*]} ---"
# Install into the newly created environment
conda install -n "$ENV_NAME" --freeze-installed -y "${TOOLS_TO_INSTALL[@]}"

# --- 5. Finalization ---
echo ""
echo "✅ Script completed successfully!"
echo "A Conda environment named '$ENV_NAME' has been created and populated with all necessary tools."
echo ""
echo "============================= IMPORTANT NEXT STEPS ============================="
echo ""
echo "1. Activate the environment before running your analysis scripts:"
echo "   conda activate $ENV_NAME"
echo ""
echo "2. Once activated, you can run all the tools directly. The pipeline script will now be able to automatically handle indexing."
echo ""
echo "3. When you are finished with your analysis, deactivate the environment:"
echo "   conda deactivate"
echo ""
echo "=============================================================================="
echo ""
