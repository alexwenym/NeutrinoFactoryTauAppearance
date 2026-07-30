#!/bin/bash
#SBATCH --job-name=nufac_ch_nubar
#SBATCH --partition=arguelles_delgado
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --output=slurm_channels_nubar_%j.out
#SBATCH --error=slurm_channels_nubar_%j.err

cd /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/awen/analysis/NeutrinoFactoryTauAppearance
. .venv/bin/activate
python scan.py --channels-nubar
