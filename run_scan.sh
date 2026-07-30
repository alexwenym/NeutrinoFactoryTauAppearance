#!/bin/bash
#SBATCH --job-name=nufac_scan
#SBATCH --partition=arguelles_delgado
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=slurm_scan_%j.out
#SBATCH --error=slurm_scan_%j.err

cd /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/awen/analysis/NeutrinoFactoryTauAppearance
. .venv/bin/activate

python - <<'EOF'
from scan import run_scan, plot_scan
import numpy as np
res = run_scan()
np.savez("scan_LE.npz", **res)
plot_scan(res)
EOF
