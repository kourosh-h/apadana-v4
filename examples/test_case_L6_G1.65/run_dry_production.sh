#!/bin/bash
#SBATCH --job-name=L6_G1.65_dry_prod
#SBATCH --account=ACCOUNT
#SBATCH --partition=PARTITION
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2

set -e
# environment modules of the cluster (skipped where the module command does not exist)
if command -v module > /dev/null 2>&1; then
    module purge
    module load gromacs
fi
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

GMX=${GMX:-gmx_mpi}
LAUNCH=${LAUNCH-srun}
MDRUN_OPTIONS=${MDRUN_OPTIONS--dlb yes}
GROMPP_ONLY=${GROMPP_ONLY:-0}

# folder of this system (SLURM runs a copy of the script, so the submit folder is used)
HERE=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}
cd "$HERE"
source "$HERE/workflow_functions.sh"

run_stage "$HERE/dry/6_production_nvt" "$HERE/mdp/dry_6_production_nvt.mdp" nvt "$HERE/dry/5_nvt_equilibration/confout.gro" "$HERE/dry/5_nvt_equilibration/state.cpt" "$HERE/system/topol.top" "$HERE/system/index.ndx"
