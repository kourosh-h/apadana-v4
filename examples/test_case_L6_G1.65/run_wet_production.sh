#!/bin/bash
#SBATCH --job-name=L6_G1.65_wet_prod
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

run_stage "$HERE/wet/8_production_npt" "$HERE/mdp/wet_8_production_npt.mdp" npt "$HERE/wet/7_nvt_equilibration/confout.gro" "$HERE/wet/7_nvt_equilibration/state.cpt" "$HERE/system/topol_wet.top" "$HERE/system/index_wet.ndx"
