#!/bin/bash
#SBATCH --job-name=L6_G1.65_dry_relax
#SBATCH --account=ACCOUNT
#SBATCH --partition=PARTITION
#SBATCH --time=24:00:00
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

run_stage "$HERE/dry/1_minimization" "$HERE/mdp/dry_1_minimization.mdp" em "$HERE/system/func.gro" "-" "$HERE/system/topol.top" "$HERE/system/index.ndx"
run_stage "$HERE/dry/2_nvt_dt0.25fs" "$HERE/mdp/dry_2_nvt_dt0.25fs.mdp" nvt "$HERE/dry/1_minimization/confout.gro" "-" "$HERE/system/topol.top" "$HERE/system/index.ndx"
run_stage "$HERE/dry/3_nvt_dt0.5fs" "$HERE/mdp/dry_3_nvt_dt0.5fs.mdp" nvt "$HERE/dry/2_nvt_dt0.25fs/confout.gro" "$HERE/dry/2_nvt_dt0.25fs/state.cpt" "$HERE/system/topol.top" "$HERE/system/index.ndx"
run_stage "$HERE/dry/4_nvt_dt1fs" "$HERE/mdp/dry_4_nvt_dt1fs.mdp" nvt "$HERE/dry/3_nvt_dt0.5fs/confout.gro" "$HERE/dry/3_nvt_dt0.5fs/state.cpt" "$HERE/system/topol.top" "$HERE/system/index.ndx"
run_stage "$HERE/dry/5_nvt_equilibration" "$HERE/mdp/dry_5_nvt_equilibration.mdp" nvt "$HERE/dry/4_nvt_dt1fs/confout.gro" "$HERE/dry/4_nvt_dt1fs/state.cpt" "$HERE/system/topol.top" "$HERE/system/index.ndx"
