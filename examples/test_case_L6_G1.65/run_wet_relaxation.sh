#!/bin/bash
#SBATCH --job-name=L6_G1.65_wet_relax
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

solvate_system "$HERE/dry/5_nvt_equilibration/confout.gro" "$HERE/wet/1_solvation"
run_stage "$HERE/wet/2_minimization" "$HERE/mdp/wet_2_minimization.mdp" em "$HERE/wet/1_solvation/solvated.gro" "-" "$HERE/system/topol_wet.top" "$HERE/system/index_wet.ndx"
run_stage "$HERE/wet/3_nvt_dt0.25fs" "$HERE/mdp/wet_3_nvt_dt0.25fs.mdp" nvt "$HERE/wet/2_minimization/confout.gro" "-" "$HERE/system/topol_wet.top" "$HERE/system/index_wet.ndx"
run_stage "$HERE/wet/4_nvt_dt0.5fs" "$HERE/mdp/wet_4_nvt_dt0.5fs.mdp" nvt "$HERE/wet/3_nvt_dt0.25fs/confout.gro" "$HERE/wet/3_nvt_dt0.25fs/state.cpt" "$HERE/system/topol_wet.top" "$HERE/system/index_wet.ndx"
run_stage "$HERE/wet/5_nvt_dt1fs" "$HERE/mdp/wet_5_nvt_dt1fs.mdp" nvt "$HERE/wet/4_nvt_dt0.5fs/confout.gro" "$HERE/wet/4_nvt_dt0.5fs/state.cpt" "$HERE/system/topol_wet.top" "$HERE/system/index_wet.ndx"
run_stage "$HERE/wet/6_npt_crescale" "$HERE/mdp/wet_6_npt_crescale.mdp" npt "$HERE/wet/5_nvt_dt1fs/confout.gro" "$HERE/wet/5_nvt_dt1fs/state.cpt" "$HERE/system/topol_wet.top" "$HERE/system/index_wet.ndx"
run_stage "$HERE/wet/7_nvt_equilibration" "$HERE/mdp/wet_7_nvt_equilibration.mdp" nvt "$HERE/wet/6_npt_crescale/confout.gro" "$HERE/wet/6_npt_crescale/state.cpt" "$HERE/system/topol_wet.top" "$HERE/system/index_wet.ndx"
