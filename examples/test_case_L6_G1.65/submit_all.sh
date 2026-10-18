#!/bin/bash
# Submits the four jobs of this system. Run it in this folder: bash submit_all.sh
set -e
cd "$(dirname "$0")"
dry_relax=$(sbatch --parsable run_dry_relaxation.sh)
sbatch --dependency=afterok:$dry_relax run_dry_production.sh
wet_relax=$(sbatch --parsable --dependency=afterok:$dry_relax run_wet_relaxation.sh)
sbatch --dependency=afterok:$wet_relax run_wet_production.sh
echo "submitted: dry relaxation $dry_relax, wet relaxation $wet_relax and both production runs"
