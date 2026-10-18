#!/bin/bash
# Functions used by the run scripts of this system. Not run directly.
#
# Environment variables (set in the run scripts, can be changed before running):
#   GMX            GROMACS command, e.g. gmx_mpi or gmx
#   LAUNCH         MPI launcher for mdrun, e.g. srun; empty to run without it
#   MDRUN_OPTIONS  extra options of gmx mdrun
#   GROMPP_ONLY=1  only check every stage with grompp, without running MD (for testing)

# run_stage <stage folder> <mdp file> <tpr name> <start structure> <checkpoint or -> <topology> <index file>
#   Runs one stage in its own folder. A stage with confout.gro is finished and is skipped;
#   an interrupted run continues from its state.cpt.
run_stage() {
    local folder=$1 mdp=$2 tpr=$3 start=$4 checkpoint=$5 topology=$6 index=$7

    if [ -f "$folder/confout.gro" ]; then
        echo "***** $folder: finished before, skipped"
        return 0
    fi
    echo "***** $folder: starting"
    mkdir -p "$folder"
    cd "$folder"

    local checkpoint_option=""
    if [ "$checkpoint" != "-" ] && [ -f "$checkpoint" ]; then
        checkpoint_option="-t $checkpoint"
    fi
    if [ ! -f "$tpr.tpr" ]; then
        $GMX grompp -f "$mdp" -c "$start" -r "$start" $checkpoint_option -n "$index" -p "$topology" -o "$tpr.tpr" -maxwarn 2
    fi

    if [ "$GROMPP_ONLY" = "1" ]; then
        # test mode: no MD, the start structure is passed on to the next stage
        cp "$start" confout.gro
    else
        local resume_option=""
        if [ -f state.cpt ]; then
            resume_option="-cpi state.cpt"
        fi
        $LAUNCH $GMX mdrun -s "$tpr.tpr" $resume_option $MDRUN_OPTIONS
    fi
    echo "***** $folder: finished"
    cd "$HERE"
}

# solvate_system <dry structure> <stage folder>
#   Fills the empty space of the box with water (gmx solvate with the 3-site water box of GROMACS),
#   renames the water to the CHARMM TIP3P names of TIP3.itp (TIP3: OH2, H1, H2), writes the number
#   of water molecules for topol_wet.top and adds the group TIP3 to the index file.
solvate_system() {
    local dry_structure=$1 folder=$2

    if [ -f "$folder/solvated.gro" ]; then
        echo "***** $folder: finished before, skipped"
        return 0
    fi
    echo "***** $folder: starting"
    mkdir -p "$folder"
    cd "$folder"

    $GMX solvate -cp "$dry_structure" -cs spc216.gro -o water_spc.gro
    sed -e 's/SOL     OW/TIP3   OH2/' -e 's/SOL    HW1/TIP3    H1/' -e 's/SOL    HW2/TIP3    H2/' \
        water_spc.gro > solvated.gro

    local n_solute n_total n_water
    n_solute=$(sed -n 2p "$dry_structure" | tr -d ' ')
    n_total=$(sed -n 2p solvated.gro | tr -d ' ')
    n_water=$(( (n_total - n_solute) / 3 ))
    echo "***** $n_water water molecules added"

    printf "TIP3   %10d\n" "$n_water" > "$HERE/system/water_molecules.itp"

    cp "$HERE/system/index.ndx" "$HERE/system/index_wet.ndx"
    printf "[ TIP3 ]\n" >> "$HERE/system/index_wet.ndx"
    awk -v first=$((n_solute + 1)) -v last="$n_total" 'BEGIN {
        for (i = first; i <= last; i++) {
            printf "%6d", i
            if ((i - first + 1) % 15 == 0) printf "\n"
        }
        printf "\n\n"
    }' >> "$HERE/system/index_wet.ndx"

    echo "***** $folder: finished"
    cd "$HERE"
}
