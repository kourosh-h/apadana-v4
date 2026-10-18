"""
Simulation workflow: the .mdp file of every stage and the scripts that run them.

The SAM systems are prepared stepwise before the production runs. The chains of a freshly built
system are often in strained conformations, so the first NVT stages use short time steps that are
made longer from stage to stage:

Dry systems                                     Wet systems
1 energy minimization                           1 solvation of the dry step-5 structure (TIP3P)
2 NVT  1 ns, dt 0.25 fs                         2 energy minimization
3 NVT  1 ns, dt 0.5 fs                          3 NVT  1 ns, dt 0.25 fs
4 NVT  1 ns, dt 1 fs                            4 NVT  1 ns, dt 0.5 fs
5 NVT 20 ns, dt 2 fs (equilibration)            5 NVT  1 ns, dt 1 fs
6 NVT 100 ns, dt 2 fs (production)              6 NPT 10 ns, dt 1 fs, C-rescale, semi-isotropic
                                                7 NVT 20 ns, dt 2 fs (equilibration)
                                                8 NPT 100 ns, dt 2 fs, Parrinello-Rahman (production)

All other run parameters (cut-offs, PME, thermostat, constraints, output frequency) come from the
.mdp templates in data/mdp/.

Position restraints (topol.top): -DPOSRES_SURFACE restrains the whole slab (POSRES_SUR.itp) and is
used in the first dry stages; -DPOSRES restrains the centermost 1 nm of the slab (POSRES_MIDSUR.itp)
and is used from dry stage 5 on and in all wet stages.

A new bare silica slab (surface_mode='cell', coverage 0) is relaxed before chains are grafted on it:
1 energy minimization
2 NPT 10 ns, dt 1 fs, C-rescale, semi-isotropic (only x and y change), no position restraints
"""
import os
from pathlib import Path

MDP_TEMPLATES = Path(__file__).parent / 'data' / 'mdp'

# ----------------------------------------------------------------------
# The stages. Change the lengths or time steps here.
#   folder       : folder of the stage, created by the run script
#   template     : .mdp file in data/mdp used as the starting point
#   tpr          : name of the run input file (GROMACS writes confout.gro, traj_comp.xtc, ... next to it)
#   dt_ps        : time step (ps); None for energy minimization
#   length_ns    : simulated time (ns); None for energy minimization
#   define       : position restraints, -DPOSRES_SURFACE (whole slab) or -DPOSRES (middle layer)
#   new_velocities : True for the first MD stage, which generates velocities at 298.15 K
#   extra        : other .mdp values that differ from the template
# ----------------------------------------------------------------------

DRY_RELAXATION = [
    dict(folder='1_minimization', template='minimization_dry.mdp', tpr='em', dt_ps=None, length_ns=None,
         define='-DPOSRES_SURFACE', new_velocities=False, extra={}),
    dict(folder='2_nvt_dt0.25fs', template='nvt_dry.mdp', tpr='nvt', dt_ps=0.00025, length_ns=1,
         define='-DPOSRES_SURFACE', new_velocities=True, extra={}),
    dict(folder='3_nvt_dt0.5fs', template='nvt_dry.mdp', tpr='nvt', dt_ps=0.0005, length_ns=1,
         define='-DPOSRES_SURFACE', new_velocities=False, extra={}),
    dict(folder='4_nvt_dt1fs', template='nvt_dry.mdp', tpr='nvt', dt_ps=0.001, length_ns=1,
         define='-DPOSRES_SURFACE', new_velocities=False, extra={}),
    dict(folder='5_nvt_equilibration', template='nvt_dry.mdp', tpr='nvt', dt_ps=0.002, length_ns=20,
         define='-DPOSRES', new_velocities=False, extra={}),
]

DRY_PRODUCTION = [
    dict(folder='6_production_nvt', template='nvt_dry.mdp', tpr='nvt', dt_ps=0.002, length_ns=100,
         define='-DPOSRES', new_velocities=False, extra={}),
]

WET_RELAXATION = [
    dict(folder='2_minimization', template='minimization_wet.mdp', tpr='em', dt_ps=None, length_ns=None,
         define='-DPOSRES', new_velocities=False, extra={}),
    dict(folder='3_nvt_dt0.25fs', template='nvt_wet.mdp', tpr='nvt', dt_ps=0.00025, length_ns=1,
         define='-DPOSRES', new_velocities=True, extra={}),
    dict(folder='4_nvt_dt0.5fs', template='nvt_wet.mdp', tpr='nvt', dt_ps=0.0005, length_ns=1,
         define='-DPOSRES', new_velocities=False, extra={}),
    dict(folder='5_nvt_dt1fs', template='nvt_wet.mdp', tpr='nvt', dt_ps=0.001, length_ns=1,
         define='-DPOSRES', new_velocities=False, extra={}),
    # C-rescale barostat with a time constant of 0.2 ps, coupled every 10 steps. GROMACS does not
    # accept a coupling interval of 100 steps with tau_p = 0.2 ps (tau_p >= 5 nstpcouple dt).
    dict(folder='6_npt_crescale', template='npt_crescale_wet.mdp', tpr='npt', dt_ps=0.001, length_ns=10,
         define='-DPOSRES', new_velocities=False, extra={}),
    dict(folder='7_nvt_equilibration', template='nvt_wet.mdp', tpr='nvt', dt_ps=0.002, length_ns=20,
         define='-DPOSRES', new_velocities=False, extra={}),
]

WET_PRODUCTION = [
    dict(folder='8_production_npt', template='npt_parrinello_rahman_wet.mdp', tpr='npt', dt_ps=0.002,
         length_ns=100, define='-DPOSRES', new_velocities=False, extra={}),
]

# Relaxation of a new bare silica slab. The templates are the .mdp files with which the slab of the
# package was relaxed; extra sets the time constants of the thermostat and the barostat (0.2 ps)
# and the Fourier spacing (0.12 nm) used in all other stages. The barostat is coupled every
# 10 steps, as in the wet C-rescale stage (see above).
# The thermostat group is Surface (all atoms of the bare slab) instead of system, because the
# run scripts give GROMACS the index file, which has no group System.
SURFACE_RELAXATION = [
    dict(folder='1_minimization', template='minimization_surface.mdp', tpr='em', dt_ps=None, length_ns=None,
         define=None, new_velocities=False, extra={}),
    dict(folder='2_npt_crescale', template='npt_crescale_surface.mdp', tpr='npt', dt_ps=0.001, length_ns=10,
         define=None, new_velocities=True,
         extra={'tau_t': '0.2', 'tau_p': '0.2', 'nstpcouple': '10', 'fourierspacing': '0.12',
                'tc_grps': 'Surface', 'comm_grps': 'Surface'}),
]

# Name of the chain group in the index file and the .mdp files; removed for a slab without chains.
CHAIN_GROUP = 'Ots'
# Width of a column in the tc_grps, tau_t and ref_t lines of the .mdp files.
GROUP_COLUMN_WIDTH = 13

# Settings of the HPC jobs (SLURM). The account and the partition have to be set by the user,
# through the argument hpc of build_sam_system or in the run scripts.
HPC_DEFAULTS = {
    'account': 'ACCOUNT',
    'partition': 'PARTITION',
    'nodes': 1,
    'ntasks_per_node': 64,
    'cpus_per_task': 2,
    'time_relaxation': '24:00:00',
    'time_production': '36:00:00',
    'modules': 'gromacs',
    'gmx': 'gmx_mpi',
    'launcher': 'srun',
    'mdrun_options': '-dlb yes',
}


# ----------------------------------------------------------------------
# .mdp files
# ----------------------------------------------------------------------

def set_mdp_value(lines, key, value, comment_out=False):
    """
    Set 'key = value' in the lines of an .mdp file, keeping the alignment of the file.
    With comment_out=True the line of the key is turned into a comment instead.
    A key that is not in the file is added at the end.
    """
    # the line that sets the key; a commented-out line is used only when there is no active one
    found = None
    for n, line in enumerate(lines):
        if '=' not in line:
            continue
        name = line.split('=')[0].strip()
        if name == key:
            found = n
            break
        if found is None and name.lstrip(';').strip() == key:
            found = n

    if found is None:
        if not comment_out:
            lines.append(f'{key:<24}= {value}\n')
    elif comment_out:
        if not lines[found].startswith(';'):
            lines[found] = ';' + lines[found]
    else:
        # an active line that already has this value is kept as it is,
        # so that unchanged settings stay identical to the template
        old_value = lines[found].split('=', 1)[1].split(';')[0].strip()
        if lines[found].startswith(';') or old_value != str(value):
            lines[found] = f'{key:<24}= {value}\n'


def get_mdp_value(lines, key):
    """Value of an active (not commented-out) 'key = value' line, or None if there is none."""
    for line in lines:
        if '=' not in line or line.startswith(';'):
            continue
        name, value = line.split('=', 1)
        if name.strip() == key:
            return value.split(';')[0].strip()
    return None


def remove_group(lines, group):
    """
    Remove a group from the thermostat lines (tc_grps with its tau_t and ref_t) and from
    comm_grps. Used for the pristine slab without chains, which has no Ots group; the wet
    systems then have tc_grps = Surface TIP3.
    """
    tc_groups = get_mdp_value(lines, 'tc_grps')
    if tc_groups is not None and group in tc_groups.split():
        position = tc_groups.split().index(group)
        for key in ['tc_grps', 'tau_t', 'ref_t']:
            values = get_mdp_value(lines, key).split()
            del values[position]
            columns = ''
            for value in values[:-1]:
                columns += f'{value:<{GROUP_COLUMN_WIDTH}}'
            columns += values[-1]
            set_mdp_value(lines, key, columns)

    comm_groups = get_mdp_value(lines, 'comm_grps')
    if comm_groups is not None and group in comm_groups.split():
        values = comm_groups.split()
        values.remove(group)
        set_mdp_value(lines, 'comm_grps', ' '.join(values))


def number_of_steps(length_ns, dt_ps):
    """Number of MD steps for a run of length_ns nanoseconds with a time step of dt_ps picoseconds."""
    return int(round(length_ns * 1000 / dt_ps))


def write_stage_mdp(stage, mdp_file, with_chains=True):
    """
    Write the .mdp file of one stage. with_chains=False removes the chain group (Ots)
    from the thermostat and center-of-mass groups, for a slab without chains.
    """
    with open(MDP_TEMPLATES / stage['template']) as template:
        lines = template.readlines()

    # define None: no position restraints, the template line is left as it is
    if stage['define'] is not None:
        set_mdp_value(lines, 'define', stage['define'])
    if stage['dt_ps'] is not None:
        set_mdp_value(lines, 'dt', stage['dt_ps'])
        set_mdp_value(lines, 'nsteps', number_of_steps(stage['length_ns'], stage['dt_ps']))
        if stage['new_velocities']:
            set_mdp_value(lines, 'continuation', 'no')
            set_mdp_value(lines, 'gen-vel', 'yes')
            set_mdp_value(lines, 'gen-temp', '298.15')
            set_mdp_value(lines, 'gen-seed', '-1')
        else:
            set_mdp_value(lines, 'continuation', 'yes')
            set_mdp_value(lines, 'gen-vel', 'no')
            set_mdp_value(lines, 'gen-temp', None, comment_out=True)
            set_mdp_value(lines, 'gen-seed', None, comment_out=True)
    for key, value in stage['extra'].items():
        set_mdp_value(lines, key, value)
    if not with_chains:
        remove_group(lines, CHAIN_GROUP)

    with open(mdp_file, 'w') as mdp:
        mdp.writelines(lines)


def all_stages():
    """All stages, as (state, stage) pairs."""
    pairs = []
    for stage in DRY_RELAXATION + DRY_PRODUCTION:
        pairs.append(('dry', stage))
    for stage in WET_RELAXATION + WET_PRODUCTION:
        pairs.append(('wet', stage))
    return pairs


# ----------------------------------------------------------------------
# Run scripts
# ----------------------------------------------------------------------

FUNCTIONS_SCRIPT = r'''#!/bin/bash
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
'''


def slurm_header(job_name, time, hpc):
    """SLURM lines and environment of one run script."""
    lines = [
        '#!/bin/bash',
        f'#SBATCH --job-name={job_name}',
        f'#SBATCH --account={hpc["account"]}',
        f'#SBATCH --partition={hpc["partition"]}',
        f'#SBATCH --time={time}',
        f'#SBATCH --nodes={hpc["nodes"]}',
        f'#SBATCH --ntasks-per-node={hpc["ntasks_per_node"]}',
        f'#SBATCH --cpus-per-task={hpc["cpus_per_task"]}',
        '',
        'set -e',
    ]
    if hpc['modules']:
        lines += ['# environment modules of the cluster (skipped where the module command does not exist)',
                  'if command -v module > /dev/null 2>&1; then',
                  '    module purge',
                  f'    module load {hpc["modules"]}',
                  'fi']
    lines += [
        'export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}',
        '',
        f'GMX=${{GMX:-{hpc["gmx"]}}}',
        f'LAUNCH=${{LAUNCH-{hpc["launcher"]}}}',
        f'MDRUN_OPTIONS=${{MDRUN_OPTIONS-{hpc["mdrun_options"]}}}',
        'GROMPP_ONLY=${GROMPP_ONLY:-0}',
        '',
        '# folder of this system (SLURM runs a copy of the script, so the submit folder is used)',
        'HERE=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}',
        'cd "$HERE"',
        'source "$HERE/workflow_functions.sh"',
        '',
    ]
    return lines


def stage_commands(state, stages, start_structure, start_checkpoint):
    """run_stage lines for consecutive stages; each stage starts from the output of the one before."""
    if state == 'wet':
        topology, index = '$HERE/system/topol_wet.top', '$HERE/system/index_wet.ndx'
    else:
        topology, index = '$HERE/system/topol.top', '$HERE/system/index.ndx'

    lines = []
    structure, checkpoint = start_structure, start_checkpoint
    for stage in stages:
        folder = f'$HERE/{state}/{stage["folder"]}'
        mdp = f'$HERE/mdp/{state}_{stage["folder"]}.mdp'
        lines.append(f'run_stage "{folder}" "{mdp}" {stage["tpr"]} "{structure}" "{checkpoint}" "{topology}" "{index}"')
        structure = f'{folder}/confout.gro'
        checkpoint = f'{folder}/state.cpt' if stage['dt_ps'] is not None else '-'
    return lines


def write_script(path, lines):
    with open(path, 'w') as script:
        script.write('\n'.join(lines) + '\n')
    os.chmod(path, 0o755)


def write_workflow(output_folder, system_name, job_label, hpc=None, with_chains=True):
    """
    Write the .mdp files and the run scripts of the dry and wet workflow into output_folder.

    output_folder : folder of the system (contains system/)
    system_name   : base name of the structure, e.g. func (system/func.gro)
    job_label     : short label for the SLURM job names, e.g. L8_G2.40
    hpc           : dictionary with SLURM settings; missing keys are taken from HPC_DEFAULTS
    with_chains   : False for a slab without chains (coverage 0); the group Ots is then left
                    out of the .mdp files
    """
    settings = dict(HPC_DEFAULTS)
    if hpc:
        settings.update(hpc)

    mdp_folder = os.path.join(output_folder, 'mdp')
    os.makedirs(mdp_folder, exist_ok=True)
    for state, stage in all_stages():
        write_stage_mdp(stage, os.path.join(mdp_folder, f'{state}_{stage["folder"]}.mdp'), with_chains)

    with open(os.path.join(output_folder, 'workflow_functions.sh'), 'w') as functions:
        functions.write(FUNCTIONS_SCRIPT)

    dry_last_relaxation = f'$HERE/dry/{DRY_RELAXATION[-1]["folder"]}'
    wet_last_relaxation = f'$HERE/wet/{WET_RELAXATION[-1]["folder"]}'

    # dry relaxation: from the built structure
    lines = slurm_header(f'{job_label}_dry_relax', settings['time_relaxation'], settings)
    lines += stage_commands('dry', DRY_RELAXATION, f'$HERE/system/{system_name}.gro', '-')
    write_script(os.path.join(output_folder, 'run_dry_relaxation.sh'), lines)

    # dry production: from the end of the dry relaxation
    lines = slurm_header(f'{job_label}_dry_prod', settings['time_production'], settings)
    lines += stage_commands('dry', DRY_PRODUCTION, f'{dry_last_relaxation}/confout.gro',
                            f'{dry_last_relaxation}/state.cpt')
    write_script(os.path.join(output_folder, 'run_dry_production.sh'), lines)

    # wet relaxation: solvation of the end structure of the dry relaxation, then the wet stages
    lines = slurm_header(f'{job_label}_wet_relax', settings['time_relaxation'], settings)
    lines.append(f'solvate_system "{dry_last_relaxation}/confout.gro" "$HERE/wet/1_solvation"')
    lines += stage_commands('wet', WET_RELAXATION, '$HERE/wet/1_solvation/solvated.gro', '-')
    write_script(os.path.join(output_folder, 'run_wet_relaxation.sh'), lines)

    # wet production: from the end of the wet relaxation
    lines = slurm_header(f'{job_label}_wet_prod', settings['time_production'], settings)
    lines += stage_commands('wet', WET_PRODUCTION, f'{wet_last_relaxation}/confout.gro',
                            f'{wet_last_relaxation}/state.cpt')
    write_script(os.path.join(output_folder, 'run_wet_production.sh'), lines)

    # submit all four jobs; each starts when the one it needs has finished
    lines = [
        '#!/bin/bash',
        '# Submits the four jobs of this system. Run it in this folder: bash submit_all.sh',
        'set -e',
        'cd "$(dirname "$0")"',
        'dry_relax=$(sbatch --parsable run_dry_relaxation.sh)',
        'sbatch --dependency=afterok:$dry_relax run_dry_production.sh',
        'wet_relax=$(sbatch --parsable --dependency=afterok:$dry_relax run_wet_relaxation.sh)',
        'sbatch --dependency=afterok:$wet_relax run_wet_production.sh',
        'echo "submitted: dry relaxation $dry_relax, wet relaxation $wet_relax and both production runs"',
    ]
    write_script(os.path.join(output_folder, 'submit_all.sh'), lines)


def write_surface_workflow(output_folder, system_name, job_label, hpc=None):
    """
    Write the .mdp files and the run script of the relaxation of a new bare silica slab:
    energy minimization, then 10 ns NPT in which only x and y change.
    The relaxed slab is surface/2_npt_crescale/confout.gro; chains are grafted on it with
    build_sam_system(..., surface_mode='user', user_surface=<that file>).
    """
    settings = dict(HPC_DEFAULTS)
    if hpc:
        settings.update(hpc)

    mdp_folder = os.path.join(output_folder, 'mdp')
    os.makedirs(mdp_folder, exist_ok=True)
    for stage in SURFACE_RELAXATION:
        write_stage_mdp(stage, os.path.join(mdp_folder, f'surface_{stage["folder"]}.mdp'), with_chains=False)

    with open(os.path.join(output_folder, 'workflow_functions.sh'), 'w') as functions:
        functions.write(FUNCTIONS_SCRIPT)

    lines = slurm_header(f'{job_label}_surface_relax', settings['time_relaxation'], settings)
    lines += stage_commands('surface', SURFACE_RELAXATION, f'$HERE/system/{system_name}.gro', '-')
    lines += ['echo "relaxed slab: $HERE/surface/' + SURFACE_RELAXATION[-1]['folder'] + '/confout.gro"']
    write_script(os.path.join(output_folder, 'run_surface_relaxation.sh'), lines)


def stage_description(stage):
    """Short description of the run of one stage, for the README."""
    if stage['dt_ps'] is None:
        return 'energy minimization (steepest descent)'
    if 'parrinello' in stage['template']:
        return 'NPT, Parrinello-Rahman, semi-isotropic'
    if 'crescale' in stage['template']:
        return 'NPT, C-rescale, semi-isotropic'
    return 'NVT, V-rescale'


def surface_workflow_table():
    """Markdown table of the stages of the slab relaxation, for the README of a new slab."""
    rows = ['| Stage folder | Run | Time step (fs) | Length (ns) | Restrained silica |',
            '| --- | --- | --- | --- | --- |']
    for stage in SURFACE_RELAXATION:
        if stage['dt_ps'] is None:
            dt, length = '-', '-'
        else:
            dt, length = f'{stage["dt_ps"] * 1000:g}', f'{stage["length_ns"]:g}'
        rows.append(f'| surface/{stage["folder"]} | {stage_description(stage)} | {dt} | {length} | none |')
    return rows


def workflow_table():
    """Markdown table of all stages, for the README of the system."""
    rows = ['| State | Stage folder | Run | Time step (fs) | Length (ns) | Restrained silica |',
            '| --- | --- | --- | --- | --- | --- |']
    for state, stage in all_stages():
        if state == 'wet' and stage is WET_RELAXATION[0]:
            rows.append('| wet | 1_solvation | gmx solvate, TIP3P water | - | - | - |')
        if stage['dt_ps'] is None:
            dt, length = '-', '-'
        else:
            dt, length = f'{stage["dt_ps"] * 1000:g}', f'{stage["length_ns"]:g}'
        if stage['define'] == '-DPOSRES_SURFACE':
            restraint = 'whole slab'
        else:
            restraint = 'middle 1 nm'
        rows.append(f'| {state} | {stage["folder"]} | {stage_description(stage)} | {dt} | {length} | {restraint} |')
    return rows
