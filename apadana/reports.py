"""
Text reports written for every built system: report_sizes.txt and README.md.
"""
from .workflow import SURFACE_RELAXATION, surface_workflow_table, workflow_table


def write_size_report(system, chain_length, coverage, n_chains):
    """Write report_sizes.txt with the box size (nm) and the grafting numbers."""
    with open("report_sizes.txt", 'w+') as report:
        report.writelines(f"boxX = {system.box.Lx}\n")
        report.writelines(f"boxY = {system.box.Ly}\n")
        report.writelines(f"boxZ = {system.box.Lz}\n")
        report.writelines(f"coverage = {coverage}\n")
        report.writelines(f"chainLength = {chain_length}\n")
        report.writelines(f"number of grafted chains = {n_chains}\n")


def write_system_readme(readme_path, system_name, system, settings, n_chains, n_surface_atoms, n_sam_atoms,
                        n_middle_atoms, slab_center_nm):
    """
    Write README.md of the system: how the system was built, its main numbers (lengths in nm,
    z measured from the center of the silica slab), the files and the simulation workflow.
    """
    area = system.box.Lx * system.box.Ly
    lines = [
        f"# SAM system L = {settings['chain_length']}, Γ = {settings['coverage']} molecules/nm²",
        "",
        "Alkylsilane SAM on both faces of a hydroxylated α-quartz (001) silica slab, built with apadana.",
        "Every molecule is bonded to a hydroxyl group of the surface chosen at random, and both faces",
        "carry the target coverage. The binding locations of the two faces are drawn independently, and",
        "the chains start in the all-trans conformation, aligned with the surface normal.",
        "",
        "![overview of the system](system_overview.png)",
        "",
        "## System",
        "",
        "| Quantity | Value |",
        "| --- | --- |",
        f"| Chain length L (carbon atoms) | {settings['chain_length']} |",
        f"| Target coverage (molecules/nm^2) | {settings['coverage']} |",
        f"| Chains per face | {n_chains} |",
        f"| Chains in total | {2 * n_chains} |",
        f"| Actual coverage per face (molecules/nm^2) | {n_chains / area:.3f} |",
        f"| Box (nm) | {system.box.Lx:.4f} x {system.box.Ly:.4f} x {system.box.Lz:.4f} |",
        f"| Face area (nm^2) | {area:.3f} |",
        f"| Silica atoms (residue {settings['surface_resname']}) | {n_surface_atoms} |",
        f"| Alkylsilane atoms (residue {settings['sam_resname']}) | {n_sam_atoms} |",
        f"| Restrained middle-layer atoms | {n_middle_atoms} (central {settings['restrained_layer_nm']} nm of the slab) |",
        f"| Center of the silica slab, z (nm) | {slab_center_nm:.4f} |",
        "",
        "The z = 0 plane of the overview sheet corresponds to the z-axial geometrical center of the silica",
        f"slab, which lies at z = {slab_center_nm:.4f} nm in the GROMACS box. Subtract this value from the z",
        "coordinates of the box to measure z from the slab center.",
        "",
        "## Build settings",
        "",
        "| Setting | Value |",
        "| --- | --- |",
    ]
    for key, value in settings.items():
        lines.append(f"| {key} | {value} |")
    lines += [
        "",
        "## Files",
        "",
        "| File | Content |",
        "| --- | --- |",
        "| `system_overview.pdf`, `system_overview.png` | one-page overview sheet: side views, density profile, binding locations on both faces, system data and workflow (lengths in nm, z = 0 at the slab center) |",
        f"| `system/{system_name}.gro` | coordinates (nm) |",
        f"| `system/{system_name}.pdb`, `system/{system_name}.mol2` | coordinates for visualization |",
        f"| `system/{system_name}.itp` | molecule definition: atoms, charges, masses and bonded terms |",
        "| `system/topol.top` | topology of the dry system; includes the force fields of `system/ForceField/` |",
        "| `system/topol_wet.top` | topology of the wet system: the same plus TIP3P water |",
        "| `system/ForceField/` | CHARMM36 (chains), INTERFACE (silica), Si-C connection at the grafting point and CHARMM TIP3P water |",
        "| `system/index.ndx` | groups Surface, Ots, Ots1 ... (one per chain) and midSurface |",
        "| `system/POSRES_SUR.itp` | restraints of the whole slab (`define = -DPOSRES_SURFACE`) |",
        "| `system/POSRES_MIDSUR.itp` | restraints of the middle layer of the slab (`define = -DPOSRES`) |",
        "| `system/POSRES_OTS.itp` | restraints of the chain atoms (not used in the workflow) |",
        "| `system/report_sizes.txt` | box size and grafting numbers |",
        f"| `system/{system_name}.top` | foyer topology from which `{system_name}.itp` was made |",
        "| `mdp/` | the .mdp file of every stage |",
        "| `run_*.sh`, `workflow_functions.sh`, `submit_all.sh` | scripts that run the workflow |",
        "",
        "Written when the scripts run: `dry/` and `wet/` with one folder per stage, and",
        "`system/index_wet.ndx` and `system/water_molecules.itp` when the system is solvated.",
        "",
        "## Simulation workflow",
        "",
        "The system is prepared stepwise before the production runs. Each stage starts from the end",
        "structure (and checkpoint) of the stage before. The wet system is made by solvating the end",
        "structure of dry stage 5 with water.",
        "",
    ]
    lines += workflow_table()
    lines += [
        "",
        "Every stage runs in its own folder, so the GROMACS output keeps its default names",
        "(`confout.gro`, `traj_comp.xtc`, `ener.edr`, `md.log`, `state.cpt`; production runs",
        "`dry/6_production_nvt/nvt.tpr` and `wet/8_production_npt/npt.tpr`).",
        "",
        "## Running on an HPC cluster",
        "",
        "1. Copy this folder to the cluster.",
        "2. Check the SLURM lines at the top of the four `run_*.sh` scripts (account, partition, time).",
        "3. In this folder: `bash submit_all.sh`. It submits",
        "   `run_dry_relaxation.sh`, then `run_dry_production.sh` and `run_wet_relaxation.sh` (both wait",
        "   for the dry relaxation), then `run_wet_production.sh` (waits for the wet relaxation).",
        "",
        "A script can also be submitted alone (`sbatch run_dry_relaxation.sh`). Stages that are",
        "finished are skipped, and an interrupted run continues from its checkpoint, so a job",
        "that hits its time limit can be submitted again.",
        "",
        "To check all input files without running MD (for example on a laptop):",
        "",
        "```",
        "GMX=gmx LAUNCH= GROMPP_ONLY=1 bash run_dry_relaxation.sh",
        "```",
        "",
    ]
    with open(readme_path, 'w+') as readme:
        readme.write("\n".join(lines))


def write_surface_readme(readme_path, system_name, system, settings, n_surface_atoms, slab_center_nm):
    """
    Write README.md of a new bare silica slab (surface_mode='cell', coverage 0): the slab, its
    files, and how it is relaxed and then used to graft chains.
    """
    relaxed_slab = f"surface/{SURFACE_RELAXATION[-1]['folder']}/confout.gro"
    lines = [
        f"# New silica slab of {settings['x_tiles']} x {settings['y_tiles']} unit cells",
        "",
        "Bare hydroxylated α-quartz (001) silica slab, built with apadana from the unit cell. The slab is",
        "relaxed in GROMACS before the alkylsilane chains are grafted on it. This folder holds the slab",
        "and the files for its relaxation.",
        "",
        "![overview of the slab](system_overview.png)",
        "",
        "## Slab",
        "",
        "| Quantity | Value |",
        "| --- | --- |",
        f"| Unit cells in x and y | {settings['x_tiles']} x {settings['y_tiles']} |",
        f"| Box (nm) | {system.box.Lx:.4f} x {system.box.Ly:.4f} x {system.box.Lz:.4f} |",
        f"| Face area (nm^2) | {system.box.Lx * system.box.Ly:.3f} |",
        f"| Silica atoms (residue {settings['surface_resname']}) | {n_surface_atoms} |",
        f"| Center of the silica slab, z (nm) | {slab_center_nm:.4f} |",
        "",
        "## Relaxation",
        "",
        "The relaxation consists of an energy minimization (steepest descent) and a 10 ns NPT",
        "simulation. In the NPT simulation the C-rescale barostat couples only the x and y lengths of",
        "the box (compressibility 4.5·10<sup>-6</sup> bar<sup>-1</sup> in x and y, 0 in z), so that the slab",
        "relaxes in its own plane while the box height stays the same. No atoms are restrained.",
        "",
    ]
    lines += surface_workflow_table()
    lines += [
        "",
        "On an HPC cluster: check the SLURM lines at the top of `run_surface_relaxation.sh`, then",
        "`sbatch run_surface_relaxation.sh`. To check the input files without running MD:",
        "",
        "```",
        "GMX=gmx LAUNCH= GROMPP_ONLY=1 bash run_surface_relaxation.sh",
        "```",
        "",
        "## Grafting chains on the relaxed slab",
        "",
        f"The relaxed slab is `{relaxed_slab}`. Copy it to the folder in which the",
        "systems are built and give it to apadana:",
        "",
        "```python",
        "from apadana import build_sam_system",
        "build_sam_system(chain_length=8, coverage=2.40, surface_mode='user',",
        "                 user_surface='confout.gro', output_folder='L8_G2.40')",
        "```",
        "",
        "## Files",
        "",
        "| File | Content |",
        "| --- | --- |",
        "| `system_overview.pdf`, `system_overview.png` | one-page overview sheet: side views, density profile, hydroxyl groups, atomic layers and slab data (lengths in nm, z = 0 at the slab center) |",
        f"| `system/{system_name}.gro` | coordinates (nm) |",
        f"| `system/{system_name}.pdb`, `system/{system_name}.mol2` | coordinates for visualization |",
        f"| `system/{system_name}.itp` | molecule definition: atoms, charges, masses and bonded terms |",
        "| `system/topol.top` | topology; includes the force fields of `system/ForceField/` |",
        "| `system/ForceField/` | INTERFACE force field of the silica, together with the force fields of the chains and the water |",
        "| `system/index.ndx` | groups Surface and midSurface |",
        "| `system/POSRES_SUR.itp`, `system/POSRES_MIDSUR.itp` | restraints of the whole slab and of its middle layer (not used in the relaxation) |",
        "| `system/report_sizes.txt` | box size |",
        f"| `system/{system_name}.top` | foyer topology from which `{system_name}.itp` was made |",
        "| `mdp/surface_*.mdp` | the .mdp file of each relaxation stage |",
        "| `run_surface_relaxation.sh`, `workflow_functions.sh` | scripts that run the relaxation |",
        "",
        "## Build settings",
        "",
        "| Setting | Value |",
        "| --- | --- |",
    ]
    for key, value in settings.items():
        lines.append(f"| {key} | {value} |")
    lines.append("")
    with open(readme_path, 'w+') as readme:
        readme.write("\n".join(lines))
