# Output files

This page describes every file that apadana writes. The test case in
[`examples/test_case_L6_G1.65`](../examples/test_case_L6_G1.65) contains all files of a SAM system
as they are written by a build, and can be used as a reference.

All structure and topology files follow the GROMACS conventions: lengths in nm (Å in the `.pdb`
and `.mol2` files), energies in kJ/mol and the origin of the box at one of its corners. The base
name of the structure files is `func` by default and is set with `system_name`.


## Folder of a SAM system

```
L6_G1.65/
    README.md
    system_overview.pdf
    system_overview.png
    system/
    mdp/
    run_dry_relaxation.sh
    run_dry_production.sh
    run_wet_relaxation.sh
    run_wet_production.sh
    workflow_functions.sh
    submit_all.sh
```

The folders `dry/` and `wet/` with one folder per simulation stage are made when the job scripts
run (see [Files written by the simulations](#files-written-by-the-simulations)).


### README.md

Summary of the system in Markdown: chain length, target and actual coverage, chains per face, box
size, face area, numbers of silica and alkylsilane atoms, number of restrained atoms, position of
the slab center in the box, all build settings, the list of files, the table of the simulation
stages and the instructions for running them.


### system_overview.pdf and system_overview.png

One A4 page with views and data of the built system. Lengths are in nm, and *z* is measured from
the geometrical center of the silica slab.

| Panel | Content |
| --- | --- |
| (a), (b) | side views of the box to scale in the *xz*- and *yz*-planes: silica slab (orange), restrained centermost layer (hatched) and the backbone of every chain (blue), with layer labels |
| (c) | mass density along *z* of the silica (filled) and of the SAM (line) |
| (d), (e) | binding locations of the chains on the top face (filled circles) and the bottom face (open circles), with the free hydroxyl oxygens (gray) |
| (f) | distance from every binding location to its nearest neighbor on the same face, compared with the mean value of a random arrangement of the same density; for a slab without chains, the atomic layers of the slab (Si, bridging O, hydroxyl O and H) |
| tables | system data and the stages of the simulation workflow |


### system/

| File | Content |
| --- | --- |
| `func.gro` | coordinates and box. Every chain has its own residue number, the same as in `func.itp` |
| `func.pdb`, `func.mol2` | the same coordinates, for visualization (for example in VMD) |
| `func.itp` | molecule definition `system1` with all atoms (type, residue, name, charge, mass) and all bonds (function 1), pairs (1), angles (5, Urey-Bradley) and dihedrals (9). Only the atom numbers and function types of the bonded terms are listed; their parameters are taken from the force-field files |
| `func.top` | topology written by foyer when the atom types are assigned, with all parameters. `func.itp` is made from it; it is not used by the simulations. Its header gives the date and the program of the build; the user and computer names are removed |
| `topol.top` | topology of the dry system. It includes the force fields from `ForceField/`, then `func.itp`, then the position restraints: `POSRES_MIDSUR.itp` with `define = -DPOSRES` and `POSRES_SUR.itp` with `define = -DPOSRES_SURFACE` |
| `topol_wet.top` | topology of the wet system: the same as `topol.top`, plus the TIP3P water model and the number of water molecules, which is read from `water_molecules.itp` |
| `index.ndx` | index groups `Surface` (all silica atoms), `Ots` (all chain atoms), `Ots1`, `Ots2`, ... (one group per chain) and `midSurface` (the restrained centermost layer). A slab without chains has only `Surface` and `midSurface` |
| `POSRES_SUR.itp` | position restraints of all silica atoms, 1000 kJ mol<sup>-1</sup> nm<sup>-2</sup> in *x*, *y* and *z* |
| `POSRES_MIDSUR.itp` | position restraints of the silica atoms within 0.5 nm of the slab center (the width is set with `restrained_layer_nm`) |
| `POSRES_OTS.itp` | position restraints of all chain atoms. Not used by the workflow; not written for a slab without chains |
| `report_sizes.txt` | box lengths (nm), coverage, chain length and number of chains per face, as plain `key = value` lines |
| `ForceField/` | the force-field folders included by the topologies: `charmm36-jul2022.ff`, `interface_silica_emami.ff`, `charmm_sam_connection_point.ff` and `charmm_water_TIP3.ff` |

Written when the wet system is solvated (stage `wet/1_solvation`):

| File | Content |
| --- | --- |
| `water_molecules.itp` | the line `TIP3 <number of water molecules>`, included in the `[ molecules ]` section of `topol_wet.top` |
| `index_wet.ndx` | the groups of `index.ndx` and the group `TIP3` with all water atoms |


### mdp/

One run-parameter file for every stage, named after the state and the folder of the stage:

| File | Stage |
| --- | --- |
| `dry_1_minimization.mdp` | energy minimization, whole slab restrained |
| `dry_2_nvt_dt0.25fs.mdp` | NVT, 1 ns, time step 0.25 fs, new velocities at 298.15 K |
| `dry_3_nvt_dt0.5fs.mdp` | NVT, 1 ns, 0.5 fs |
| `dry_4_nvt_dt1fs.mdp` | NVT, 1 ns, 1 fs |
| `dry_5_nvt_equilibration.mdp` | NVT, 20 ns, 2 fs, centermost 1 nm restrained |
| `dry_6_production_nvt.mdp` | NVT production, 100 ns, 2 fs |
| `wet_2_minimization.mdp` | energy minimization of the solvated system |
| `wet_3_nvt_dt0.25fs.mdp` | NVT, 1 ns, 0.25 fs, new velocities at 298.15 K |
| `wet_4_nvt_dt0.5fs.mdp` | NVT, 1 ns, 0.5 fs |
| `wet_5_nvt_dt1fs.mdp` | NVT, 1 ns, 1 fs |
| `wet_6_npt_crescale.mdp` | NPT, 10 ns, 1 fs, C-rescale barostat along *z* |
| `wet_7_nvt_equilibration.mdp` | NVT, 20 ns, 2 fs |
| `wet_8_production_npt.mdp` | NPT production, 100 ns, 2 fs, Parrinello-Rahman barostat along *z* |

The files are made from the templates in `apadana/data/mdp` and the stage tables in
`apadana/workflow.py`. For a slab without chains, the group `Ots` is left out of the thermostat
groups.


### Job scripts

| File | Content |
| --- | --- |
| `run_dry_relaxation.sh` | dry stages 1 to 5, from `system/func.gro` |
| `run_dry_production.sh` | dry stage 6, from the end of dry stage 5 |
| `run_wet_relaxation.sh` | solvation of the end of dry stage 5, then wet stages 2 to 7 |
| `run_wet_production.sh` | wet stage 8, from the end of wet stage 7 |
| `workflow_functions.sh` | the shell functions used by the four scripts: `run_stage` (one stage in its own folder) and `solvate_system` |
| `submit_all.sh` | submits the four scripts to SLURM with their dependencies |

Every job script starts with the SLURM settings (account, partition, run time, nodes, tasks and
CPUs per task) and the modules to load, as given with `hpc` in `build_sam_system`. The following
environment variables change the way the scripts run:

| Variable | Default | Meaning |
| --- | --- | --- |
| `GMX` | `gmx_mpi` | GROMACS command |
| `LAUNCH` | `srun` | MPI launcher of `gmx mdrun`; empty to run without it |
| `MDRUN_OPTIONS` | `-dlb yes` | additional options of `gmx mdrun` |
| `GROMPP_ONLY` | `0` | with `1`, every stage is only checked with `gmx grompp`, and the start structure is passed on to the next stage without running MD |

A stage whose folder already contains `confout.gro` is skipped, and a stage with a `state.cpt` from
an interrupted run continues from this checkpoint.


## Files written by the simulations

| Folder | Content |
| --- | --- |
| `dry/1_minimization/` ... `dry/6_production_nvt/` | one folder per dry stage |
| `wet/1_solvation/` | `water_spc.gro` (output of `gmx solvate`) and `solvated.gro` (the same with the water renamed to residue TIP3 and atoms OH2, H1, H2) |
| `wet/2_minimization/` ... `wet/8_production_npt/` | one folder per wet stage |

Each stage folder holds the run input file (`em.tpr`, `nvt.tpr` or `npt.tpr`) with the
`mdout.mdp` of `gmx grompp`, and the output of `gmx mdrun`: `confout.gro`, `ener.edr`, `md.log`,
`state.cpt`, and `traj_comp.xtc` or `traj.trr` according to the output settings of the stage. The
production trajectories are `dry/6_production_nvt/traj_comp.xtc` and
`wet/8_production_npt/traj_comp.xtc`, with the run input files `nvt.tpr` and `npt.tpr` in the same
folders.


## Folder of a new silica slab

A slab built from unit cells (`surface_mode='cell'`, `coverage=0`) is written with the files for
its relaxation instead of the dry and wet workflow:

| File | Content |
| --- | --- |
| `README.md` | summary of the slab and instructions for its relaxation and for grafting chains on it |
| `system_overview.pdf`, `system_overview.png` | overview sheet of the slab, with its atomic layers in panel (f) |
| `system/` | `func.gro`, `func.pdb`, `func.mol2`, `func.itp`, `func.top`, `topol.top`, `index.ndx` (groups `Surface` and `midSurface`), `POSRES_SUR.itp`, `POSRES_MIDSUR.itp`, `report_sizes.txt` and `ForceField/` |
| `mdp/surface_1_minimization.mdp` | energy minimization |
| `mdp/surface_2_npt_crescale.mdp` | NPT, 10 ns, 1 fs, C-rescale barostat in *x* and *y* only, no position restraints |
| `run_surface_relaxation.sh`, `workflow_functions.sh` | job script of the relaxation and its functions |

The relaxation runs in `surface/1_minimization/` and `surface/2_npt_crescale/`. The relaxed slab is
`surface/2_npt_crescale/confout.gro`, which is given to `build_sam_system` with
`surface_mode='user'` and `user_surface`.


## Atom and residue names

| Residue | Atoms |
| --- | --- |
| `SUR` (silica) | `SI` silicon, `OS` bridging oxygen, `OD` and `HD` oxygen and hydrogen of a hydroxyl group, `OP` oxygen bonded to a chain |
| `OTS` (alkylsilane, one residue per chain) | `SI1` silicon, `O1`, `O2`, `H1`, `H2` hydroxyl groups of the silane, `C1` ... `CL` carbon atoms numbered from the silicon to the chain end, `H3` ... hydrogen atoms of the chain |
| `TIP3` (water) | `OH2`, `H1`, `H2` |
