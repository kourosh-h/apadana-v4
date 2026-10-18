"""
Reading and writing of GROMACS text files.

The files written here are the molecule file (.itp), the topology (topol.top), the index groups
(index.ndx) and the position restraints (POSRES_*.itp). Their names and layout are kept fixed,
because the run scripts and the analysis scripts read them.
"""

# Section headers that can appear in a GROMACS topology file.
TOPOLOGY_SECTIONS = [
    '[ molecules ]',
    '[ system ]',
    '[ angles ]',
    '[ dihedrals ]',
    '[ pairs ]',
    '[ bonds ]',
    '[ atoms ]',
    '[ moleculetype ]',
    '[ atomtypes ]',
    '[ pairtypes ]',
    '[ angletypes ]',
    '[ dihedraltypes ]',
    '[ bondtypes ]',
    '[ defaults ]',
]

# Folder with the force-field files, relative to the system folder.
FORCEFIELD_FOLDER = 'ForceField'

# Number of atom IDs per line in the index file.
IDS_PER_LINE = 15

# Force constant of the position restraints (kJ mol^-1 nm^-2), the same in x, y and z.
POSRES_FORCE_CONSTANT = 1000


# ----------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------

def read_lines_without_comments(file_name):
    """Return the lines of a GROMACS file, skipping empty lines and lines that are only a comment (;)."""
    lines = []
    with open(file_name) as file:
        for line in file:
            content = line.partition(';')[0].strip()
            if content:
                lines.append(line.strip())
    return lines


def read_topology_sections(file_name):
    """
    Split a GROMACS topology file into its sections.

    Returns a dictionary: section header (for example '[ atoms ]') -> list of lines in that section.
    """
    lines = read_lines_without_comments(file_name)
    section_starts = [i for i, line in enumerate(lines) if line in TOPOLOGY_SECTIONS]
    section_starts.append(None)

    sections = {}
    for n, start in enumerate(section_starts[:-1]):
        end = section_starts[n + 1]
        sections[lines[start]] = lines[start + 1:end]
    return sections


# ----------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------

def write_molecule_itp(system_name):
    """
    Write {system_name}.itp, the molecule definition used by topol.top.

    The atoms, charges and masses are taken from {system_name}.top written by foyer.
    For the bonds, pairs, angles and dihedrals only the atom numbers and the function type
    are written; the parameters come from the force-field files included in topol.top.
    """
    sections = read_topology_sections(f'{system_name}.top')

    with open(f'{system_name}.itp', 'w+') as itp:
        itp.writelines('\n')
        itp.writelines('[ moleculetype ]\n')
        itp.writelines('; Name       nrexcl\n')
        for line in sections['[ moleculetype ]']:
            s = line.split()
            itp.writelines('{:>6s}{:>10s}\n'.format(s[0], s[1]))
        itp.writelines('\n')
        itp.writelines('\n')

        itp.writelines('[ atoms ]\n')
        itp.writelines(';   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB      massB\n')
        itp.writelines('; residue    1 SUR rtp SUR q -0.0\n')
        for line in sections['[ atoms ]']:
            s = line.split()
            itp.writelines('{:>6s}{:>14s}{:>7s}{:>7s}{:>7s}{:>7s}{:>14s}{:>14s}{:>6s}{:>5s}{:>12s}\n'.format(
                s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10]))
        itp.writelines('\n')

        itp.writelines('[ bonds ]\n')
        itp.writelines(';    ai     aj funct         c0         c1         c2         c3\n')
        for line in sections['[ bonds ]']:
            s = line.split()
            itp.writelines('{:>6s}{:>8s}{:>7s}\n'.format(s[0], s[1], s[2]))
        itp.writelines('\n')
        itp.writelines('\n')

        itp.writelines('[ pairs ]\n')
        itp.writelines(';    ai     aj funct         c0         c1         c2         c3\n')
        for line in sections['[ pairs ]']:
            s = line.split()
            itp.writelines('{:>6s}{:>8s}{:>7s}\n'.format(s[0], s[1], s[2]))
        itp.writelines('\n')

        # Angles use function type 5 (Urey-Bradley), the CHARMM form.
        itp.writelines('[ angles ]\n')
        itp.writelines(';    ai     aj     ak funct         c0         c1         c2         c3\n')
        for line in sections['[ angles ]']:
            s = line.split()
            itp.writelines('{:>6s}{:>8s}{:>8s}{:>7s}\n'.format(s[0], s[1], s[2], '5'))
        itp.writelines('\n')

        # Dihedrals use function type 9 (multiple proper dihedrals), the CHARMM form.
        itp.writelines('[ dihedrals ]\n')
        itp.writelines(';    ai     aj     ak     al funct         c0         c1         c2         c3         c4         c5\n')
        for line in sections['[ dihedrals ]']:
            s = line.split()
            itp.writelines('{:>6s}{:>8s}{:>8s}{:>8s}{:>8s}\n'.format(s[0], s[1], s[2], s[3], '9'))
        itp.writelines('\n')


def write_topol_top(system_name, file_name='topol.top', with_water=False):
    """
    Write the topology that GROMACS reads (topol.top, or topol_wet.top with water).

    It includes the force fields (CHARMM36 for the chains, INTERFACE for the silica and the
    parameters of the Si-C connection at the grafting point) from the ForceField folder next to it,
    then the molecule file and the position restraints of the silica:
      -DPOSRES          restrains the centermost layer of the slab (POSRES_MIDSUR.itp)
      -DPOSRES_SURFACE  restrains the whole slab (POSRES_SUR.itp, first relaxation stages)
    With with_water=True, the CHARMM TIP3P water is included, and the number of water molecules
    is read from water_molecules.itp, which is written when the system is solvated.
    """
    sections = read_topology_sections(f'{system_name}.top')

    with open(file_name, 'w+') as top:
        top.writelines(';\n')
        top.writelines(';\n')
        top.writelines(';--------------------------------------------------------\n')
        top.writelines('; Topoloy file \n')
        top.writelines(';--------------------------------------------------------\n')
        top.writelines(';\n')
        top.writelines(';\n')
        top.writelines(f'#include "{FORCEFIELD_FOLDER}/charmm36-jul2022.ff/forcefield.itp"\n')
        top.writelines(f'#include "{FORCEFIELD_FOLDER}/interface_silica_emami.ff/ff_interface_emami.itp"\n')
        top.writelines(f'#include "{FORCEFIELD_FOLDER}/charmm_sam_connection_point.ff/ff_charmm_connection_point.itp"\n')
        top.writelines('\n')
        top.writelines(';--------------------------------------------------------\n')
        top.writelines('\n')
        top.writelines(f'#include "{system_name}.itp"\n')
        top.writelines('\n')
        top.writelines('\n')
        top.writelines(';--------------------------------------------------------\n')
        top.writelines('\n')
        top.writelines('#ifdef POSRES\n')
        top.writelines('#include "POSRES_MIDSUR.itp"\n')
        top.writelines('#endif\n')
        top.writelines('\n')
        top.writelines('#ifdef POSRES_SURFACE\n')
        top.writelines('#include "POSRES_SUR.itp"\n')
        top.writelines('#endif\n')
        top.writelines('\n')
        if with_water:
            top.writelines(f'#include "{FORCEFIELD_FOLDER}/charmm_water_TIP3.ff/TIP3.itp"\n')
            top.writelines('\n')
        top.writelines(';--------------------------------------------------------\n')
        top.writelines('\n')
        top.writelines('[ system ]\n')
        top.writelines('; Name\n')
        for line in sections['[ system ]']:
            s = line.split()
            if with_water:
                # the name given to the wet systems
                top.writelines('{:>6s}{:>6s} in water\n'.format(s[0], s[1]))
            else:
                top.writelines('{:>6s}{:>6s}\n'.format(s[0], s[1]))
        top.writelines('\n')
        top.writelines(';--------------------------------------------------------\n')
        top.writelines('\n')
        top.writelines('[ molecules ]\n')
        top.writelines('; Compound    #mols\n')
        for line in sections['[ molecules ]']:
            s = line.split()
            top.writelines('{:>6s}{:>10s}\n'.format(s[0], s[1]))
        if with_water:
            top.writelines('#include "water_molecules.itp"\n')


def number_gro_residues(system_name):
    """
    Give the residues in {system_name}.gro the numbers of {system_name}.itp.

    mbuild writes all chains of the .gro file with the same residue number, while the .itp,
    .pdb and .mol2 files give every chain its own number (silica 1, chains 2, 3, ...).
    GROMACS takes the numbers from the topology, so this only matters for programs that read
    the .gro file alone (VMD, MDAnalysis). Only the residue-number column is changed.
    """
    itp = read_topology_sections(f'{system_name}.itp')
    residue_numbers = []
    for line in itp['[ atoms ]']:
        residue_numbers.append(int(line.split()[2]))

    with open(f'{system_name}.gro') as gro:
        lines = gro.readlines()
    n_atoms = int(lines[1])
    if n_atoms != len(residue_numbers):
        raise ValueError(f'{system_name}.gro has {n_atoms} atoms, {system_name}.itp has {len(residue_numbers)}')

    # the atoms are on lines 3 to n_atoms + 2; the residue number is in the first 5 columns
    for n in range(n_atoms):
        line = lines[n + 2]
        lines[n + 2] = '{:>5d}'.format(residue_numbers[n] % 100000) + line[5:]

    with open(f'{system_name}.gro', 'w') as gro:
        gro.writelines(lines)


def write_atom_ids(file, atom_ids):
    """Write atom IDs into an open index file, IDS_PER_LINE numbers per line."""
    for n, atom_id in enumerate(atom_ids, start=1):
        if n % IDS_PER_LINE == 0:
            file.writelines('{:>6d}\n'.format(atom_id))
        else:
            file.writelines('{:>6d}'.format(atom_id))


def write_index_file(surface_ids, sam_ids, ids_per_chain, middle_layer_ids):
    """
    Write index.ndx with the groups used in the simulations and the analysis:
    Surface (silica), Ots (all chains), Ots1, Ots2, ... (one group per chain)
    and midSurface (the restrained middle layer of the silica).
    A slab without chains (coverage 0) has only Surface and midSurface.
    """
    with open('index.ndx', 'w+') as ndx:
        ndx.writelines('[ Surface ]\n')
        write_atom_ids(ndx, surface_ids)

        if sam_ids:
            ndx.writelines('\n\n[ Ots ]\n')
            write_atom_ids(ndx, sam_ids)

        for chain_number, chain_ids in enumerate(ids_per_chain, start=1):
            ndx.writelines(f'\n\n[ Ots{chain_number} ]\n')
            write_atom_ids(ndx, chain_ids)

        ndx.writelines('\n\n[ midSurface ]\n')
        write_atom_ids(ndx, middle_layer_ids)
        ndx.writelines('\n\n')


def write_posres_file(file_name, atom_ids):
    """Write a position-restraint file for the given atom IDs."""
    k = POSRES_FORCE_CONSTANT
    with open(file_name, 'w+') as posres:
        posres.writelines('[ position_restraints ]\n')
        for atom_id in atom_ids:
            posres.writelines('{:>6d}{:>5d}{:>7d}{:>6d}{:>6d}\n'.format(atom_id, 1, k, k, k))


def write_posres_files(surface_ids, sam_ids, middle_layer_ids):
    """
    Write the three position-restraint files:
    POSRES_SUR.itp (whole silica slab, used during relaxation),
    POSRES_OTS.itp (all chain atoms) and
    POSRES_MIDSUR.itp (centermost 1 nm of the slab, used from the 20 ns equilibration on).
    POSRES_OTS.itp is not written for a slab without chains.
    """
    write_posres_file('POSRES_SUR.itp', surface_ids)
    if sam_ids:
        write_posres_file('POSRES_OTS.itp', sam_ids)
    write_posres_file('POSRES_MIDSUR.itp', middle_layer_ids)
