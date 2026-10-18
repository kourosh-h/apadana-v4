"""
The hydroxylated alpha-quartz (001) silica slab.

The relaxed slab of the package is 7.86 x 7.66 x 2.23 nm^3, periodic in x and y, with an OH group
area density of 9.56 nm^-2 on each of its two exposed (001) faces and the surface normal along z.
Atom names: SI (silicon), OS (bridging oxygen), OD (hydroxyl oxygen), HD (hydroxyl hydrogen).

Lengths: PDB files are in Angstrom, mbuild works in nm.
"""
import os
from pathlib import Path

import mbuild as mb
import numpy as np
from mbuild.lib.recipes import TiledCompound

SILICA_FOLDER = Path(__file__).parent / 'data' / 'silica'

# Temporary PDB file with the renamed atoms; it is read back by mbuild.
RENAMED_PDB = 'unitcell.pdb'
# Temporary PDB file made from a .gro file of the user ('user' mode).
USER_SURFACE_PDB = 'user_surface.pdb'

# Unit cells in x and y in 'cell' mode when none are given: 8 x 9 cells of 0.98260 x 0.85096 nm
# give a slab of 7.86 x 7.66 nm with 12096 atoms, the size of the relaxed slab of the package.
DEFAULT_CELL_TILES = (8, 9)

# Largest distances (nm) at which bonds are added between silica atoms.
SI_O_BOND_CUTOFF = 0.19
O_H_BOND_CUTOFF = 0.120


class SilicaSlab:
    """
    Loads or builds the silica slab.

    Surface modes:
      'slab' : the relaxed slab included in the package (default). The slab was relaxed before
               the chains are grafted on it.
      'cell' : the alpha-quartz unit cell of the package (CHARMM-GUI), repeated x_tiles x y_tiles
               times, to make a new slab. A new slab is relaxed in GROMACS before chains are
               grafted on it (build_sam_system with coverage 0 writes the files for that).
      'user' : a .pdb or .gro file of a slab given by the user, for example a relaxed new slab.

    x_tiles, y_tiles : number of repeats of the structure in x and y
    vacuum_space_nm  : empty space (nm) added above and below the slab in z
    resname          : residue name of the slab in the output files (SUR)
    """

    def __init__(self, name='SilicaSUR',
                 x_tiles=1,
                 y_tiles=1,
                 vacuum_space_nm=6,
                 resname='SUR',
                 silicon_name='SI',
                 bridging_oxygen_name='OS',
                 hydroxyl_oxygen_name='OD',
                 hydroxyl_hydrogen_name='HD'):
        self.name = name
        self.x_tiles = x_tiles
        self.y_tiles = y_tiles
        self.resname = resname
        # total empty space in z: vacuum_space_nm on each side of the slab
        self.total_vacuum_nm = vacuum_space_nm * 2
        self.silicon_name = silicon_name
        self.bridging_oxygen_name = bridging_oxygen_name
        self.hydroxyl_oxygen_name = hydroxyl_oxygen_name
        self.hydroxyl_hydrogen_name = hydroxyl_hydrogen_name

        self.unit_cell_pdb = SILICA_FOLDER / 'unitcell.pdb'
        self.relaxed_slab_pdb = SILICA_FOLDER / 'slab.pdb'

    def read_pdb_dimensions(self, pdb_path):
        """
        Read the size of a PDB structure.

        Returns (box_x, box_y, box_z, empty_space_nm, slab_thickness), where box_x, box_y, box_z
        are from the CRYST1 line (Angstrom), slab_thickness is the z range of the atoms (Angstrom)
        and empty_space_nm is half of the box height not filled by atoms (nm).
        """
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"PDB file not found at: {pdb_path}")

        box_x = box_y = box_z = 0.0
        z_values = []

        with open(pdb_path, 'r') as pdb:
            for line in pdb:
                if line.startswith("CRYST1"):
                    parts = line.split()
                    box_x = float(parts[1])
                    box_y = float(parts[2])
                    box_z = float(parts[3])
                elif line.startswith("ATOM") or line.startswith("HETATM"):
                    # z coordinate: columns 47-54 of a PDB atom line
                    z_values.append(float(line[46:54].strip()))

        if not z_values:
            raise ValueError(f"No coordinate data found in {pdb_path}")

        slab_thickness = round(abs(max(z_values) - min(z_values)), 3)
        empty_space_nm = (box_z - slab_thickness) / 2 / 10

        return box_x, box_y, box_z, empty_space_nm, slab_thickness

    def write_renamed_pdb(self, input_pdb):
        """
        Copy the input PDB to RENAMED_PDB with the residue name and the atom names
        used by the force field (SI, OS, OD, HD) and the element of every atom.
        """
        box_x, box_y, box_z, empty_space_nm, slab_thickness = self.read_pdb_dimensions(input_pdb)

        with open(input_pdb, 'r') as pdb_in:
            with open(RENAMED_PDB, 'w+') as pdb_out:
                # CRYST1 in the fixed columns of the PDB format (box lengths in columns 7-33)
                pdb_out.writelines(
                    'CRYST1{:9.3f}{:9.3f}{:9.3f}  90.00  90.00  90.00 P 1           1 \n'.format(
                        box_x, box_y, box_z))
                for line in pdb_in:
                    columns = line.split()
                    if columns[0] == 'ATOM':
                        columns[4] = '1'
                        columns[3] = self.resname
                        chain_id = 'X'
                        # Every oxygen that is not a hydroxyl oxygen (OD) is a bridging oxygen. In the
                        # relaxed slab they are called OS..., in the CHARMM-GUI unit cell O3, O12, ...
                        if columns[2][0] == "S":
                            columns[2] = self.silicon_name
                            element = 'Si'
                        elif columns[2][:2] == "OD":
                            columns[2] = self.hydroxyl_oxygen_name
                            element = 'O'
                        elif columns[2][0] == "O":
                            columns[2] = self.bridging_oxygen_name
                            element = 'O'
                        elif columns[2][:2] == "HD":
                            columns[2] = self.hydroxyl_hydrogen_name
                            element = 'H'
                        else:
                            raise ValueError(f"unknown silica atom name {columns[2]} in {input_pdb}")

                        pdb_out.writelines(
                            "{:>s}{:>7s}{:>5s}{:>4s}{:>2s}{:>4s}{:>12s}{:>8s}{:>8s}{:>6s}{:>6s}{:>12s}\n".format(
                                columns[0], columns[1], columns[2], columns[3], chain_id, columns[4], columns[5],
                                columns[6], columns[7], columns[8], columns[9], element, ))
                pdb_out.writelines("END\n")

    def load_slab(self, input_pdb):
        """
        Load a silica PDB as an mbuild compound: place it in a periodic box with the empty space
        in z, repeat it x_tiles x y_tiles times and add the Si-O and O-H bonds.
        """
        self.write_renamed_pdb(input_pdb)
        box_x, box_y, box_z, empty_space_nm, slab_thickness = self.read_pdb_dimensions(input_pdb)

        slab = mb.Compound()
        slab.name = self.resname
        slab.resname = self.resname
        mb.load(RENAMED_PDB, compound=slab)
        slab.periodicity = (True, True, True)

        # Move the slab down by the empty space of the input file, set the new box (nm) and
        # move the slab up to the middle of the new empty space.
        slab.translate(np.array([0.0, 0.0, -empty_space_nm]))
        slab.box = mb.Box([box_x / 10, box_y / 10, (slab_thickness / 10) + self.total_vacuum_nm])
        slab.translate(np.array([0.0, 0.0, (self.total_vacuum_nm / 2) + (slab_thickness / 10 / 20)]))

        slab = TiledCompound(slab, n_tiles=(self.x_tiles, self.y_tiles, 1))

        slab.generate_bonds(name_a=self.silicon_name, name_b=self.bridging_oxygen_name,
                            dmin=0.0, dmax=SI_O_BOND_CUTOFF)
        slab.generate_bonds(name_a=self.silicon_name, name_b=self.hydroxyl_oxygen_name,
                            dmin=0.0, dmax=SI_O_BOND_CUTOFF)
        slab.generate_bonds(name_a=self.hydroxyl_oxygen_name, name_b=self.hydroxyl_hydrogen_name,
                            dmin=0.0, dmax=O_H_BOND_CUTOFF)
        return slab

    def make_slab(self, mode='slab', user_surface=None):
        """
        Return the silica slab for the chosen surface mode ('slab', 'cell' or 'user').
        For 'user', user_surface is a .pdb or .gro file of the slab, for example the confout.gro
        of a surface relaxation made with apadana.
        """
        if mode == 'slab':
            return self.load_slab(self.relaxed_slab_pdb)
        if mode == 'cell':
            return self.load_slab(self.unit_cell_pdb)
        if mode == 'user':
            if user_surface is None:
                raise ValueError("surface mode 'user' needs user_surface, the path to a .pdb or .gro file of the slab")
            if str(user_surface).endswith('.gro'):
                gro_to_pdb(user_surface, USER_SURFACE_PDB)
                return self.load_slab(USER_SURFACE_PDB)
            return self.load_slab(user_surface)
        raise ValueError(f"unknown surface mode {mode!r}; use 'slab', 'cell' or 'user'")


def default_tiles(mode):
    """
    Repeats of the silica structure in x and y when none are given: 8 x 9 unit cells in 'cell'
    mode (7.86 x 7.66 nm, the size of the relaxed slab of the package), and 1 x 1 for a whole slab.
    """
    if mode == 'cell':
        return DEFAULT_CELL_TILES
    return 1, 1


def gro_to_pdb(gro_file, pdb_file):
    """
    Write a .gro file (nm) as a PDB file (Angstrom) with a CRYST1 line, in the column layout
    that SilicaSlab reads. Used for slabs that come from a GROMACS run.
    """
    with open(gro_file) as gro:
        lines = gro.read().splitlines()
    n_atoms = int(lines[1])
    box = [float(value) * 10 for value in lines[2 + n_atoms].split()[:3]]

    with open(pdb_file, 'w') as pdb:
        pdb.write('CRYST1{:9.4f}{:9.4f}{:9.4f}  90.00  90.00  90.00 P 1           1\n'.format(*box))
        for n, line in enumerate(lines[2:2 + n_atoms], start=1):
            # fixed columns of the .gro format
            residue_number = int(line[0:5])
            residue_name = line[5:10].strip()
            atom_name = line[10:15].strip()
            x, y, z = float(line[20:28]) * 10, float(line[28:36]) * 10, float(line[36:44]) * 10
            pdb.write('ATOM  {:>5d} {:<4s} {:<4s}{:>5d}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00\n'.format(
                n % 100000, atom_name, residue_name, residue_number, x, y, z))
        pdb.write('END\n')
