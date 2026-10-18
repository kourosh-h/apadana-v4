"""
Building blocks of the alkylsilane molecules.

An alkylsilane of chain length L (number of carbon atoms) is built from
  - the silane head SiOH2 (atoms SI1, O1, O2, H1, H2),
  - L - 1 CH2 groups (C1 ... C_{L-1}),
  - a terminal CH3 group (C_L).
The carbon atoms are labeled C1 ... C_L from the silicon to the chain end, so the atom names of
one molecule are SI1, O1, O2, H1, H2, C1, ..., C_L and its hydrogens.

All lengths are in nm (mbuild units).
"""
import mbuild as mb
import numpy as np

# Half of the angle between a tetrahedral angle (109.5 deg) and a straight line, in radians.
# Used to place hydrogens and ports around a tetrahedral atom.
HALF_TETRAHEDRAL_OFFSET = 0.5 * (180 - 109.5) * np.pi / 180


# ----------------------------------------------------------------------
# Chemical groups
# ----------------------------------------------------------------------

class Hydroxyl(mb.Compound):
    """OH group with one port ('up'); used to cap the silane head of a free molecule."""

    def __init__(self, i=1, j=1):
        super(Hydroxyl, self).__init__()
        oxygen = mb.Particle(name=f'O{i}', element='O', pos=[0, 0, 0])
        self.add(oxygen, label='OH')
        hydrogen = mb.Particle(name=f'H{j}', element='H', pos=[0.120, 0, 0])
        self.add(hydrogen, label='HO')
        self.add_bond((self[0], self['HO'][0]))
        port = mb.Port(anchor=self[0], orientation=[0, 1, 0], )
        self.add(port, label='up')
        self['up'].translate([0, 0.06, 0])


class SiOH2(mb.Compound):
    """
    Silane head group: Si with two hydroxyls (SI1, O1, H1, O2, H2).
    Port 'up' bonds to the first carbon, port 'down' to the oxygen of the silica surface.
    """

    def __init__(self):
        super(SiOH2, self).__init__()
        self.add(mb.Particle(name='SI1', element='Si', charge=0.78, pos=[0, 0, 0]), label='SI[$]')
        self.add(mb.Particle(name='O1', element='O', pos=[0.180, 0, 0]), label='OH[$]')
        self.add(mb.Particle(name='O2', element='O', pos=[-0.180, 0, 0]), label='OH[$]')
        self.add(mb.Particle(name='H1', element='H', pos=[0.192, 0, -0.11]), label='HO[$]')
        self.add(mb.Particle(name='H2', element='H', pos=[-0.192, 0, -0.11]), label='HO[$]')

        self.add_bond((self[0], self['OH'][0]))
        self.add_bond((self[0], self['OH'][1]))
        self.add_bond((self['OH'][0], self['HO'][0]))
        self.add_bond((self['OH'][1], self['HO'][1]))

        theta = HALF_TETRAHEDRAL_OFFSET
        self['OH'][0].rotate(theta, around=[0, 1, 0.4])
        self['OH'][1].rotate(-theta, around=[0, 1, 0.4])
        self['HO'][0].rotate(theta, around=[0, 1, 0.4])
        self['HO'][1].rotate(-theta, around=[0, 1, 0.4])

        self.add(mb.Port(anchor=self[0]), label='up')
        self['up'].translate([0, 0.098, 0])
        self['up'].rotate(theta, around=[1, 0, 0])

        self.add(mb.Port(anchor=self[0]), label='down')
        self['down'].translate([0, -0.095, 0])
        self['down'].rotate(np.pi, [0, 1, 0])


class CH2(mb.Compound):
    """Methylene group C{i} with hydrogens H{j}, H{j+1}; ports 'up' and 'down' along the chain (C-C 0.154 nm)."""

    def __init__(self, i=1, j=1):
        super(CH2, self).__init__()
        self.add(mb.Particle(name=f'C{i}', element='C', pos=[0, 0, 0]), label='CT[$]')
        self.add(mb.Particle(name=f'H{j}', element='H', pos=[-0.109, 0, 0.0]), label='HC[$]')
        self.add(mb.Particle(name=f'H{j + 1}', element='H', pos=[0.109, 0, 0.0]), label='HC[$]')

        theta = HALF_TETRAHEDRAL_OFFSET
        self['HC'][0].rotate(-theta, around=[0, 1, 0])
        self['HC'][1].rotate(theta, around=[0, 1, 0])

        self.add_bond((self['CT'][0], self['HC'][0]))
        self.add_bond((self['CT'][0], self['HC'][1]))

        self.add(mb.Port(anchor=self[0]), label='up')
        self['up'].translate([0, 0.154 / 2, 0])
        self['up'].rotate(theta, around=[1, 0, 0])

        self.add(mb.Port(anchor=self[0]), label='down')
        self['down'].translate([0, -0.154 / 2, 0])
        self['down'].rotate(np.pi, [0, 1, 0])
        self['down'].rotate(-theta, around=[1, 0, 0])


class CF2(mb.Compound):
    """Difluoromethylene group, the fluorinated version of CH2."""

    def __init__(self):
        super(CF2, self).__init__()
        self.add(mb.Particle(name='CF', element='C', pos=[0, 0, 0]), label='CT[$]')
        self.add(mb.Particle(name='F', element='F', pos=[-0.109, 0, 0.0]), label='FC[$]')
        self.add(mb.Particle(name='F', element='F', pos=[0.109, 0, 0.0]), label='FC[$]')

        theta = HALF_TETRAHEDRAL_OFFSET
        self['FC'][0].rotate(-theta, around=[0, 1, 0])
        self['FC'][1].rotate(theta, around=[0, 1, 0])

        self.add_bond((self['CF'][0], self['FC'][0]))
        self.add_bond((self['CF'][0], self['FC'][1]))

        self.add(mb.Port(anchor=self[0]), label='up')
        self['up'].translate([0, 0.154 / 2, 0])
        self['up'].rotate(theta, around=[1, 0, 0])

        self.add(mb.Port(anchor=self[0]), label='down')
        self['down'].translate([0, -0.154 / 2, 0])
        self['down'].rotate(np.pi, [0, 1, 0])
        self['down'].rotate(-theta, around=[1, 0, 0])


class CH3(mb.Compound):
    """Terminal methyl group C{i} with hydrogens H{j}, H{j+1}, H{j+2}; one port 'down'."""

    def __init__(self, i=1, j=1):
        super(CH3, self).__init__()
        theta = HALF_TETRAHEDRAL_OFFSET

        carbon = mb.Particle(name=f'C{i}', element='C', pos=[0, 0, 0])
        self.add(carbon, label='CT[$]')

        self.add(mb.Particle(name=f'H{j}', element='H', pos=[0.1, 0, -0.07]), label='HC[$]')
        self.add(mb.Particle(name=f'H{j + 1}', element='H', pos=[-0.1, 0, -0.07]), label='HC[$]')
        self.add(mb.Particle(name=f'H{j + 2}', element='H', pos=[0, 0.1, 0.07]), label='HC[$]')

        self.add_bond((self[0], self['HC'][0]))
        self.add_bond((self[0], self['HC'][1]))
        self.add_bond((self[0], self['HC'][2]))

        self.add(mb.Port(anchor=self[0]), label='down')
        self['down'].translate([0, -0.154 / 2, 0])
        self['down'].rotate(np.pi, [0, 1, 0])
        self['down'].rotate(-theta, around=[1, 0, 0])


class CF3(mb.Compound):
    """Trifluoromethyl group, the fluorinated version of CH3."""

    def __init__(self):
        super(CF3, self).__init__()
        theta = HALF_TETRAHEDRAL_OFFSET

        carbon = mb.Particle(name='CF', element='C', pos=[0, 0, 0])
        self.add(carbon, label='CF[$]')

        self.add(mb.Particle(name='F', element='F', pos=[0.1, 0, -0.07]), label='FC[$]')
        self.add(mb.Particle(name='F', element='F', pos=[-0.1, 0, -0.07]), label='FC[$]')
        self.add(mb.Particle(name='F', element='F', pos=[0, 0.1, 0.07]), label='FC[$]')

        self.add_bond((self[0], self['FC'][0]))
        self.add_bond((self[0], self['FC'][1]))
        self.add_bond((self[0], self['FC'][2]))

        self.add(mb.Port(anchor=self[0]), label='down')
        self['down'].translate([0, -0.154 / 2, 0])
        self['down'].rotate(np.pi, [0, 1, 0])
        self['down'].rotate(-theta, around=[1, 0, 0])


# ----------------------------------------------------------------------
# Chains
# ----------------------------------------------------------------------

class AlkylChain(mb.Compound):
    """
    All-trans chain of n_groups CH2 groups. Carbons are numbered from i, hydrogens from j.
    Each new group is rotated by 180 deg so that the chain forms a zigzag.
    """

    def __init__(self, n_groups, i=1, j=4):
        super(AlkylChain, self).__init__()
        half_turn = 180 * np.pi / 180

        last_group = CH2(i, j)
        self.add(last_group)
        self.add(last_group['down'], 'down', containment=False)

        for n in range(n_groups - 1):
            i += 1
            j += 2
            new_group = CH2(i, j)
            last_group.rotate(half_turn, around=[0, 1, 0])
            mb.force_overlap(move_this=new_group,
                             from_positions=new_group['down'],
                             to_positions=last_group['up'])
            new_group.rotate(half_turn, around=[0, 1, 0])
            self.add(new_group)
            last_group = new_group
        last_group.rotate(half_turn, around=[0, 1, 0])

        self.add(last_group['up'], 'up', containment=False)


class FluoroalkylChain(mb.Compound):
    """Chain with one CH2 group next to the head followed by CF2 groups."""

    def __init__(self, n_groups):
        super(FluoroalkylChain, self).__init__()
        half_turn = 180 * np.pi / 180

        last_group = CH2()
        self.add(last_group)
        self.add(last_group['down'], 'down', containment=False)

        for n in range(n_groups - 1):
            if n < 1:
                new_group = CH2()
            else:
                new_group = CF2()
            last_group.rotate(half_turn, around=[0, 1, 0])
            mb.force_overlap(move_this=new_group,
                             from_positions=new_group['down'],
                             to_positions=last_group['up'])
            new_group.rotate(half_turn, around=[0, 1, 0])
            self.add(new_group)
            last_group = new_group
        last_group.rotate(half_turn, around=[0, 1, 0])

        self.add(last_group['up'], 'up', containment=False)


# ----------------------------------------------------------------------
# Alkylsilane molecule
# ----------------------------------------------------------------------

class AlkylsilaneMolecule(mb.Compound):
    """
    Alkylsilane with n_ch2 CH2 groups and a terminal CH3, i.e. chain length L = n_ch2 + 1.

    capped=True : free molecule, the head has a third OH group (O3, H3).
    capped=False: molecule for grafting; the head keeps an open port 'down' that is
                  bonded to an oxygen of the silica surface.
    """

    def __init__(self, n_ch2=6, capped=True):
        super(AlkylsilaneMolecule, self).__init__()
        head = SiOH2()
        head.name = 'SIT'
        tail_carbon = n_ch2 + 1

        if capped:
            backbone = AlkylChain(n_groups=n_ch2, i=1, j=4)
            backbone.name = 'BAK'
            tail = CH3(tail_carbon, (n_ch2 * 2) + 1 + 3)
            tail.name = 'TIL'
            cap = Hydroxyl(i=3, j=3)
            cap.name = 'CAP'
            self.add(head)
            mb.force_overlap(move_this=cap, from_positions=cap['up'], to_positions=head['down'])
            self.add(cap)
            mb.force_overlap(move_this=backbone, from_positions=backbone['down'], to_positions=head['up'])
            self.add(backbone)
            mb.force_overlap(move_this=tail, from_positions=tail['down'], to_positions=backbone['up'])
            self.add(tail)
        else:
            backbone = AlkylChain(n_groups=n_ch2, i=1, j=3)
            backbone.name = 'BAK'
            tail = CH3(tail_carbon, (n_ch2 * 2) + 1 + 2)
            tail.name = 'TIL'
            self.add(head)
            mb.force_overlap(move_this=backbone, from_positions=backbone['down'], to_positions=head['up'])
            self.add(backbone)
            mb.force_overlap(move_this=tail, from_positions=tail['down'], to_positions=backbone['up'])
            self.add(tail)
            self.add(head['down'], 'down', containment=False)


class Alkylsilane:
    """
    Makes alkylsilane molecules of one chain length for grafting on the silica slab.

    n_ch2   : number of CH2 groups, chain length L = n_ch2 + 1
    resname : residue name of the molecules in the output files (OTS)
    """

    def __init__(self, n_ch2=6, resname='OTS'):
        self.n_ch2 = n_ch2
        self.resname = resname

    def capped_molecule(self):
        """Free alkylsilane with three OH groups on the silicon."""
        molecule = AlkylsilaneMolecule(self.n_ch2, capped=True)
        molecule.name = self.resname
        molecule.resname = self.resname
        return molecule

    def for_top_face(self):
        """Molecule for grafting on the top face of the slab (pointing up, +z)."""
        molecule = AlkylsilaneMolecule(self.n_ch2, capped=False)
        molecule.name = self.resname
        molecule.resname = self.resname
        return molecule

    def for_bottom_face(self):
        """Molecule for grafting on the bottom face of the slab (turned by 180 deg, pointing down)."""
        molecule = self.for_top_face()
        molecule.spin(np.pi, [1, 0, 0])
        return molecule

    def save_files(self, molecule, file_name):
        """Save a molecule as .pdb, .mol2 and .xyz."""
        molecule.save(f'{file_name}.pdb', overwrite=True, residues=self.resname)
        molecule.save(f'{file_name}.mol2', overwrite=True, residues=self.resname)
        molecule.save(f'{file_name}.xyz', overwrite=True, residues=self.resname)
