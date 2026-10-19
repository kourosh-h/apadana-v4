"""
apadana: builds alkylsilane self-assembled monolayers (SAMs) on a hydroxylated
alpha-quartz silica slab and writes the GROMACS input files.

Main functions: build_sam_system (see builder.py) and build_many_systems, which builds several
systems at the same time (see parallel.py).
"""
from .builder import build_sam_system, functionalizer
from .molecules import Alkylsilane, AlkylsilaneMolecule
from .parallel import build_many_systems
from .surface import SilicaSlab
from .visualize import plot_system

__version__ = '4.0'
