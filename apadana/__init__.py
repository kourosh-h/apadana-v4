"""
apadana: builds alkylsilane self-assembled monolayers (SAMs) on a hydroxylated
alpha-quartz silica slab and writes the GROMACS input files.

Main function: build_sam_system (see builder.py).
"""
from .builder import build_sam_system, functionalizer
from .molecules import Alkylsilane, AlkylsilaneMolecule
from .surface import SilicaSlab
from .visualize import plot_system

__version__ = '4.0'
