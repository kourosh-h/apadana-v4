"""
Build one SAM system with apadana.

Change the values below and run:  python build_example.py
"""
from apadana import build_sam_system

build_sam_system(
    chain_length=8,          # number of carbon atoms in the alkyl chain, L
    coverage=2.40,           # grafting density on each face, molecules/nm^2
    surface_mode='slab',     # 'slab' (relaxed slab of the package), 'cell' (new slab from unit
                             # cells, with coverage=0) or 'user' (your relaxed slab, see user_surface)
    user_surface=None,       # .gro or .pdb file of your slab for surface_mode='user'
    x_tiles=None,            # repeats of the silica structure in x (None: 8 in 'cell' mode, else 1)
    y_tiles=None,            # repeats of the silica structure in y (None: 9 in 'cell' mode, else 1)
    vacuum_space_nm=6,       # empty space above and below the slab, nm
    restrained_layer_nm=1,   # width of the restrained middle layer of the slab, nm
    surface_resname='SUR',   # residue name of the silica
    sam_resname='OTS',       # residue name of the alkylsilanes
    system_name='func',      # base name of the output files
    output_folder='L8_G2.40',
    hpc={'account': 'my_account', 'partition': 'batch'},  # SLURM settings (see HPC_DEFAULTS in workflow.py)
    random_seed=None,        # an integer gives the same grafting sites every time
)
