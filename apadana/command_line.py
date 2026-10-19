"""
Command-line tool of apadana.

    apadana build    build one system, or a series of systems at the same time
    apadana sheet    draw the overview sheet of a built system again

The same commands work as  python -m apadana ...  Run  apadana --help  or  apadana build --help
for all options. Options that are not given keep the defaults of build_sam_system.
"""
import argparse
import sys
from pathlib import Path

from . import __version__
from .builder import build_sam_system
from .parallel import build_many_systems
from .surface import default_tiles
from .visualize import plot_system
from .workflow import HPC_DEFAULTS

BUILD_EXAMPLES = """examples:
  apadana build --chain-length 8 --coverage 2.40 --seed 1
      one system, written to L8_G2.40/
  apadana build --chain-length 4 5 6 7 8 --coverage 0.45 2.40 --seed 1 --processes 4
      ten systems (every chain length with every coverage), four at the same time
  apadana build --coverage 0 --surface-mode cell --x-tiles 10 --y-tiles 12
      a new silica slab of 10 x 12 unit cells, with the files for its relaxation
  apadana build --chain-length 8 --coverage 2.40 --surface-mode user --user-surface confout.gro
      chains grafted on a relaxed slab of the user
  apadana build --chain-length 8 --coverage 2.40 --hpc account=my_account --hpc partition=batch
      with the SLURM settings of the job scripts
"""


# ----------------------------------------------------------------------
# Reading the command line
# ----------------------------------------------------------------------

def read_hpc_settings(pairs):
    """SLURM settings of the job scripts from pairs key=value, with the keys of HPC_DEFAULTS."""
    hpc = {}
    for pair in pairs:
        if '=' not in pair:
            raise ValueError(f'--hpc expects key=value, not {pair!r}')
        key, value = pair.split('=', 1)
        if key not in HPC_DEFAULTS:
            raise ValueError(f'unknown --hpc key {key!r}; the keys are: {", ".join(HPC_DEFAULTS)}')
        # numbers stay numbers (nodes, ntasks_per_node, cpus_per_task)
        if isinstance(HPC_DEFAULTS[key], int):
            value = int(value)
        hpc[key] = value
    return hpc


def coverage_number(coverage_text):
    """
    The coverage as a number, as it would be written in Python: a whole number when it is typed
    without a decimal point (0), else a float (2.40). The files then read the same as those of a
    build_sam_system call with the same value.
    """
    if coverage_text.isdigit():
        return int(coverage_text)
    return float(coverage_text)


def folder_name(chain_length, coverage_text, surface_mode, x_tiles, y_tiles):
    """
    Default folder of a system: L<chain length>_G<coverage as typed>, for example L8_G2.40,
    or slab_<x>x<y> for a new slab built from unit cells.
    """
    if surface_mode == 'cell' and float(coverage_text) == 0:
        default_x, default_y = default_tiles('cell')
        if x_tiles is None:
            x_tiles = default_x
        if y_tiles is None:
            y_tiles = default_y
        return f'slab_{x_tiles}x{y_tiles}'
    return f'L{chain_length}_G{coverage_text}'


def common_build_settings(arguments):
    """
    The build_sam_system arguments that are the same for all systems of one command. Only the
    options that were given are passed on, so that all other values keep their defaults.
    """
    options = {
        'surface_mode': arguments.surface_mode,
        'user_surface': arguments.user_surface,
        'x_tiles': arguments.x_tiles,
        'y_tiles': arguments.y_tiles,
        'vacuum_space_nm': arguments.vacuum,
        'restrained_layer_nm': arguments.restrained_layer,
        'surface_resname': arguments.surface_resname,
        'sam_resname': arguments.sam_resname,
        'system_name': arguments.system_name,
        'random_seed': arguments.seed,
    }
    settings = {}
    for key, value in options.items():
        if value is not None:
            settings[key] = value
    if arguments.hpc:
        settings['hpc'] = read_hpc_settings(arguments.hpc)
    if arguments.no_sheet:
        settings['make_plot'] = False
    return settings


# ----------------------------------------------------------------------
# apadana build
# ----------------------------------------------------------------------

def run_build(arguments):
    """Build one system, or every combination of the given chain lengths and coverages."""
    for coverage_text in arguments.coverage:
        if float(coverage_text) < 0:
            raise ValueError(f'the coverage cannot be negative ({coverage_text})')
    chain_lengths = arguments.chain_length
    if chain_lengths is None:
        # a slab without chains needs no chain length
        for coverage_text in arguments.coverage:
            if float(coverage_text) > 0:
                raise ValueError('give --chain-length for a coverage above 0')
        chain_lengths = [0]
    if arguments.surface_mode == 'user' and arguments.user_surface is None:
        raise ValueError("--surface-mode user needs --user-surface, the .gro or .pdb file of the slab")

    settings = common_build_settings(arguments)
    surface_mode = settings.get('surface_mode', 'slab')

    systems = []
    for chain_length in chain_lengths:
        for coverage_text in arguments.coverage:
            name = folder_name(chain_length, coverage_text, surface_mode, arguments.x_tiles, arguments.y_tiles)
            systems.append({'chain_length': chain_length, 'coverage': coverage_number(coverage_text), 'name': name})

    # one system: --output is its folder
    if len(systems) == 1:
        system = systems[0]
        output_folder = arguments.output or system['name']
        build_sam_system(chain_length=system['chain_length'], coverage=system['coverage'],
                         output_folder=output_folder, **settings)
        return 0

    # several systems: --output is the folder in which their folders are made
    parent = Path(arguments.output or '.')
    series = []
    for system in systems:
        series.append({'chain_length': system['chain_length'], 'coverage': system['coverage'],
                       'output_folder': str(parent / system['name'])})
    results = build_many_systems(series, processes=arguments.processes, **settings)
    for output_folder, finished, minutes in results:
        if not finished:
            return 1
    return 0


# ----------------------------------------------------------------------
# apadana sheet
# ----------------------------------------------------------------------

def read_build_settings(readme_file):
    """The rows of the table 'Build settings' in the README of a built system, as strings."""
    settings = {}
    in_table = False
    with open(readme_file) as readme:
        for line in readme:
            if line.startswith('## '):
                in_table = line.strip() == '## Build settings'
                continue
            columns = [column.strip() for column in line.strip().strip('|').split('|')]
            if in_table and len(columns) == 2 and columns[0] not in ['Setting', '---']:
                settings[columns[0]] = columns[1]
    return settings


def run_sheet(arguments):
    """Draw the overview sheet of a built system again, from its files."""
    folder = Path(arguments.folder)
    readme_file = folder / 'README.md'
    if not readme_file.exists():
        raise ValueError(f'{folder} is not the folder of a system built by apadana (no README.md)')
    settings = read_build_settings(readme_file)
    gro_file = folder / 'system' / f"{settings.get('system_name', 'func')}.gro"
    if not gro_file.exists():
        raise ValueError(f'{gro_file} is missing')

    chain_length = int(settings['chain_length'])
    coverage = float(settings['coverage'])
    # a new slab of unit cells has the relaxation script instead of the dry and wet workflow
    if (folder / 'run_surface_relaxation.sh').exists():
        workflow = 'surface'
    else:
        workflow = 'sam'

    if arguments.output:
        output_files = [arguments.output]
    else:
        output_files = [folder / 'system_overview.pdf', folder / 'system_overview.png']
    for output_file in output_files:
        plot_system(gro_file, output_file, chain_length, coverage,
                    sam_resname=settings.get('sam_resname', 'OTS'),
                    surface_resname=settings.get('surface_resname', 'SUR'),
                    restrained_layer_nm=float(settings.get('restrained_layer_nm', 1)),
                    settings=settings, workflow=workflow)
        print(f'written: {output_file}')
    return 0


# ----------------------------------------------------------------------
# The parser
# ----------------------------------------------------------------------

def make_parser():
    parser = argparse.ArgumentParser(
        prog='apadana',
        description='Builds alkylsilane SAMs on hydroxylated alpha-quartz silica and prepares their '
                    'GROMACS simulations.')
    parser.add_argument('--version', action='version', version=f'apadana {__version__}')
    commands = parser.add_subparsers(dest='command', title='commands')

    build = commands.add_parser(
        'build', help='build one system, or a series of systems at the same time',
        description='Build SAM systems. Several values of --chain-length or --coverage build every\n'
                    'combination, one system per process.',
        epilog=BUILD_EXAMPLES, formatter_class=argparse.RawDescriptionHelpFormatter)
    build.add_argument('--chain-length', type=int, nargs='+', metavar='L',
                       help='number of carbon atoms in the alkyl chain (one or more values)')
    build.add_argument('--coverage', nargs='+', required=True, metavar='GAMMA',
                       help='coverage on each face in molecules/nm^2 (one or more values); 0 for a slab without chains')
    build.add_argument('--seed', type=int, metavar='N', help='random seed: the same seed gives the same binding locations')
    build.add_argument('--output', metavar='FOLDER', help='folder of the system (one system), or the folder in which the '
                                        'folders of a series are made (default: L<L>_G<coverage> here)')
    build.add_argument('--processes', type=int, metavar='N',
                       help='number of systems built at the same time (default: one per CPU core)')
    build.add_argument('--surface-mode', choices=['slab', 'cell', 'user'],
                       help="silica slab: relaxed slab of the package ('slab', default), new slab from unit "
                            "cells ('cell') or a slab of the user ('user')")
    build.add_argument('--user-surface', metavar='FILE', help=".gro or .pdb file of the slab for --surface-mode user")
    build.add_argument('--x-tiles', type=int, metavar='N', help='repeats of the silica structure in x (default 8 for cell, else 1)')
    build.add_argument('--y-tiles', type=int, metavar='N', help='repeats of the silica structure in y (default 9 for cell, else 1)')
    build.add_argument('--vacuum', type=float, metavar='NM', help='empty space above and below the slab (default 6 nm)')
    build.add_argument('--restrained-layer', type=float, metavar='NM',
                       help='width of the restrained centermost layer of the slab (default 1 nm)')
    build.add_argument('--surface-resname', metavar='NAME', help='residue name of the silica (default SUR)')
    build.add_argument('--sam-resname', metavar='NAME', help='residue name of the alkylsilanes (default OTS)')
    build.add_argument('--system-name', metavar='NAME', help='base name of the structure and topology files (default func)')
    build.add_argument('--hpc', action='append', metavar='KEY=VALUE',
                       help=f'SLURM setting of the job scripts, repeated for several; keys: {", ".join(HPC_DEFAULTS)}')
    build.add_argument('--no-sheet', action='store_true', help='do not draw the overview sheet')
    build.set_defaults(function=run_build)

    sheet = commands.add_parser('sheet', help='draw the overview sheet of a built system again',
                                description='Draw the overview sheet of a system built by apadana again, '
                                            'from the files in its folder.')
    sheet.add_argument('folder', help='folder of the system')
    sheet.add_argument('--output', metavar='FILE', help='file to write (.pdf or .png); default: system_overview.pdf and '
                                        '.png in the folder of the system')
    sheet.set_defaults(function=run_sheet)
    return parser


def main(argv=None):
    """Run the command-line tool; returns the exit code (0 when everything worked)."""
    parser = make_parser()
    arguments = parser.parse_args(argv)
    if arguments.command is None:
        parser.print_help()
        return 0
    try:
        return arguments.function(arguments)
    except ValueError as error:
        # a wrong input: a short message instead of a traceback
        print(f'apadana: {error}', file=sys.stderr)
        return 2
