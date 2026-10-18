"""
Make the images of the README: the banner (a light and a dark version) and the upper parts of
two overview sheets, one of a SAM system and one of a new silica slab.

The systems are built with apadana in work_folder (default: readme_images/ in the current folder);
systems that are already there are not built again. Building all of them takes about 20 minutes.

Run:  python docs/make_readme_images.py [work_folder]
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use('Agg')  # write files only, no window
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from PIL import Image

from apadana import build_sam_system, plot_system
from apadana.visualize import CHAIN_COLOR, SILICA_FILL, STYLE, chain_backbones, read_gro

IMAGES = Path(__file__).parent / 'images'

# Banner: chains of length L = 8 at three coverages (molecules/nm^2), side view in the yz-plane.
BANNER_CHAIN_LENGTH = 8
BANNER_COVERAGES = [0.45, 1.65, 3.90]
BANNER_Z_RANGE = 2.4  # nm above and below the slab center

# Colors of the banner on a light and on a dark page (GitHub light and dark themes).
BANNER_COLORS = {
    'light': {'background': '#FFFFFF', 'text': '#555555', 'chain': CHAIN_COLOR, 'anchor': '#0C3C84',
              'silica_fill': SILICA_FILL, 'silicon': '#E4840C', 'oxygen': '#9C540C'},
    'dark': {'background': '#0D1117', 'text': '#C9D1D9', 'chain': '#5A9BEA', 'anchor': '#A8CBF5',
             'silica_fill': '#3A2410', 'silicon': '#F0A030', 'oxygen': '#B4651A'},
}

# Height (inches from the top of the A4 sheet) of the part of the overview sheets shown in the README.
SHEET_PART_INCHES = 8.45
PAGE_HEIGHT_INCHES = 11.69


def build(folder, **settings):
    """Build a system with apadana unless it is already in folder. Returns its .gro file."""
    gro_file = folder / 'system' / 'func.gro'
    if not gro_file.exists():
        build_sam_system(output_folder=folder, make_plot=False, **settings)
    return gro_file


def draw_banner(gro_files, coverages, output_file, theme):
    """
    Side views (yz-plane) of the silica slab with its SAMs, one panel per coverage. The silica atoms
    are drawn as small dots (silicon and oxygen), the chains as lines through their backbones.
    """
    colors = BANNER_COLORS[theme]
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(9.0, 2.3), facecolor=colors['background'])
        width = 0.29
        for n, (gro_file, coverage) in enumerate(zip(gro_files, coverages)):
            atoms = read_gro(gro_file)
            box_x, box_y, box_z = atoms['box']
            residue_names = np.array(atoms['residue_name'])
            atom_names = np.array(atoms['atom_name'])
            xyz = atoms['xyz']
            is_silica = residue_names == 'SUR'
            center = (xyz[is_silica, 2].max() + xyz[is_silica, 2].min()) / 2
            xyz[:, 2] -= center
            slab_bottom, slab_top = xyz[is_silica, 2].min(), xyz[is_silica, 2].max()
            backbones = chain_backbones(atoms, 'OTS', BANNER_CHAIN_LENGTH)

            ax = fig.add_axes([0.035 + n * (width + 0.035), 0.17, width, 0.8])
            ax.set_facecolor(colors['background'])
            ax.add_patch(Rectangle((0, slab_bottom), box_y, slab_top - slab_bottom, facecolor=colors['silica_fill'],
                                   lw=0))
            for element, color in [('O', colors['oxygen']), ('SI', colors['silicon'])]:
                chosen = is_silica & np.char.startswith(atom_names, element)
                ax.scatter(xyz[chosen, 1], xyz[chosen, 2], s=0.8, color=color, lw=0, rasterized=True)
            for backbone in backbones:
                ax.plot(backbone[:, 1], backbone[:, 2], color=colors['chain'], lw=0.75, solid_capstyle='round')
            anchors = np.array([backbone[0] for backbone in backbones])
            ax.scatter(anchors[:, 1], anchors[:, 2], s=3, color=colors['anchor'], lw=0, zorder=3)
            ax.set_xlim(0, box_y)
            ax.set_ylim(-BANNER_Z_RANGE, BANNER_Z_RANGE)
            ax.set_aspect('equal')
            ax.axis('off')
            fig.text(0.035 + n * (width + 0.035) + width / 2, 0.08,
                     f'$\\Gamma$ = {coverage:.2f} molecules/nm²', ha='center', va='bottom', fontsize=9,
                     color=colors['text'])
        fig.savefig(output_file, dpi=200, facecolor=colors['background'])
        plt.close(fig)


def sheet_part(gro_file, output_file, chain_length, coverage, settings, workflow):
    """Overview sheet of a system, cut below the maps of the two faces, as a compact .png."""
    full_sheet = output_file.with_name('full_' + output_file.name)
    plot_system(gro_file, full_sheet, chain_length, coverage, settings=settings, workflow=workflow)
    image = Image.open(full_sheet)
    width, height = image.size
    part = image.crop((0, 0, width, int(height * SHEET_PART_INCHES / PAGE_HEIGHT_INCHES)))
    part = part.convert('P', palette=Image.ADAPTIVE, colors=256)
    part.save(output_file, optimize=True)
    full_sheet.unlink()


def main(work_folder):
    work_folder = Path(work_folder)
    IMAGES.mkdir(exist_ok=True)

    # banner
    gro_files = []
    for coverage in BANNER_COVERAGES:
        folder = work_folder / f'L{BANNER_CHAIN_LENGTH}_G{coverage:.2f}'
        gro_files.append(build(folder, chain_length=BANNER_CHAIN_LENGTH, coverage=coverage, random_seed=1))
    for theme in ['light', 'dark']:
        draw_banner(gro_files, BANNER_COVERAGES, IMAGES / f'banner_{theme}.png', theme)

    # overview sheet of the test case
    test_case = Path(__file__).parent.parent / 'examples' / 'test_case_L6_G1.65' / 'system' / 'func.gro'
    settings = {'surface_mode': 'slab', 'x_tiles': 1, 'y_tiles': 1, 'random_seed': 1, 'system_name': 'func'}
    sheet_part(test_case, IMAGES / 'overview_L6_G1.65.png', 6, 1.65, settings, 'sam')

    # overview sheet of a new slab of 8 x 9 unit cells
    slab = build(work_folder / 'new_slab_8x9', chain_length=0, coverage=0, surface_mode='cell')
    settings = {'surface_mode': 'cell', 'x_tiles': 8, 'y_tiles': 9, 'system_name': 'func'}
    sheet_part(slab, IMAGES / 'overview_new_slab_8x9.png', 0, 0, settings, 'surface')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'readme_images')
