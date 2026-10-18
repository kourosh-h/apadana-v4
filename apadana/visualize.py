"""
Overview sheet of a built system: one A4 page (system_overview.pdf, and a .png for the README).

Lengths are in nm, and z is measured from the z-axial geometrical center of the silica slab
(the GROMACS box has its origin at a corner, so the z of the slab center is subtracted from all
z coordinates). The silica is drawn in orange and the alkylsilane chains in blue.

  (a) side view, xz-plane      (b) side view, yz-plane      (c) mass density along z
  (d) binding locations,       (e) binding locations,       (f) nearest-neighbor distance
      top face                     bottom face                  of the chains on each face
  table of the system data                  table of the simulation workflow

The two faces are told apart by the markers and lines: top face filled markers and solid lines,
bottom face open markers and dashed lines. For a slab without chains, panel (f) shows the atomic
layers of the slab instead.
"""
import matplotlib

matplotlib.use('Agg')  # write files only, no window
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch, Rectangle

from .workflow import (DRY_PRODUCTION, DRY_RELAXATION, SURFACE_RELAXATION, WET_PRODUCTION, WET_RELAXATION,
                       stage_description)

# Atomic masses (g/mol), chosen from the first letters of the atom name.
MASSES = {'Si': 28.0855, 'O': 15.9994, 'C': 12.011, 'H': 1.008}

# Conversion of a mass density from g/mol per nm^3 to kg/m^3.
GRAM_PER_MOL_PER_NM3_TO_KG_PER_M3 = 1.66054

# Colors: silica in orange and chains in blue, as in the snapshots of the systems.
SILICA_FILL = '#F8DCB4'
SILICA_LINE = '#C9690C'
SILICA_TEXT = '#8A4A0C'
CHAIN_COLOR = '#1257B0'
ANCHOR_COLOR = '#0C3C84'
HYDROXYL_COLOR = '#B5B5B5'
TEXT_COLOR = '#000000'
TEXT_SOFT = '#555555'
RULE_COLOR = '#A0A0A0'

# Kinds of silica atoms in the layer diagram of a slab without chains: names in the .gro file,
# label and color.
SILICA_LAYERS = [(['SI'], 'Si', '#FF7F0E'),
                 (['OS'], 'O\nbridging', '#D62728'),
                 (['OD', 'OP'], 'O\nhydroxyl', '#9467BD'),
                 (['HD'], 'H\nhydroxyl', '#7F7F7F')]

# Page: A4 portrait, in inches.
PAGE_WIDTH = 8.27
PAGE_HEIGHT = 11.69

# Text sizes (pt).
FS = 8.0
FS_SMALL = 7.0
FS_PANEL = 9.0
FS_TITLE = 14.0

STYLE = {
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold', 'mathtext.default': 'it',
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.unicode_minus': False,
    'font.size': FS, 'axes.labelsize': FS, 'axes.titlesize': FS, 'xtick.labelsize': FS,
    'ytick.labelsize': FS, 'legend.fontsize': FS_SMALL, 'text.color': TEXT_COLOR, 'axes.labelcolor': TEXT_COLOR,
    'axes.edgecolor': TEXT_COLOR, 'xtick.color': TEXT_COLOR, 'ytick.color': TEXT_COLOR, 'axes.linewidth': 1.0,
    'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.major.size': 4.0, 'ytick.major.size': 4.0,
    'xtick.minor.size': 2.0, 'ytick.minor.size': 2.0, 'xtick.major.width': 0.9, 'ytick.major.width': 0.9,
    'xtick.minor.width': 0.6, 'ytick.minor.width': 0.6, 'xtick.minor.visible': True, 'ytick.minor.visible': True,
    'xtick.major.pad': 3.0, 'ytick.major.pad': 3.0, 'axes.labelpad': 2.5, 'axes.titlepad': 4.0,
    'legend.frameon': False, 'legend.handlelength': 1.8, 'legend.handletextpad': 0.4,
    'legend.labelspacing': 0.28, 'savefig.dpi': 200, 'figure.facecolor': 'white', 'axes.facecolor': 'white',
}


# ----------------------------------------------------------------------
# Reading the structure
# ----------------------------------------------------------------------

def element_of(atom_name):
    """Element of an atom from its name (SI, SI1 -> Si; O... -> O; C... -> C; H... -> H)."""
    if atom_name.startswith('SI'):
        return 'Si'
    return atom_name[0]


def read_gro(gro_file):
    """
    Read a .gro file.

    Returns a dictionary with lists 'residue_number', 'residue_name', 'atom_name' and an
    (n_atoms, 3) array 'xyz' in nm, and the box lengths 'box' (nm).
    """
    with open(gro_file) as gro:
        lines = gro.read().splitlines()
    n_atoms = int(lines[1])
    atoms = {'residue_number': [], 'residue_name': [], 'atom_name': [], 'xyz': np.zeros((n_atoms, 3))}
    for n, line in enumerate(lines[2:2 + n_atoms]):
        # fixed columns of the .gro format
        atoms['residue_number'].append(int(line[0:5]))
        atoms['residue_name'].append(line[5:10].strip())
        atoms['atom_name'].append(line[10:15].strip())
        atoms['xyz'][n] = [float(line[20:28]), float(line[28:36]), float(line[36:44])]
    atoms['box'] = [float(value) for value in lines[2 + n_atoms].split()[:3]]
    return atoms


def chain_backbones(atoms, sam_resname, chain_length):
    """
    Backbone of every grafted chain: an array of the positions of SI1, C1, ..., C_L.
    The atoms of one chain follow each other in the file, and every chain starts with SI1.
    Returns a list with one (L + 1, 3) array per chain.
    """
    backbone_names = ['SI1'] + [f'C{i}' for i in range(1, chain_length + 1)]
    chains = []
    for n, residue_name in enumerate(atoms['residue_name']):
        name = atoms['atom_name'][n]
        if residue_name != sam_resname or name not in backbone_names:
            continue
        if name == 'SI1':
            chains.append({})
        chains[-1][name] = atoms['xyz'][n]

    backbones = []
    for chain in chains:
        backbones.append(np.array([chain[name] for name in backbone_names]))
    return backbones


# ----------------------------------------------------------------------
# Quantities computed from the structure
# ----------------------------------------------------------------------

def mass_density_profile(z_values, masses, area, z_edges):
    """Mass density (kg/m^3) along z of the given atoms, in slices between z_edges (nm)."""
    mass_per_bin, edges = np.histogram(z_values, bins=z_edges, weights=masses)
    bin_volume = area * (z_edges[1] - z_edges[0])
    return mass_per_bin / bin_volume * GRAM_PER_MOL_PER_NM3_TO_KG_PER_M3


def nearest_neighbor_distances(xy, box_x, box_y):
    """
    Distance (nm) from every point to its nearest neighbor in the xy-plane, with the periodic
    images of the box (minimum-image convention).
    """
    distances = []
    for n in range(len(xy)):
        dx = xy[:, 0] - xy[n, 0]
        dy = xy[:, 1] - xy[n, 1]
        dx = dx - box_x * np.round(dx / box_x)
        dy = dy - box_y * np.round(dy / box_y)
        r = np.sqrt(dx ** 2 + dy ** 2)
        r[n] = np.inf
        distances.append(r.min())
    return np.array(distances)


# ----------------------------------------------------------------------
# Drawing
# ----------------------------------------------------------------------

def axes_in_inches(fig, left, bottom, width, height):
    """New axes at a position given in inches from the lower left corner of the page."""
    return fig.add_axes([left / PAGE_WIDTH, bottom / PAGE_HEIGHT, width / PAGE_WIDTH, height / PAGE_HEIGHT])


def label_panel(ax, letter, title=None):
    """Panel letter (a), (b), ... at the upper left corner, and a short title above the panel."""
    ax.text(-0.02, 1.035, f'({letter})', transform=ax.transAxes, fontsize=FS_PANEL, va='bottom', ha='right')
    if title:
        ax.set_title(title, fontsize=FS)


def show_box(ax):
    """Draw all four sides of the axes, for the panels that show the simulation box."""
    for side in ['top', 'right']:
        ax.spines[side].set_visible(True)


def draw_side_view(ax, horizontal, box_width, z_low, z_high, slab_bottom, slab_top, restrained_half_width,
                   chains):
    """
    Side view to scale. horizontal is 0 (x) or 1 (y). The silica slab is an orange band with the
    restrained centermost layer hatched; every chain is a line through its backbone (SI1, C1 ... C_L)
    and its silicon (SI1) a dark blue dot.
    """
    ax.add_patch(Rectangle((0, slab_bottom), box_width, slab_top - slab_bottom, facecolor=SILICA_FILL,
                           edgecolor=SILICA_LINE, lw=0.8, zorder=1))
    ax.add_patch(Rectangle((0, -restrained_half_width), box_width, 2 * restrained_half_width, facecolor='none',
                           edgecolor=SILICA_LINE, lw=0.0, hatch='////', alpha=0.5, zorder=1.5))
    for backbone in chains:
        ax.plot(backbone[:, horizontal], backbone[:, 2], color=CHAIN_COLOR, lw=0.6, solid_capstyle='round', zorder=3)
    if chains:
        anchors = np.array([backbone[0] for backbone in chains])
        ax.scatter(anchors[:, horizontal], anchors[:, 2], s=3.5, color=ANCHOR_COLOR, lw=0, zorder=4)
    ax.axhline(0, color=TEXT_COLOR, lw=0.5, ls=':', zorder=5)
    ax.set_xlim(0, box_width)
    ax.set_ylim(z_low, z_high)
    ax.set_aspect('equal')
    show_box(ax)


def draw_face_map(ax, box_x, box_y, hydroxyl_xy, anchor_xy, filled):
    """
    Top view of one face: free hydroxyl oxygens (small gray dots) and the binding locations of the
    chains, as filled (top face) or open (bottom face) blue circles.
    """
    if len(hydroxyl_xy):
        ax.scatter(hydroxyl_xy[:, 0], hydroxyl_xy[:, 1], s=2.5, color=HYDROXYL_COLOR, lw=0, zorder=1)
    if len(anchor_xy):
        if filled:
            ax.scatter(anchor_xy[:, 0], anchor_xy[:, 1], s=15, color=CHAIN_COLOR, edgecolor='white', lw=0.5, zorder=2)
        else:
            ax.scatter(anchor_xy[:, 0], anchor_xy[:, 1], s=13, facecolor='white', edgecolor=CHAIN_COLOR, lw=0.9,
                       zorder=2)
    ax.set_xlim(0, box_x)
    ax.set_ylim(0, box_y)
    ax.set_aspect('equal')
    show_box(ax)
    ax.set_xlabel('$x$ (nm)')
    ax.set_ylabel('$y$ (nm)')


def draw_layer_diagram(ax, z_values, atom_names, z_from, z_to):
    """
    Atomic layers of the silica along z: one column per kind of atom (Si, bridging O, hydroxyl O, H),
    with a bar at every height at which atoms of that kind are found. The width of a bar grows with
    the number of atoms at that height.
    """
    edges = np.arange(z_from, z_to + 0.02, 0.02)
    counts_per_kind = []
    for names, label, color in SILICA_LAYERS:
        chosen = np.isin(atom_names, names)
        counts, edges = np.histogram(z_values[chosen], bins=edges)
        counts_per_kind.append(counts)
    largest = max(counts.max() for counts in counts_per_kind)
    centers = (edges[:-1] + edges[1:]) / 2
    for column, (counts, (names, label, color)) in enumerate(zip(counts_per_kind, SILICA_LAYERS)):
        filled = counts > 0
        half_width = 0.42 * counts[filled] / largest
        ax.barh(centers[filled], 2 * half_width, left=column - half_width, height=0.018, color=color, lw=0)
    ax.axhline(0, color=TEXT_COLOR, lw=0.5, ls=':')
    ax.set_xticks(range(len(SILICA_LAYERS)))
    ax.set_xticklabels([label for names, label, color in SILICA_LAYERS], fontsize=FS_SMALL)
    ax.tick_params(axis='x', which='both', length=0)
    ax.set_xlim(-0.55, len(SILICA_LAYERS) - 0.45)
    ax.set_ylim(z_from, z_to)
    ax.set_ylabel('$z$ (nm)')


def draw_table(ax, title, header, rows, column_x, line_height, title_note=None):
    """
    A table in the style of the tables of a journal article: a title, a top rule, the header,
    a thin rule, the rows and a bottom rule. column_x are the left edges of the columns (0-1).
    title_note is a short text at the right end of the title line. Returns the height (axes
    units) below the bottom rule.
    """
    ax.axis('off')
    y = 1.0
    ax.text(0, y, title, fontsize=FS + 1, va='top', transform=ax.transAxes)
    if title_note:
        ax.text(1, y - 0.15 * line_height, title_note, fontsize=FS_SMALL, color=TEXT_SOFT, va='top', ha='right',
                transform=ax.transAxes)
    y -= 1.5 * line_height
    ax.plot([0, 1], [y, y], color=TEXT_COLOR, lw=0.9, transform=ax.transAxes, clip_on=False)
    y -= 0.25 * line_height
    for x, text in zip(column_x, header):
        ax.text(x, y, text, fontsize=FS_SMALL, color=TEXT_SOFT, va='top', transform=ax.transAxes)
    y -= line_height
    ax.plot([0, 1], [y + 0.1 * line_height] * 2, color=RULE_COLOR, lw=0.5, transform=ax.transAxes, clip_on=False)
    y -= 0.15 * line_height
    for row in rows:
        for x, text in zip(column_x, row):
            ax.text(x, y, text, fontsize=FS_SMALL, va='top', transform=ax.transAxes)
        y -= line_height
    ax.plot([0, 1], [y + 0.1 * line_height] * 2, color=TEXT_COLOR, lw=0.9, transform=ax.transAxes, clip_on=False)
    return y


def workflow_rows(workflow):
    """Rows of the workflow table: 'sam' (dry and wet stages) or 'surface' (relaxation of a new slab)."""
    if workflow == 'surface':
        stages = [('slab', stage) for stage in SURFACE_RELAXATION]
    else:
        stages = []
        for stage in DRY_RELAXATION + DRY_PRODUCTION:
            stages.append(('dry', stage))
        stages.append(('wet', None))
        for stage in WET_RELAXATION + WET_PRODUCTION:
            stages.append(('wet', stage))

    rows = []
    for state, stage in stages:
        if stage is None:
            rows.append(['wet', '1_solvation', 'gmx solvate, TIP3P', '-', '-', '-'])
            continue
        if stage['dt_ps'] is None:
            dt, length = '-', '-'
        else:
            dt, length = f'{stage["dt_ps"] * 1000:g}', f'{stage["length_ns"]:g}'
        if stage['define'] == '-DPOSRES_SURFACE':
            restraint = 'whole slab'
        elif stage['define'] == '-DPOSRES':
            restraint = 'middle 1 nm'
        else:
            restraint = 'none'
        run = stage_description(stage).replace('energy minimization (steepest descent)', 'minimization')
        run = run.replace(', semi-isotropic', '').replace('Parrinello-Rahman', 'P-R')
        rows.append([state, stage['folder'], run, dt, length, restraint])
    return rows


# ----------------------------------------------------------------------
# The sheet
# ----------------------------------------------------------------------

def plot_system(gro_file, output_file, chain_length, coverage, sam_resname='OTS', surface_resname='SUR',
                restrained_layer_nm=1.0, settings=None, workflow='sam'):
    """
    Save the overview sheet of a built system as output_file (.pdf or .png), one A4 page.

    gro_file            : the built structure (system/func.gro)
    chain_length        : number of carbon atoms in the chains, L
    coverage            : target coverage (molecules/nm^2)
    restrained_layer_nm : width of the restrained centermost layer of the slab (nm)
    settings            : build settings shown in the header (surface_mode, random_seed, ...), optional
    workflow            : 'sam' (dry and wet stages), 'surface' (relaxation of a new slab) or None
    """
    if settings is None:
        settings = {}
    atoms = read_gro(gro_file)
    box_x, box_y, box_z = atoms['box']
    area = box_x * box_y
    residue_names = np.array(atoms['residue_name'])
    atom_names = np.array(atoms['atom_name'])
    xyz = atoms['xyz'].copy()

    # Move z = 0 to the center of the silica slab (midpoint of its lowest and highest atom).
    is_silica = residue_names == surface_resname
    slab_center = (xyz[is_silica, 2].max() + xyz[is_silica, 2].min()) / 2
    xyz[:, 2] -= slab_center
    atoms['xyz'] = xyz
    z_low, z_high = -slab_center, box_z - slab_center
    slab_bottom, slab_top = xyz[is_silica, 2].min(), xyz[is_silica, 2].max()

    # chains of each face, their binding locations (SI1) and the free hydroxyl oxygens (OD) of each face
    backbones = chain_backbones(atoms, sam_resname, chain_length)
    top_chains = [b for b in backbones if b[0, 2] > 0]
    bottom_chains = [b for b in backbones if b[0, 2] <= 0]
    n_chains = len(backbones)
    is_hydroxyl = is_silica & (atom_names == 'OD')
    hydroxyl_top = xyz[is_hydroxyl & (xyz[:, 2] > 0)]
    hydroxyl_bottom = xyz[is_hydroxyl & (xyz[:, 2] <= 0)]
    anchors_top = np.array([b[0] for b in top_chains]).reshape(-1, 3)
    anchors_bottom = np.array([b[0] for b in bottom_chains]).reshape(-1, 3)

    # SAM thickness on each face: highest (lowest) chain carbon minus the top (bottom) silica atom
    sam_top = max([b[:, 2].max() for b in top_chains], default=slab_top)
    sam_bottom = min([b[:, 2].min() for b in bottom_chains], default=slab_bottom)

    # restrained centermost layer: silica atoms within half the layer width of the slab center,
    # the same rule as in builder.find_middle_layer
    half_width = restrained_layer_nm / 2
    n_restrained = int(np.sum(is_silica & (xyz[:, 2] > -half_width) & (xyz[:, 2] < half_width)))

    is_sam = residue_names == sam_resname
    masses = np.array([MASSES[element_of(name)] for name in atom_names])

    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(PAGE_WIDTH, PAGE_HEIGHT))
        left_edge = 0.6 / PAGE_WIDTH
        right_edge = 1 - 0.45 / PAGE_WIDTH

        # ------------------------------------------------------------ header
        if n_chains > 0:
            title = f'SAM system,  $L$ = {chain_length},  $\\Gamma$ = {coverage:.2f} molecules/nm²'
            subtitle = ('Alkylsilane SAMs on both exposed (001) faces of a hydroxylated α-quartz slab. Both faces have '
                        'the target\nsurface coverage; the binding locations of the two faces were chosen independently.')
        elif workflow == 'surface':
            title = 'New silica slab, before relaxation'
            subtitle = ('Hydroxylated α-quartz (001) slab built from unit cells. The slab is relaxed in GROMACS before '
                        'alkylsilane\nchains are grafted on it (surface_mode=\'user\').')
        else:
            title = 'Pristine silica slab,  $\\Gamma$ = 0'
            subtitle = 'Hydroxylated α-quartz (001) slab without chains, prepared for the dry and the wet simulations.'
        fig.text(left_edge, 11.22 / PAGE_HEIGHT, title, fontsize=FS_TITLE, va='baseline')
        fig.text(left_edge, 11.04 / PAGE_HEIGHT, subtitle, fontsize=FS, color=TEXT_SOFT, va='top', linespacing=1.4)
        details = []
        for key in ['surface_mode', 'x_tiles', 'y_tiles', 'random_seed', 'system_name']:
            if key in settings:
                details.append(f'{key} = {settings[key]}')
        if details:
            fig.text(left_edge, 10.6 / PAGE_HEIGHT, 'build settings:   ' + '   '.join(details), fontsize=FS_SMALL,
                     color=TEXT_SOFT, va='baseline')
        fig.add_artist(plt.Line2D([left_edge, right_edge], [10.48 / PAGE_HEIGHT] * 2, color=TEXT_COLOR, lw=0.9))

        # ------------------------------------------------------------ row 1: side views and density
        # both side views use the same scale (inches per nm), as large as the space allows
        row_top = 10.05
        scale = min(3.55 / box_z, 4.2 / (box_x + box_y))
        view_height = box_z * scale
        view_bottom = row_top - view_height
        ax_xz = axes_in_inches(fig, 0.75, view_bottom, box_x * scale, view_height)
        ax_yz = axes_in_inches(fig, 0.75 + box_x * scale + 0.15, view_bottom, box_y * scale, view_height)
        labels_left = 0.75 + (box_x + box_y) * scale + 0.15 + 0.06
        density_left = labels_left + 0.95
        ax_rho = axes_in_inches(fig, density_left, view_bottom, 7.82 - density_left, view_height)

        draw_side_view(ax_xz, 0, box_x, z_low, z_high, slab_bottom, slab_top, half_width, backbones)
        draw_side_view(ax_yz, 1, box_y, z_low, z_high, slab_bottom, slab_top, half_width, backbones)
        ax_xz.set_xlabel('$x$ (nm)')
        ax_xz.set_ylabel('$z$ (nm)')
        ax_yz.set_xlabel('$y$ (nm)')
        ax_yz.tick_params(labelleft=False)
        label_panel(ax_xz, 'a', 'side view, $xz$-plane')
        label_panel(ax_yz, 'b', 'side view, $yz$-plane')
        ax_xz.text(0.97, 0.015, f'{box_x:.2f} × {box_y:.2f} × {box_z:.2f} nm³', transform=ax_xz.transAxes,
                   ha='right', va='bottom', fontsize=FS_SMALL - 0.5, color=TEXT_SOFT)

        # labels of the layers, to the right of the side views
        def layer_label(z, text, color, size=FS_SMALL):
            fig.text(labels_left / PAGE_WIDTH, (view_bottom + (z - z_low) * scale) / PAGE_HEIGHT, text,
                     fontsize=size, color=color, va='center', ha='left')

        layer_label(0.0, 'silica\nslab', SILICA_TEXT, FS)
        layer_label(-half_width - 0.35, f'restrained\n{restrained_layer_nm:g} nm', SILICA_LINE)
        if top_chains:
            layer_label((slab_top + sam_top) / 2 + 0.25, 'SAM', CHAIN_COLOR, FS)
        if bottom_chains:
            layer_label((slab_bottom + sam_bottom) / 2 - 0.25, 'SAM', CHAIN_COLOR, FS)
        layer_label((sam_top + z_high) / 2, 'vacuum', TEXT_SOFT)
        layer_label((sam_bottom + z_low) / 2, 'vacuum', TEXT_SOFT)

        # (c) mass density along z: silica as a filled band, the SAM as a line
        z_edges = np.arange(z_low, z_high + 0.1, 0.1)
        z_centers = (z_edges[:-1] + z_edges[1:]) / 2
        rho_silica = mass_density_profile(xyz[is_silica, 2], masses[is_silica], area, z_edges)
        ax_rho.fill_betweenx(z_centers, 0, rho_silica, color=SILICA_FILL, lw=0, zorder=1)
        ax_rho.plot(rho_silica, z_centers, color=SILICA_LINE, lw=0.8, zorder=2)
        key = [Patch(facecolor=SILICA_FILL, edgecolor=SILICA_LINE, lw=0.8, label='silica')]
        if n_chains > 0:
            rho_sam = mass_density_profile(xyz[is_sam, 2], masses[is_sam], area, z_edges)
            ax_rho.plot(rho_sam, z_centers, color=CHAIN_COLOR, lw=1.3, zorder=3)
            key.append(plt.Line2D([], [], color=CHAIN_COLOR, lw=1.3, label='SAM'))
        ax_rho.legend(handles=key, loc='upper right')
        ax_rho.axhline(0, color=TEXT_COLOR, lw=0.5, ls=':')
        ax_rho.set_ylim(z_low, z_high)
        ax_rho.set_xlim(0, rho_silica.max() * 1.12)
        ax_rho.set_xlabel('$\\rho$ (kg/m³)')
        ax_rho.set_ylabel('$z$ (nm)')
        label_panel(ax_rho, 'c', 'mass density')

        # ------------------------------------------------------------ row 2: the two faces
        row_top = view_bottom - 0.75
        map_scale = min(1.85 / box_y, 2.2 / box_x)
        map_height = box_y * map_scale
        map_bottom = row_top - map_height
        ax_top = axes_in_inches(fig, 0.75, map_bottom, box_x * map_scale, map_height)
        ax_bottom = axes_in_inches(fig, 0.75 + box_x * map_scale + 0.55, map_bottom, box_x * map_scale, map_height)
        draw_face_map(ax_top, box_x, box_y, hydroxyl_top, anchors_top, filled=True)
        draw_face_map(ax_bottom, box_x, box_y, hydroxyl_bottom, anchors_bottom, filled=False)
        if n_chains > 0:
            label_panel(ax_top, 'd', f'top face, {len(anchors_top)} chains')
            label_panel(ax_bottom, 'e', f'bottom face, {len(anchors_bottom)} chains')
        else:
            label_panel(ax_top, 'd', f'top face, {len(hydroxyl_top)} OH groups')
            label_panel(ax_bottom, 'e', f'bottom face, {len(hydroxyl_bottom)} OH groups')

        # key of the two face maps, below them
        key = [plt.Line2D([], [], ls='none', marker='o', ms=4.4, color=CHAIN_COLOR, mec='white', mew=0.5,
                          label='binding location, top face'),
               plt.Line2D([], [], ls='none', marker='o', ms=4.0, mfc='white', mec=CHAIN_COLOR, mew=0.9,
                          label='binding location, bottom face'),
               plt.Line2D([], [], ls='none', marker='o', ms=2.2, color=HYDROXYL_COLOR, mew=0,
                          label='free hydroxyl oxygen')]
        if n_chains == 0:
            key = key[2:]
        fig.legend(handles=key, loc='upper left', ncol=3, fontsize=FS_SMALL, handletextpad=0.2, columnspacing=1.4,
                   bbox_to_anchor=(0.62 / PAGE_WIDTH, (map_bottom - 0.42) / PAGE_HEIGHT))

        # (f) nearest-neighbor distance of the binding locations on each face, or the atomic layers
        # of a slab without chains
        ax_nn = axes_in_inches(fig, density_left, map_bottom, 7.82 - density_left, map_height)
        if n_chains > 0:
            all_distances = []
            for anchors in [anchors_top, anchors_bottom]:
                if len(anchors) > 1:
                    all_distances.append(nearest_neighbor_distances(anchors[:, :2], box_x, box_y))
            if all_distances:
                largest = max(d.max() for d in all_distances)
                bins = np.linspace(0, largest * 1.05, 26)
                for distances, line_style, label in zip(all_distances, ['-', '--'], ['top face', 'bottom face']):
                    ax_nn.hist(distances, bins=bins, histtype='step', color=CHAIN_COLOR, ls=line_style, lw=1.2,
                               label=label)
                # room above the bars for the legend
                ax_nn.set_ylim(0, ax_nn.get_ylim()[1] * 1.45)
                # mean nearest-neighbor distance of points placed at random with the same density
                random_mean = 1 / (2 * np.sqrt(len(anchors_top) / area))
                ax_nn.axvline(random_mean, color=TEXT_SOFT, lw=0.8, ls=':')
                ax_nn.text(random_mean, ax_nn.get_ylim()[1] * 0.97, 'random\narrangement ', fontsize=FS_SMALL - 0.5,
                           color=TEXT_SOFT, va='top', ha='right', linespacing=1.2)
                ax_nn.legend(loc='upper right')
            ax_nn.set_xlabel('nearest-neighbor distance (nm)')
            ax_nn.set_ylabel('chains')
            label_panel(ax_nn, 'f', 'spacing of the chains')
        else:
            draw_layer_diagram(ax_nn, xyz[is_silica, 2], atom_names[is_silica], slab_bottom - 0.15, slab_top + 0.15)
            label_panel(ax_nn, 'f', 'atomic layers of the slab')

        # ------------------------------------------------------------ row 3: tables
        table_top = map_bottom - 0.78
        table_height = table_top - 0.6
        ax_numbers = axes_in_inches(fig, 0.6, 0.6, 3.35, table_height)
        ax_stages = axes_in_inches(fig, 4.3, 0.6, 3.52, table_height)
        line_height = 0.148 / table_height

        rows = []
        if n_chains > 0:
            rows += [
                ['Chain length $L$', f'{chain_length} carbon atoms'],
                ['Target coverage $\\Gamma$', f'{coverage:.2f} molecules/nm²'],
                ['Chains, top / bottom face', f'{len(anchors_top)} / {len(anchors_bottom)}'],
                ['Coverage, top / bottom face',
                 f'{len(anchors_top) / area:.3f} / {len(anchors_bottom) / area:.3f} molecules/nm²'],
            ]
        rows += [
            ['Free OH groups, top / bottom face',
             f'{len(hydroxyl_top) / area:.2f} / {len(hydroxyl_bottom) / area:.2f} per nm²'],
            ['Simulation box', f'{box_x:.4f} × {box_y:.4f} × {box_z:.4f} nm³'],
            ['Face area', f'{area:.3f} nm²'],
            ['Silica slab thickness', f'{slab_top - slab_bottom:.2f} nm'],
        ]
        if n_chains > 0:
            rows.append(['SAM thickness, top / bottom face',
                         f'{sam_top - slab_top:.2f} / {slab_bottom - sam_bottom:.2f} nm (all-trans)'])
        rows += [
            ['Atoms: silica / chains / total', f'{is_silica.sum()} / {is_sam.sum()} / {len(atom_names)}'],
            ['Restrained centermost layer', f'{n_restrained} atoms ({restrained_layer_nm:g} nm)'],
            ['Slab center in the GROMACS box', f'$z$ = {slab_center:.4f} nm'],
        ]
        draw_table(ax_numbers, 'System', ['Quantity', 'Value'], rows, [0.0, 0.5], line_height)

        if workflow is not None:
            stage_rows = workflow_rows(workflow)
            title = 'Relaxation of the slab' if workflow == 'surface' else 'Simulation workflow'
            header = ['State', 'Stage folder', 'Run', 'dt (fs)', 'ns', 'Restrained']
            column_x = [0.0, 0.09, 0.37, 0.64, 0.74, 0.82]
            if workflow == 'surface':
                y = draw_table(ax_stages, title, header, stage_rows, column_x, line_height,
                               title_note='run: sbatch run_surface_relaxation.sh')
                note = ('Relaxed slab: surface/2_npt_crescale/confout.gro, used with\n'
                        "build_sam_system(..., surface_mode='user', user_surface=...)")
                ax_stages.text(0, y - 0.5 * line_height, note, fontsize=FS_SMALL, color=TEXT_SOFT, va='top',
                               transform=ax_stages.transAxes, linespacing=1.5)
            else:
                draw_table(ax_stages, title, header, stage_rows, column_x, line_height,
                           title_note='run: bash submit_all.sh')
        else:
            ax_stages.axis('off')

        # ------------------------------------------------------------ footer
        fig.add_artist(plt.Line2D([left_edge, right_edge], [0.48 / PAGE_HEIGHT] * 2, color=RULE_COLOR, lw=0.6))
        fig.text(left_edge, 0.32 / PAGE_HEIGHT, 'apadana 4.0', fontsize=FS_SMALL, color=TEXT_SOFT, va='baseline')
        fig.text(right_edge, 0.32 / PAGE_HEIGHT,
                 f'Lengths in nm; $z$ = 0 at the center of the silica slab ($z$ in the GROMACS files = $z$ + '
                 f'{slab_center:.4f} nm).', fontsize=FS_SMALL, color=TEXT_SOFT, va='baseline', ha='right')

        fig.savefig(output_file)
        plt.close(fig)
