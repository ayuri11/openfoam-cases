import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import get_cmap

power_3d = np.load('power_density_3d.npy')
midplane = power_3d[:, :, 14] / 1e6

# reconstruct full core from 1/12 wedge
full_core = np.zeros_like(midplane)
for i in range(12):
    rotated = rotate(midplane, angle=i*30, reshape=False, order=1)
    full_core += rotated

# =============================================================================
# CUSTOM COLORMAP — inferno but starting from warm white instead of black
# zero power = soft warm white, high power = orange/yellow
# =============================================================================
inferno = get_cmap('inferno')
colors_inferno = inferno(np.linspace(0, 1, 256))
# blend the bottom 40% toward a warm white
for i in range(100):
    t = i / 100
    colors_inferno[i] = (
        1.0 - t * (1.0 - colors_inferno[i][0]),   # R
        1.0 - t * (1.0 - colors_inferno[i][1]),   # G
        0.95 - t * (0.95 - colors_inferno[i][2]), # B — slightly warm white
        1.0
    )
custom_cmap = LinearSegmentedColormap.from_list('inferno_light', colors_inferno)

# =============================================================================
# THEME
# =============================================================================
bg_outer    = '#3a3f4a'   # dark grey outer — matches axial plot
bg_panel    = '#f0f2f5'   # very light grey-white inner panel
text_color  = '#ffffff'
spine_color = '#6a7a8a'
title_color = '#ffffff'

fig, ax = plt.subplots(figsize=(7, 6))
fig.patch.set_facecolor(bg_outer)
ax.set_facecolor(bg_panel)

# =============================================================================
# MAIN PLOT — custom colormap, light bg
# =============================================================================
im = ax.imshow(
    full_core.T,
    origin='lower',
    extent=[-60, 60, -60, 60],
    cmap=custom_cmap,
    interpolation='bilinear',
    vmin=0
)

# =============================================================================
# COLORBAR
# =============================================================================
cbar = plt.colorbar(im, ax=ax, pad=0.02)
cbar.set_label('Power Density (MW/m³)',
               fontsize=13, color=text_color, labelpad=10)
cbar.ax.yaxis.set_tick_params(color=text_color, labelcolor=text_color)
cbar.outline.set_edgecolor(spine_color)
plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'),
         color=text_color, fontsize=11)

# =============================================================================
# BOUNDARY CIRCLES — dark navy, thick, clearly visible on light bg
# =============================================================================
theta = np.linspace(0, 2*np.pi, 300)
ax.plot(45*np.cos(theta), 45*np.sin(theta),
        '--', color='#1a2535', linewidth=2.2,
        label='Core boundary (r = 45 cm)')
ax.plot(60*np.cos(theta), 60*np.sin(theta),
        ':', color='#1a2535', linewidth=2.0,
        label='Reflector boundary (r = 60 cm)')

# =============================================================================
# LABELS
# =============================================================================
ax.set_xlabel('x (cm)', fontsize=14, color=text_color, labelpad=10)
ax.set_ylabel('y (cm)', fontsize=14, color=text_color, labelpad=10)
ax.set_title(
    'TRAVIS Core — Reconstructed Full-Core Radial Power\n(from 1/12 symmetry, at z = 0)',
    fontsize=13, color=title_color, pad=13
)

# =============================================================================
# TICKS AND SPINES
# =============================================================================
ax.tick_params(axis='both', colors=text_color,
               labelsize=12, width=1.8, length=5)
ax.grid(False)
for spine in ax.spines.values():
    spine.set_edgecolor(spine_color)
    spine.set_linewidth(2.0)

# =============================================================================
# LEGEND — light box with dark text, sits on light panel
# =============================================================================
legend = ax.legend(
    fontsize=10,
    loc='lower right',
    facecolor='#e8edf2',
    edgecolor='#1a2535',
    framealpha=0.92
)
for text in legend.get_texts():
    text.set_color('#1a2535')

plt.tight_layout()
plt.savefig('radial_power_white_theme.png', dpi=300, facecolor=bg_outer)
plt.savefig('radial_power_white_theme.pdf', facecolor=bg_outer)
print("Saved: radial_power_white_theme.png and .pdf")
plt.show()
