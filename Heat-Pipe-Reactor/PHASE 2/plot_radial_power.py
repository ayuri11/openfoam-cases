import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate

power_3d = np.load('power_density_3d.npy')
midplane = power_3d[:, :, 14] / 1e6

# reconstruct full core from 1/12 wedge
full_core = np.zeros_like(midplane)
for i in range(12):
    rotated = rotate(midplane, angle=i*30, reshape=False, order=1)
    full_core += rotated

# =============================================================================
# THEME — dark grey outer matching ParaView, light inner panel for contrast
# =============================================================================
bg_outer   = '#3a3f4a'   # dark grey outer — same as axial plot
bg_panel   = '#d0d8e0'   # light blue-grey inner panel — so inferno colors pop
                          # and boundary circles are clearly visible
text_color = '#ffffff'    # white for outer labels
tick_color = '#ffffff'    # white ticks
spine_color= '#6a7a8a'    # medium grey spines
title_color= '#ffffff'

fig, ax = plt.subplots(figsize=(7, 6))
fig.patch.set_facecolor(bg_outer)
ax.set_facecolor(bg_panel)

# =============================================================================
# MAIN PLOT — inferno colormap on light background
# =============================================================================
im = ax.imshow(
    full_core.T,
    origin='lower',
    extent=[-60, 60, -60, 60],
    cmap='inferno',
    interpolation='bilinear'
)

# =============================================================================
# COLORBAR — styled to match outer theme
# =============================================================================
cbar = plt.colorbar(im, ax=ax, pad=0.02)
cbar.set_label('Power Density (MW/m³)',
               fontsize=13, color=text_color, labelpad=10)
cbar.ax.yaxis.set_tick_params(color=text_color, labelcolor=text_color)
cbar.outline.set_edgecolor(spine_color)
plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'),
         color=text_color, fontsize=11)

# =============================================================================
# BOUNDARY CIRCLES — dark colors now visible on light panel
# =============================================================================
theta = np.linspace(0, 2*np.pi, 300)
ax.plot(45*np.cos(theta), 45*np.sin(theta),
        '--', color='#1a2a3a', linewidth=2.0,
        label='Core boundary (r = 45 cm)')
ax.plot(60*np.cos(theta), 60*np.sin(theta),
        ':', color='#3a4a5a', linewidth=1.8,
        label='Reflector boundary (r = 60 cm)')

# =============================================================================
# LABELS — large white text on dark outer
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
ax.tick_params(axis='both', colors=tick_color,
               labelsize=12, width=1.8, length=5)
ax.grid(False)
for spine in ax.spines.values():
    spine.set_edgecolor(spine_color)
    spine.set_linewidth(2.0)

# =============================================================================
# LEGEND — dark box on light panel — use dark text inside
# =============================================================================
legend = ax.legend(
    fontsize=10,
    loc='lower right',
    facecolor='#e8edf2',
    edgecolor='#3a4a5a',
    framealpha=0.92
)
for text in legend.get_texts():
    text.set_color('#1a2a3a')   # dark text inside legend on light bg

plt.tight_layout()
plt.savefig('radial_power_paraview_theme.png', dpi=300, facecolor=bg_outer)
plt.savefig('radial_power_paraview_theme.pdf', facecolor=bg_outer)
print("Saved: radial_power_paraview_theme.png and .pdf")
plt.show()
