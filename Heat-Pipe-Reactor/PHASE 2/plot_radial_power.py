import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate

power_3d = np.load('power_density_3d.npy')
midplane = power_3d[:, :, 14] / 1e6

# The tally captured only the 30° wedge (1/12 of core)
# Reconstruct full core by rotating and adding 11 copies
full_core = np.zeros_like(midplane)
for i in range(12):
    rotated = rotate(midplane, angle=i*30, reshape=False, order=1)
    full_core += rotated

# =============================================================================
# DARK THEME SETUP — matches axial profile aesthetic
# =============================================================================
bg_color    = '#2e2e2e'
text_color  = '#cccccc'
spine_color = '#666666'
title_color = '#ffffff'

fig, ax = plt.subplots(figsize=(7, 6))
fig.patch.set_facecolor(bg_color)
ax.set_facecolor(bg_color)

# =============================================================================
# MAIN PLOT
# =============================================================================
im = ax.imshow(
    full_core.T,
    origin='lower',
    extent=[-60, 60, -60, 60],
    cmap='inferno',
    interpolation='bilinear'
)

# =============================================================================
# COLORBAR — dark styled
# =============================================================================
cbar = plt.colorbar(im, ax=ax, pad=0.02)
cbar.set_label('Power Density (MW/m³)', fontsize=12, color=text_color, labelpad=10)
cbar.ax.yaxis.set_tick_params(color=text_color, labelcolor=text_color)
cbar.outline.set_edgecolor(spine_color)
plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color=text_color, fontsize=10)

# =============================================================================
# BOUNDARY CIRCLES
# =============================================================================
theta = np.linspace(0, 2*np.pi, 300)
ax.plot(45*np.cos(theta), 45*np.sin(theta),
        '--', color='#aaaaaa', linewidth=1.2, label='Core boundary (r = 45 cm)')
ax.plot(60*np.cos(theta), 60*np.sin(theta),
        ':', color='#777777', linewidth=1.0, label='Reflector boundary (r = 60 cm)')

# =============================================================================
# LABELS, TICKS, SPINES
# =============================================================================
ax.set_xlabel('x (cm)', fontsize=12, color=text_color, labelpad=8)
ax.set_ylabel('y (cm)', fontsize=12, color=text_color, labelpad=8)
ax.set_title(
    'AYURI Core — Reconstructed Full-Core Radial Power\n(from 1/12 symmetry, at z = 0)',
    fontsize=12, color=title_color, pad=12
)
ax.tick_params(axis='both', colors=text_color, labelsize=10, which='both')
ax.grid(False)
for spine in ax.spines.values():
    spine.set_edgecolor(spine_color)
    spine.set_linewidth(1.5)

# =============================================================================
# LEGEND
# =============================================================================
ax.legend(
    fontsize=9,
    loc='lower right',
    facecolor='#3a3a3a',
    edgecolor='#555555',
    labelcolor=text_color,
    framealpha=0.85
)

plt.tight_layout()
plt.savefig('radial_power_map_dark.png', dpi=300, facecolor=bg_color)
plt.savefig('radial_power_map_dark.pdf', facecolor=bg_color)
print("Saved: radial_power_map_dark.png and radial_power_map_dark.pdf")
plt.show()
