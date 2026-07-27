import numpy as np
import matplotlib.pyplot as plt

power_3d = np.load('power_density_3d.npy')   # shape (50, 50, 28)

# Take midplane slice (z index 14 = center)
midplane = power_3d[:, :, 14] / 1e6          # W/m³ → MW/m³

# =============================================================================
# DARK THEME SETUP — matches axial profile aesthetic
# =============================================================================
bg_color    = '#2e2e2e'   # outer figure and axes background
text_color  = '#cccccc'   # axis labels and tick labels
spine_color = '#666666'   # axis spine lines
title_color = '#ffffff'

fig, ax = plt.subplots(figsize=(7, 6))
fig.patch.set_facecolor(bg_color)
ax.set_facecolor(bg_color)

# =============================================================================
# MAIN PLOT — same inferno colormap, same data, restyled
# =============================================================================
im = ax.imshow(
    midplane.T,
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
# BOUNDARY CIRCLES — same as before, lighter to sit on dark bg
# =============================================================================
theta = np.linspace(0, 2*np.pi, 300)
ax.plot(45*np.cos(theta), 45*np.sin(theta),
        '--', color='#aaaaaa', linewidth=1.2, label='Core boundary (r = 45 cm)')
ax.plot(60*np.cos(theta), 60*np.sin(theta),
        ':', color='#777777', linewidth=1.0, label='Reflector boundary (r = 60 cm)')

# =============================================================================
# AXES LABELS AND TITLE
# =============================================================================
ax.set_xlabel('x (cm)', fontsize=12, color=text_color, labelpad=8)
ax.set_ylabel('y (cm)', fontsize=12, color=text_color, labelpad=8)
ax.set_title(
    'AYURI Core — Reconstructed Full-Core Radial Power\n(from 1/12 symmetry, at z = 0)',
    fontsize=12, color=title_color, pad=12
)

# =============================================================================
# TICKS — light color, no inner grid
# =============================================================================
ax.tick_params(axis='both', colors=text_color, labelsize=10, which='both')
ax.grid(False)

# =============================================================================
# SPINES — styled to match axial profile (clean lines, dark color)
# =============================================================================
for spine in ax.spines.values():
    spine.set_edgecolor(spine_color)
    spine.set_linewidth(1.5)

# =============================================================================
# LEGEND — dark box matching background
# =============================================================================
legend = ax.legend(
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
