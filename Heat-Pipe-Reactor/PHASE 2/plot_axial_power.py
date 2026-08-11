import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load the CSV saved by Phase 2
df = pd.read_csv('axial_power_profile.csv')

# Expected columns: z_cm, power_density_W_per_m3
z     = df.iloc[:, 0].values   # z position in cm
power = df.iloc[:, 1].values   # power density in W/m³

# Convert to more readable units
power_MW = power / 1e6         # W/m³ → MW/m³

# =============================================================================
# DARK THEME SETUP
# =============================================================================
bg_color    = '#2e2e2e'
text_color  = '#cccccc'
spine_color = '#666666'
title_color = '#ffffff'
line_color  = '#5b9bd5'   # blue matching poster theme

fig, ax = plt.subplots(figsize=(5, 8))
fig.patch.set_facecolor(bg_color)
ax.set_facecolor(bg_color)

# =============================================================================
# MAIN PLOT — same structure as original, just recolored
# =============================================================================
ax.plot(power_MW, z, 'o-', color=line_color, linewidth=2, markersize=5,
        markerfacecolor=line_color, markeredgecolor=bg_color, markeredgewidth=1.5)

# Shade under the curve for visual weight
ax.fill_betweenx(z, power_MW, alpha=0.15, color=line_color)

# =============================================================================
# REFERENCE LINES
# =============================================================================
ax.axhline(y=0,   color='#888888', linestyle='--',
           linewidth=0.8, label='Core midplane')
ax.axhline(y=80,  color='#666666', linestyle=':',
           linewidth=0.8, label='Fuel boundary')
ax.axhline(y=-80, color='#666666', linestyle=':',
           linewidth=0.8)

# =============================================================================
# LABELS, TICKS, SPINES
# =============================================================================
ax.set_xlabel('Power Density (MW/m³)', fontsize=12,
              color=text_color, labelpad=8)
ax.set_ylabel('Axial Position z (cm)', fontsize=12,
              color=text_color, labelpad=8)
ax.set_title('AYURI Core — Axial Power Distribution\n(normalized to 2.5 MWth)',
             fontsize=12, color=title_color, pad=12)

ax.tick_params(axis='both', colors=text_color, labelsize=10)
ax.grid(False)

for spine in ax.spines.values():
    spine.set_edgecolor(spine_color)
    spine.set_linewidth(1.5)

# =============================================================================
# LEGEND
# =============================================================================
ax.legend(fontsize=10,
          facecolor='#3a3a3a',
          edgecolor='#555555',
          labelcolor=text_color,
          framealpha=0.85)

plt.tight_layout()
plt.savefig('axial_power_profile_dark.png', dpi=300, facecolor=bg_color)
plt.savefig('axial_power_profile_dark.pdf', facecolor=bg_color)
print("Saved: axial_power_profile_dark.png and axial_power_profile_dark.pdf")
plt.show()
