import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load the CSV saved by Phase 2
df = pd.read_csv('axial_power_profile.csv')
z     = df.iloc[:, 0].values
power = df.iloc[:, 1].values
power_MW = power / 1e6

# =============================================================================
# DARK THEME — poster-optimized for print readability
# =============================================================================
bg_outer    = '#1a2535'   # outer figure background — dark navy
bg_panel    = '#243044'   # axes panel — slightly lighter so fill stands out
text_color  = '#e8edf2'   # bright light grey — high contrast on dark bg
spine_color = '#5a7a9a'   # medium blue-grey for axis lines
ref_color   = '#8aabcc'   # lighter blue for reference lines
title_color = '#ffffff'
line_color  = '#5b9bd5'   # blue line
fill_color  = '#3a6a9a'   # distinctly darker blue for fill — visible but not competing

fig, ax = plt.subplots(figsize=(6, 9))
fig.patch.set_facecolor(bg_outer)
ax.set_facecolor(bg_panel)

# =============================================================================
# FILL — slightly darker blue so curve stands out clearly on top
# =============================================================================
ax.fill_betweenx(z, power_MW, alpha=0.55, color=fill_color)

# =============================================================================
# MAIN LINE — thicker for print visibility
# =============================================================================
ax.plot(power_MW, z, 'o-',
        color=line_color,
        linewidth=3.0,          # thicker line
        markersize=6,           # bigger dots
        markerfacecolor=line_color,
        markeredgecolor=bg_outer,
        markeredgewidth=1.5)

# =============================================================================
# REFERENCE LINES — thicker and brighter
# =============================================================================
ax.axhline(y=0,
           color=ref_color, linestyle='--',
           linewidth=1.4, label='Core midplane', alpha=0.9)
ax.axhline(y=80,
           color=ref_color, linestyle=':',
           linewidth=1.2, label='Fuel boundary', alpha=0.7)
ax.axhline(y=-80,
           color=ref_color, linestyle=':',
           linewidth=1.2, alpha=0.7)

# =============================================================================
# LABELS — larger font sizes for poster print
# =============================================================================
ax.set_xlabel('Power Density (MW/m³)',
              fontsize=15, color=text_color, labelpad=10)
ax.set_ylabel('Axial Position z (cm)',
              fontsize=15, color=text_color, labelpad=10)
ax.set_title('AYURI Core — Axial Power Distribution\n(normalized to 2.5 MWth)',
             fontsize=14, color=title_color, pad=14)

# =============================================================================
# TICKS — larger and brighter
# =============================================================================
ax.tick_params(axis='both', colors=text_color,
               labelsize=13, width=1.5, length=5)
ax.grid(False)

# =============================================================================
# SPINES — brighter and thicker
# =============================================================================
for spine in ax.spines.values():
    spine.set_edgecolor(spine_color)
    spine.set_linewidth(2.0)

# =============================================================================
# LEGEND — readable box
# =============================================================================
ax.legend(fontsize=12,
          facecolor='#2a3d55',
          edgecolor=spine_color,
          labelcolor=text_color,
          framealpha=0.92,
          loc='upper left')

plt.tight_layout()
plt.savefig('axial_power_profile_dark_v2.png', dpi=300, facecolor=bg_outer)
plt.savefig('axial_power_profile_dark_v2.pdf', facecolor=bg_outer)
print("Saved: axial_power_profile_dark_v2.png and .pdf")
plt.show()
