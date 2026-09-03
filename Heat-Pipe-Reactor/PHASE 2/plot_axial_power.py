import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load the CSV saved by Phase 2
df = pd.read_csv('axial_power_profile.csv')
z        = df.iloc[:, 0].values
power    = df.iloc[:, 1].values
power_MW = power / 1e6

# =============================================================================
# PARAVIEW-MATCHED THEME
# sampled directly from the neutronics figure:
# outer bg: #6b7280 (medium grey)
# inner panel: #4a5260 (darker grey-blue)
# colorbar red peak → use as accent
# colorbar blue base → use as fill
# text: white
# =============================================================================
bg_outer   = '#6b7280'   # outer figure — matches ParaView window grey
bg_panel   = '#4a5260'   # axes panel — matches ParaView inner viewport
text_color = '#ffffff'   # white text — matches ParaView labels
spine_color= '#8a9aaa'   # light grey-blue — matches ParaView border lines
ref_color  = '#c0ccd8'   # light grey for reference lines
line_color = '#d04030'   # red — matches ParaView colorbar hot end
fill_color = '#1a3a6a'   # deep blue — matches ParaView colorbar cold end
title_color= '#ffffff'

fig, ax = plt.subplots(figsize=(6, 9))
fig.patch.set_facecolor(bg_outer)
ax.set_facecolor(bg_panel)

# =============================================================================
# FILL — deep blue matching ParaView cold color
# =============================================================================
ax.fill_betweenx(z, power_MW, alpha=0.60, color=fill_color)

# =============================================================================
# MAIN LINE — red matching ParaView hot color, thick for print
# =============================================================================
ax.plot(power_MW, z, 'o-',
        color=line_color,
        linewidth=3.0,
        markersize=6,
        markerfacecolor=line_color,
        markeredgecolor=bg_panel,
        markeredgewidth=1.5)

# =============================================================================
# REFERENCE LINES
# =============================================================================
ax.axhline(y=0,
           color=ref_color, linestyle='--',
           linewidth=1.5, label='Core midplane', alpha=0.9)
ax.axhline(y=80,
           color=ref_color, linestyle=':',
           linewidth=1.2, label='Fuel boundary', alpha=0.75)
ax.axhline(y=-80,
           color=ref_color, linestyle=':',
           linewidth=1.2, alpha=0.75)

# =============================================================================
# LABELS — large for poster print
# =============================================================================
ax.set_xlabel('Power Density (MW/m³)',
              fontsize=15, color=text_color, labelpad=10)
ax.set_ylabel('Axial Position z (cm)',
              fontsize=15, color=text_color, labelpad=10)
ax.set_title('TRAVIS Core — Axial Power Distribution\n(normalized to 2.5 MWth)',
             fontsize=14, color=title_color, pad=14)

# =============================================================================
# TICKS AND SPINES
# =============================================================================
ax.tick_params(axis='both', colors=text_color,
               labelsize=13, width=1.5, length=5)
ax.grid(False)

for spine in ax.spines.values():
    spine.set_edgecolor(spine_color)
    spine.set_linewidth(2.0)

# =============================================================================
# LEGEND — matching ParaView colorbar box style
# =============================================================================
ax.legend(fontsize=12,
          facecolor='#3a4250',
          edgecolor=spine_color,
          labelcolor=text_color,
          framealpha=0.92,
          loc='upper left')

plt.tight_layout()
plt.savefig('axial_power_paraview_theme.png', dpi=300, facecolor=bg_outer)
plt.savefig('axial_power_paraview_theme.pdf', facecolor=bg_outer)
print("Saved: axial_power_paraview_theme.png and .pdf")
plt.show()
