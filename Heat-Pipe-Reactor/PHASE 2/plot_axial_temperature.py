# plot_axial_temperature.py
# Plots estimated axial temperature profile of heat pipe vapor and
# graphite surface using Phase 2 axial power shape + Phase 3 IVTBC results

import numpy as np
import matplotlib.pyplot as plt
import json

# =============================================================================
# LOAD DATA
# =============================================================================
# axial power profile from Phase 2
data = np.loadtxt('axial_power_profile.csv', delimiter=',', skiprows=1)
z_cm  = data[:, 0]     # axial positions in cm
q_Wm3 = data[:, 1]     # power density in W/m³

# IVTBC results from Phase 3 — keys matched to actual JSON output
with open('../PHASE 3/phase3_ivtbc_results.json', 'r') as f:
    ivtbc = json.load(f)

T_inf_C = ivtbc['T_inf_C']                        # °C
r1      = ivtbc['r1_m2KW']                        # m²K/W evaporator resistance
r2      = ivtbc['r2_m2KW']                        # m²K/W condenser resistance
N_hp    = ivtbc['N_heat_pipes']
P_tot_W = ivtbc['P_total_MWth'] * 1e6             # MWth → W

# reference values from ring_0 (all rings identical in uniform assumption)
q_hp_avg_W  = ivtbc['results']['ring_0']['q_hp_kW'] * 1000   # kW → W
q_flux_avg  = ivtbc['results']['ring_0']['q_flux']            # W/m²
Tvap_avg_C  = ivtbc['results']['ring_0']['Tvap_C']            # °C
Tm_avg_C    = ivtbc['results']['ring_0']['Tm_C']              # °C

# =============================================================================
# RECONSTRUCT AXIAL TEMPERATURE PROFILE
# Scale from the uniform average using the axial power shape
# Tvap(z) and Tm(z) vary proportionally with local power relative to mean
# =============================================================================
q_norm    = q_Wm3 / q_Wm3.mean()         # axial shape factor, mean = 1.0

# axial vapor temperature: base T_inf + ΔT scaled by local power
dT_vap_avg = Tvap_avg_C - T_inf_C        # average ΔT from coolant to vapor
Tvap_z_C   = T_inf_C + dT_vap_avg * q_norm

# axial graphite surface temperature: Tvap + r1*q'' scaled by local power
dT_tm_avg  = Tm_avg_C - Tvap_avg_C       # average ΔT across evaporator wall
Tm_z_C     = Tvap_z_C + dT_tm_avg * q_norm

# =============================================================================
# DARK THEME PLOT
# =============================================================================
bg_color    = '#2e2e2e'
text_color  = '#cccccc'
spine_color = '#666666'
title_color = '#ffffff'

fig, ax = plt.subplots(figsize=(5, 8))
fig.patch.set_facecolor(bg_color)
ax.set_facecolor(bg_color)

# Tvap axial profile — blue matching axial power plot
ax.plot(Tvap_z_C, z_cm, '-o', color='#5b9bd5', linewidth=2.5,
        markersize=3.5, markerfacecolor='#5b9bd5',
        markeredgecolor=bg_color, markeredgewidth=1.5,
        label='$T_{vap}$ — vapor temperature')

# Tm axial profile — orange for contrast
ax.plot(Tm_z_C, z_cm, '-o', color='#e8833a', linewidth=2.5,
        markersize=3.5, markerfacecolor='#e8833a',
        markeredgecolor=bg_color, markeredgewidth=1.5,
        label='$T_{m}$ — graphite surface')

# shaded gap between Tvap and Tm showing r1 thermal resistance
ax.fill_betweenx(z_cm, Tvap_z_C, Tm_z_C,
                 alpha=0.12, color='#e8833a',
                 label='ΔT across evaporator wall')

# safety limit vertical lines
ax.axvline(x=900,  color='#e74c3c', linewidth=1.5,
           linestyle='--', label='Cladding limit 900°C')
ax.axvline(x=1600, color='#e67e22', linewidth=1.2,
           linestyle=':', label='TRISO limit 1600°C')

# reference horizontal lines — midplane and fuel boundaries
ax.axhline(y=0,   color='#888888', linewidth=1.2,
           linestyle='--', alpha=0.6)
ax.axhline(y=80,  color='#666666', linewidth=1.0,
           linestyle=':', alpha=0.6)
ax.axhline(y=-80, color='#666666', linewidth=1.0,
           linestyle=':', alpha=0.6)

# annotate peak values at midplane
ax.annotate(f'{Tvap_z_C[len(z_cm)//2]:.1f}°C',
            xy=(Tvap_z_C[len(z_cm)//2], z_cm[len(z_cm)//2]),
            xytext=(Tvap_z_C[len(z_cm)//2] - 80, z_cm[len(z_cm)//2] + 8),
            color='#5b9bd5', fontsize=9,
            arrowprops=dict(arrowstyle='->', color='#5b9bd5', lw=1.0))

ax.annotate(f'{Tm_z_C[len(z_cm)//2]:.1f}°C',
            xy=(Tm_z_C[len(z_cm)//2], z_cm[len(z_cm)//2]),
            xytext=(Tm_z_C[len(z_cm)//2] + 20, z_cm[len(z_cm)//2] + 8),
            color='#e8833a', fontsize=9,
            arrowprops=dict(arrowstyle='->', color='#e8833a', lw=1.0))

# axes labels and title
ax.set_xlabel('Temperature (°C)', fontsize=12, color=text_color, labelpad=8)
ax.set_ylabel('Axial position z (cm)', fontsize=12, color=text_color, labelpad=8)
ax.set_title('AYURI Core — Axial Temperature Profile\n(IVTBC, normalized to 2.5 MWth)',
             fontsize=12, color=title_color, pad=12)

# ticks and spines
ax.tick_params(axis='both', colors=text_color, labelsize=10)
ax.grid(False)
for spine in ax.spines.values():
    spine.set_edgecolor(spine_color)
    spine.set_linewidth(1.5)

# legend
ax.legend(fontsize=9, loc='lower right',
          facecolor='#3a3a3a', edgecolor='#555555',
          labelcolor=text_color, framealpha=0.85)

plt.tight_layout()
plt.savefig('axial_temperature_profile_dark.png', dpi=300, facecolor=bg_color)
plt.savefig('axial_temperature_profile_dark.pdf', facecolor=bg_color)
print("Saved: axial_temperature_profile_dark.png and .pdf")
plt.show()
