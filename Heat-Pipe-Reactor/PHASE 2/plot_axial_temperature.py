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
z_cm   = data[:, 0]          # axial positions in cm
q_Wm3  = data[:, 1]          # power density in W/m³

# IVTBC results from Phase 3
with open('../PHASE 3/phase3_ivtbc_results.json', 'r') as f:
    ivtbc = json.load(f)

# =============================================================================
# RECONSTRUCT AXIAL TEMPERATURE PROFILE
# Tvap(z) = T_inf + (r2 × q_pipe(z))
# q_pipe(z) scales with local axial power shape
# Tm(z)   = Tvap(z) + r1 × q''(z)
# =============================================================================
T_inf  = ivtbc['T_inf_K'] - 273.15          # K → °C
r1     = ivtbc['r1_m2KW']                   # evaporator resistance m²K/W
r2     = ivtbc['r2_m2KW']                   # condenser resistance m²K/W
N_hp   = ivtbc['N_hp']
D_hp   = ivtbc['D_hp_m']
lc     = ivtbc['lc_m']
P_tot  = ivtbc['P_total_W']

# axial peaking factor: normalize power profile to its mean
q_norm = q_Wm3 / q_Wm3.mean()              # shape factor, mean = 1.0

# average heat pipe load
q_avg_W = P_tot / N_hp                     # W per pipe (average)

# axial heat pipe load varies with local power
q_pipe_z = q_avg_W * q_norm               # W per pipe at each z

# heat flux at condenser surface at each z
A_evap = np.pi * D_hp * lc                # evaporator area m²
q_flux_z = q_pipe_z / A_evap             # W/m²

# vapor temperature profile
Tvap_z_C = T_inf + r2 * (q_pipe_z / (np.pi * D_hp * lc))

# graphite surface temperature profile  
Tm_z_C = Tvap_z_C + r1 * q_flux_z

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

# plot Tvap and Tm profiles
ax.plot(Tvap_z_C, z_cm, '-o', color='#5b9bd5', linewidth=2.5,
        markersize=3.5, markerfacecolor='#5b9bd5',
        markeredgecolor=bg_color, markeredgewidth=1.5,
        label='T$_{vap}$ — vapor temperature')

ax.plot(Tm_z_C, z_cm, '-o', color='#e8833a', linewidth=2.5,
        markersize=3.5, markerfacecolor='#e8833a',
        markeredgecolor=bg_color, markeredgewidth=1.5,
        label='T$_{m}$ — graphite surface')

# fill between Tvap and Tm to show ΔT gap
ax.fill_betweenx(z_cm, Tvap_z_C, Tm_z_C,
                 alpha=0.12, color='#e8833a', label='ΔT (r₁ resistance)')

# safety limit lines
ax.axvline(x=900,  color='#e74c3c', linewidth=1.5,
           linestyle='--', label='Cladding limit 900°C')
ax.axvline(x=1600, color='#e67e22', linewidth=1.2,
           linestyle=':', label='TRISO limit 1600°C')

# midplane and fuel boundary reference lines
ax.axhline(y=0,   color='#888888', linewidth=1.2, linestyle='--')
ax.axhline(y=80,  color='#666666', linewidth=1.0, linestyle=':')
ax.axhline(y=-80, color='#666666', linewidth=1.0, linestyle=':')

# labels
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
