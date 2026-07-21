import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load the CSV saved by Phase 2
df = pd.read_csv('axial_power_profile.csv')

# Expected columns: z_cm, power_density_W_per_m3
# Adjust column names if yours differ
z    = df.iloc[:, 0].values   # z position in cm
power = df.iloc[:, 1].values  # power density in W/m³

# Convert to more readable units
power_MW = power / 1e6        # W/m³ → MW/m³

fig, ax = plt.subplots(figsize=(5, 8))  # tall and narrow — suits axial data

ax.plot(power_MW, z, 'o-', color='#c0392b', linewidth=2, markersize=5)

# Shade under the curve for visual weight
ax.fill_betweenx(z, power_MW, alpha=0.15, color='#c0392b')

# Reference lines
ax.axhline(y=0,  color='gray', linestyle='--', linewidth=0.8, label='Core midplane')
ax.axhline(y=80, color='black', linestyle=':', linewidth=0.8, label='Fuel boundary')
ax.axhline(y=-80, color='black', linestyle=':', linewidth=0.8)

ax.set_xlabel('Power Density (MW/m³)', fontsize=12)
ax.set_ylabel('Axial Position z (cm)', fontsize=12)
ax.set_title('AYURI Core — Axial Power Distribution\n(normalized to 2.5 MWth)', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('axial_power_profile.png', dpi=300)
plt.savefig('axial_power_profile.pdf')          # vector format — best for slides/thesis
print("Saved: axial_power_profile.png and .pdf")
plt.show()
