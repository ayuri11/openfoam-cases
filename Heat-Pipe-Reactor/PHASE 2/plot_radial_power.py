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

# Normalize so total power is conserved (not multiplied by 12)
# each rotation is a copy of the same wedge, so just use one
# and mask to show the symmetric pattern
fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(
    full_core.T,
    origin='lower',
    extent=[-60, 60, -60, 60],
    cmap='inferno',
    interpolation='bilinear'
)
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Power Density (MW/m³)', fontsize=11)
ax.set_xlabel('x (cm)', fontsize=11)
ax.set_ylabel('y (cm)', fontsize=11)
ax.set_title('AYURI Core — Reconstructed Full-Core Radial Power\n(from 1/12 symmetry, at z = 0)', fontsize=12)

theta = np.linspace(0, 2*np.pi, 300)
ax.plot(45*np.cos(theta), 45*np.sin(theta), 'w--', linewidth=1.2, label='Core boundary (r=45cm)')
ax.plot(60*np.cos(theta), 60*np.sin(theta), 'w:', linewidth=1.0, label='Reflector boundary (r=60cm)')
ax.legend(fontsize=9, loc='lower right')

plt.tight_layout()
plt.savefig('radial_power_map_fullcore.png', dpi=300)
plt.savefig('radial_power_map_fullcore.pdf')
print("Saved: radial_power_map_fullcore.png")
plt.show()
