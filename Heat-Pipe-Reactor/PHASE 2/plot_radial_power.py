import numpy as np
import matplotlib.pyplot as plt

power_3d = np.load('power_density_3d.npy')   # shape (50, 50, 28)

# Take midplane slice (z index 14 = center)
midplane = power_3d[:, :, 14] / 1e6          # W/m³ → MW/m³

fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(
    midplane.T,
    origin='lower',
    extent=[-60, 60, -60, 60],               # reflector_radius = 60 cm
    cmap='inferno',
    interpolation='bilinear'
)
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Power Density (MW/m³)', fontsize=11)
ax.set_xlabel('x (cm)', fontsize=11)
ax.set_ylabel('y (cm)', fontsize=11)
ax.set_title('AYURI Core — Radial Power Distribution\nat Core Midplane (z = 0)', fontsize=12)

# Draw core boundary circle
theta = np.linspace(0, 2*np.pi, 300)
ax.plot(45*np.cos(theta), 45*np.sin(theta), 'w--', linewidth=1.2, label='Core boundary (r=45cm)')
ax.plot(60*np.cos(theta), 60*np.sin(theta), 'w:',  linewidth=1.0, label='Reflector boundary (r=60cm)')
ax.legend(fontsize=9, loc='lower right')

plt.tight_layout()
plt.savefig('radial_power_map.png', dpi=300)
plt.savefig('radial_power_map.pdf')
print("Saved: radial_power_map.png and .pdf")
plt.show()
