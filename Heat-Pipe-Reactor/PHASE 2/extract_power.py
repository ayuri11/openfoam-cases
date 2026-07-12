import numpy as np
import openmc
import json

# =============================================================================
# PHASE 2: EXTRACT PIN POWER DISTRIBUTION FROM OPENMC STATEPOINT
# reads ../PHASE 1/statepoint.100.h5 → extracts heating tally → converts to W/m³
# output feeds into OpenFOAM fvModels as heat source per mesh cell
# =============================================================================

# =============================================================================
# STEP 1: LOAD STATEPOINT AND EXTRACT TALLY
# =============================================================================
print("Loading ../PHASE 1/statepoint.100.h5...")
sp = openmc.StatePoint('../PHASE 1/statepoint.100.h5')

# get the power distribution tally defined in Phase 1
tally = sp.get_tally(name='power_distribution')
print(f"Tally loaded: {tally.name}")
print(f"Tally scores: {tally.scores}")

# =============================================================================
# STEP 2: EXTRACT RAW HEATING DATA
# heating score in OpenMC = energy deposited in eV per source neutron
# shape = (mesh_x, mesh_y, mesh_z, 1) flattened
# =============================================================================
# get mean values and standard deviations
heating_mean = tally.get_values(scores=['heating'], value='mean').flatten()
heating_std  = tally.get_values(scores=['heating'], value='std_dev').flatten()
flux_mean    = tally.get_values(scores=['flux'],    value='mean').flatten()

print(f"\nRaw heating array shape: {heating_mean.shape}")
print(f"Mesh dimensions: 50 x 50 x 28 = {50*50*28} cells")
print(f"Total heating (sum): {heating_mean.sum():.4e} eV/source-neutron")

# =============================================================================
# STEP 3: CONVERT TO PHYSICAL UNITS (W/m³)
# OpenMC heating is in eV per source neutron per cm³
# Need: target thermal power P_total (W) to normalize
# =============================================================================
# design target from methodology: ~2.5 MWth
P_total_W   = 2.5e6  # W — target thermal power
eV_to_J     = 1.60218e-19     # J per eV

# mesh geometry (must match what was defined in Phase 1)
core_radius      = 45.0        # cm
core_height      = 160.0       # cm
mesh_nx, mesh_ny, mesh_nz = 50, 50, 28

# mesh cell volume
# total mesh spans 2*core_radius × 2*core_radius × core_height
mesh_x_span = 2 * core_radius  # cm
mesh_y_span = 2 * core_radius  # cm
mesh_z_span = core_height       # cm

cell_vol_cm3 = (mesh_x_span/mesh_nx) * (mesh_y_span/mesh_ny) * (mesh_z_span/mesh_nz)
cell_vol_m3  = cell_vol_cm3 * 1e-6  # convert cm³ to m³

print(f"\nMesh cell volume: {cell_vol_cm3:.4f} cm³ = {cell_vol_m3:.6e} m³")

# normalize: heating_mean is per source neutron
# sum of heating × eV_to_J = total power per source neutron (in J/source-n)
# normalization factor = P_total / (sum_heating × eV_to_J)
total_heating_J = heating_mean.sum() * eV_to_J  # J per source neutron
norm_factor = P_total_W / total_heating_J        # source neutrons per second

print(f"Total heating per source neutron: {total_heating_J:.4e} J/source-n")
print(f"Normalization factor: {norm_factor:.4e} source-n/s")

# convert to W/m³ per mesh cell
# heating_mean[i] eV/source-n × norm_factor source-n/s × eV_to_J J/eV / cell_vol_m3 m³
q_Wm3 = heating_mean * norm_factor * eV_to_J / cell_vol_m3

print(f"\nPower density (W/m³):")
print(f"  Maximum: {q_Wm3.max():.4e} W/m³")
print(f"  Minimum (nonzero): {q_Wm3[q_Wm3>0].min():.4e} W/m³")
print(f"  Mean (nonzero): {q_Wm3[q_Wm3>0].mean():.4e} W/m³")

# verify total power integrates to ~2.5 MWth
total_power_check = (q_Wm3 * cell_vol_m3).sum()
print(f"Verification: integrated total power = {total_power_check/1e6:.4f} MWth")
print(f"  Note: normalized to 15 MWth reference; actual design power = 2.5 MWth (see Phase 3)")

# =============================================================================
# STEP 4: RESHAPE INTO 3D ARRAY (x, y, z)
# OpenMC mesh ordering: z varies slowest, x varies fastest
# =============================================================================
q_3d = q_Wm3.reshape((mesh_nz, mesh_ny, mesh_nx))
# reorder to (x, y, z) for OpenFOAM convention
q_3d_xyz = q_3d.transpose((2, 1, 0))  # now shape = (nx, ny, nz)

print(f"\n3D power density array shape (x,y,z): {q_3d_xyz.shape}")

# =============================================================================
# STEP 5: EXTRACT AXIAL PROFILE
# average radially → 1D axial power profile
# this shows the axial peaking factor for IVTBC validation
# =============================================================================
axial_profile = q_3d_xyz.mean(axis=(0,1))  # average over x,y → shape (nz,)
axial_z = np.linspace(-core_height/2, core_height/2, mesh_nz)

print(f"\nAxial power profile (W/m³):")
for i, (z, q) in enumerate(zip(axial_z, axial_profile)):
    bar = '█' * int(q / axial_profile.max() * 30)
    print(f"  z={z:+7.1f} cm: {q:.3e} W/m³  {bar}")

# =============================================================================
# STEP 6: SAVE OUTPUTS
# saves normalized power density for OpenFOAM input
# saves summary statistics for methodology documentation
# =============================================================================

# save full 3D array as numpy binary
np.save('power_density_3d.npy', q_3d_xyz)
print(f"\nSaved: power_density_3d.npy  shape={q_3d_xyz.shape}")

# save axial profile as csv
np.savetxt('axial_power_profile.csv',
           np.column_stack([axial_z, axial_profile]),
           delimiter=',',
           header='z_cm,q_Wm3',
           comments='')
print("Saved: axial_power_profile.csv")

# save summary statistics as json for methodology documentation
summary = {
    'k_effective': float(sp.keff.n),
    'k_uncertainty': float(sp.keff.s),
    'P_total_MWth_reference': P_total_W / 1e6,
    'norm_factor': float(norm_factor),
    'q_max_Wm3': float(q_Wm3.max()),
    'q_mean_Wm3': float(q_Wm3[q_Wm3>0].mean()),
    'cell_vol_cm3': cell_vol_cm3,
    'mesh_dims': [mesh_nx, mesh_ny, mesh_nz],
    'total_power_check_MWth': float(total_power_check/1e6)
}

with open('phase1_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)
print("Saved: phase1_summary.json")

# =============================================================================
# STEP 7: GENERATE OPENFOAM fvModels INPUT
# writes the heat source definition for each axial slice
# this is what gets pasted into OpenFOAM constant/fvModels file
# =============================================================================
print("\nGenerating OpenFOAM fvModels heat source entries...")

fv_lines = []
fv_lines.append("// OpenFOAM fvModels — heat source from OpenMC Phase 1")
fv_lines.append("// Generated by extract_power.py")
fv_lines.append("// Units: W/m³ normalized to 2.5 MWth reference power")
fv_lines.append("")

# write one entry per axial slice (averaged radially as starting approximation)
# in full coupling, each fuel pin × axial slice gets its own entry
for k in range(mesh_nz):
    q_slice = q_3d_xyz[:,:,k].mean()  # radial average for this axial slice
    z_lo = -core_height/2 + k * (core_height/mesh_nz)
    z_hi = z_lo + (core_height/mesh_nz)
    fv_lines.append(f"heatSource_axial_{k:02d}")
    fv_lines.append("{")
    fv_lines.append(f"    type            heatSource;")
    fv_lines.append(f"    selectionMode   cellZone;")
    fv_lines.append(f"    cellZone        fuelPin_axial_{k:02d};  // z: {z_lo:.1f} to {z_hi:.1f} cm")
    fv_lines.append(f"    q               {q_slice:.4e};  // W/m³")
    fv_lines.append("}")
    fv_lines.append("")

with open('fvModels_heatSource.txt', 'w') as f:
    f.write('\n'.join(fv_lines))
print("Saved: fvModels_heatSource.txt — paste into OpenFOAM constant/fvModels")

print("\n=== PHASE 2 COMPLETE ===")
print("Next: use power_density_3d.npy and fvModels_heatSource.txt")
print("      to set up OpenFOAM thermal-hydraulic case for IVTBC")
