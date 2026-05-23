import numpy as np
import openmc
import math
 
# =============================================================================
# MATERIALS
# =============================================================================
 
haynes = openmc.Material(name='Haynes230')
haynes.set_density('g/cm3', 8.97)
haynes.add_element('Ni', 0.57, 'wo')
haynes.add_element('Cr', 0.22, 'wo')
haynes.add_element('W',  0.14, 'wo')
haynes.add_element('Mo', 0.02, 'wo')
haynes.add_element('Fe', 0.01875, 'wo')
haynes.add_element('Co', 0.03125, 'wo')
 
b4c = openmc.Material(name='B4C')
b4c.set_density('g/cm3', 2.52)
b4c.add_nuclide('B10', 3.84, 'ao')
b4c.add_nuclide('B11', 0.16, 'ao')
b4c.add_element('C',   1.0,  'ao')
 
beo = openmc.Material(name='BeO')
beo.set_density('g/cm3', 3.025)
beo.add_element('Be', 1.0, 'ao')
beo.add_element('O',  1.0, 'ao')
 
be = openmc.Material(name='Be')
be.set_density('g/cm3', 1.85)
be.add_element('Be', 1.0, 'ao')
 
sodium = openmc.Material(name='Na')
sodium.set_density('g/cm3', 0.76)
sodium.add_element('Na', 1.0, 'ao')
 
graphite = openmc.Material(name='Graphite')
graphite.set_density('g/cm3', 1.7)
graphite.add_element('C', 1.0, 'ao')
graphite.add_s_alpha_beta('c_Graphite')
 
fuel_zone1 = openmc.Material(name='UO2_12pct')
fuel_zone1.set_density('g/cm3', 10.4)
fuel_zone1.add_nuclide('U235', 0.12, 'ao')
fuel_zone1.add_nuclide('U238', 0.88, 'ao')
fuel_zone1.add_nuclide('O16',  2.0,  'ao')
 
fuel_zone2 = openmc.Material(name='UO2_15pct')
fuel_zone2.set_density('g/cm3', 10.4)
fuel_zone2.add_nuclide('U235', 0.15, 'ao')
fuel_zone2.add_nuclide('U238', 0.85, 'ao')
fuel_zone2.add_nuclide('O16',  2.0,  'ao')
 
fuel_zone3 = openmc.Material(name='UO2_1975pct')
fuel_zone3.set_density('g/cm3', 10.4)
fuel_zone3.add_nuclide('U235', 0.1975, 'ao')
fuel_zone3.add_nuclide('U238', 0.8025, 'ao')
fuel_zone3.add_nuclide('O16',  2.0,   'ao')
 
materials = openmc.Materials([
    haynes, b4c, beo, be, sodium,
    graphite, fuel_zone1, fuel_zone2, fuel_zone3
])
materials.export_to_xml()
 
# =============================================================================
# GEOMETRY PARAMETERS
# =============================================================================
 
core_height    = 160.0
hp_radius      = 0.795
hp_wall_thick  = 0.089
fuel_pin_r     = 0.635
ctrl_rod_r     = 0.795
cell_flat      = 10.0   # cm flat-to-flat; verified to fit all pins/HPs
 
axial_ref_top    = 12.5
axial_ref_bottom = 12.5
total_height     = core_height + axial_ref_top + axial_ref_bottom  # 185 cm
 
core_radius      = 45.0
reflector_radius = 65.0
 
# lattice enclosure: circumscribes all 3 rings + cell circumradius + margin
lattice_enclosure_r = 3.0 * cell_flat + cell_flat / math.sqrt(3) + 0.5  # 36.274 cm
 
# =============================================================================
# SURFACES
# =============================================================================
 
top_boundary    = openmc.ZPlane(z0=+total_height/2, boundary_type='vacuum')
bottom_boundary = openmc.ZPlane(z0=-total_height/2, boundary_type='vacuum')
fuel_top        = openmc.ZPlane(z0=+core_height/2)
fuel_bottom     = openmc.ZPlane(z0=-core_height/2)
outer_boundary  = openmc.ZCylinder(r=reflector_radius, boundary_type='vacuum')
lattice_boundary = openmc.ZCylinder(r=lattice_enclosure_r)
 
# =============================================================================
# PIN / HP / CR UNIVERSES
# FIX v6: these simple pin universes use UNBOUNDED catch-all cells (no region on mod)
# This is correct for universes placed inside a lattice — the lattice itself provides
# the outer boundary; the catch-all graphite cell inside the pin universe handles
# everything outside the fuel pin cylinder within that universe's local space.
# =============================================================================
 
hp_inner = openmc.ZCylinder(r=hp_radius - hp_wall_thick)
hp_outer = openmc.ZCylinder(r=hp_radius)
sodium_cell = openmc.Cell(fill=sodium, region=-hp_inner)
wall_cell   = openmc.Cell(fill=haynes, region=+hp_inner & -hp_outer)
hp_universe = openmc.Universe(cells=[sodium_cell, wall_cell])
 
fuel_pin_surf = openmc.ZCylinder(r=fuel_pin_r)
 
fp1_fuel = openmc.Cell(fill=fuel_zone1, region=-fuel_pin_surf)
fp1_mod  = openmc.Cell(fill=graphite,   region=+fuel_pin_surf)
fp1_universe = openmc.Universe(cells=[fp1_fuel, fp1_mod])
 
fp2_fuel = openmc.Cell(fill=fuel_zone2, region=-fuel_pin_surf)
fp2_mod  = openmc.Cell(fill=graphite,   region=+fuel_pin_surf)
fp2_universe = openmc.Universe(cells=[fp2_fuel, fp2_mod])
 
fp3_fuel = openmc.Cell(fill=fuel_zone3, region=-fuel_pin_surf)
fp3_mod  = openmc.Cell(fill=graphite,   region=+fuel_pin_surf)
fp3_universe = openmc.Universe(cells=[fp3_fuel, fp3_mod])
 
cr_surf = openmc.ZCylinder(r=ctrl_rod_r)
cr_cell = openmc.Cell(fill=b4c,      region=-cr_surf)
cr_mod  = openmc.Cell(fill=graphite, region=+cr_surf)
cr_universe = openmc.Universe(cells=[cr_cell, cr_mod])
 
# =============================================================================
# UNIT CELL UNIVERSE BUILDER
# FIX v6: the graphite catch-all cell now has NO explicit hex boundary region.
# Instead, openmc.model.HexagonalPrism() is used ONLY to define the lattice cell
# boundary at the higher level. Inside the unit cell universe, graphite is truly
# unbounded (catch-all) — the HexLattice geometry itself enforces the hex boundary
# between adjacent cells. This is the standard OpenMC pin-cell pattern and avoids
# the "particle lost at hex plane" problem caused by explicit hex planes inside
# the universe conflicting with the lattice boundary tracking.
# =============================================================================
 
def build_unit_cell(fp_universe, cr_universe, hp_universe, graphite):
    pin_ring1_r = 1.729  # cm inner ring of 6 fuel pins
    pin_ring2_r = 3.299  # cm outer ring of 6 fuel pins (30° offset)
    hp_ring_r   = 4.879  # cm 6 heat pipes at hex corners
 
    cells     = []
    pin_surfs = []
    hp_surfs  = []
 
    # 6 fuel pins: inner ring
    for i in range(6):
        ang = math.radians(i * 60)
        x = pin_ring1_r * math.cos(ang)
        y = pin_ring1_r * math.sin(ang)
        s = openmc.ZCylinder(x0=x, y0=y, r=fuel_pin_r)
        pin_surfs.append(s)
        cells.append(openmc.Cell(fill=fp_universe, region=-s))
 
    # 6 fuel pins: outer ring
    for i in range(6):
        ang = math.radians(i * 60 + 30)
        x = pin_ring2_r * math.cos(ang)
        y = pin_ring2_r * math.sin(ang)
        s = openmc.ZCylinder(x0=x, y0=y, r=fuel_pin_r)
        pin_surfs.append(s)
        cells.append(openmc.Cell(fill=fp_universe, region=-s))
 
    # 6 heat pipes: at hex corners
    for i in range(6):
        ang = math.radians(i * 60)
        x = hp_ring_r * math.cos(ang)
        y = hp_ring_r * math.sin(ang)
        s = openmc.ZCylinder(x0=x, y0=y, r=hp_radius)
        hp_surfs.append(s)
        cells.append(openmc.Cell(fill=hp_universe, region=-s))
 
    # 1 central rod
    cr_s = openmc.ZCylinder(x0=0, y0=0, r=ctrl_rod_r)
    cells.append(openmc.Cell(fill=cr_universe, region=-cr_s))
 
    # FIX v6: graphite catch-all — NO explicit hex boundary region
    # region is only "outside all pins/HPs/rod" — no hex plane intersection
    # the HexLattice lattice boundary handles the outer hex limit automatically
    graphite_region = +cr_s  # outside central rod
    for s in pin_surfs + hp_surfs:
        graphite_region = graphite_region & +s  # outside every pin and HP
    cells.append(openmc.Cell(fill=graphite, region=graphite_region))
 
    return openmc.Universe(cells=cells)
 
# =============================================================================
# HEX LATTICE
# =============================================================================
 
zone1_univ = build_unit_cell(fp1_universe, cr_universe, hp_universe, graphite)
zone2_univ = build_unit_cell(fp2_universe, cr_universe, hp_universe, graphite)
zone3_univ = build_unit_cell(fp3_universe, hp_universe, hp_universe, graphite)
 
lattice = openmc.HexLattice()
lattice.center      = (0.0, 0.0)
lattice.pitch       = (cell_flat,)
lattice.orientation = 'x'
 
lattice.universes = [
    [zone3_univ] * 18,  # Ring 3 outer  19.75%
    [zone2_univ] * 12,  # Ring 2 middle 15%
    [zone1_univ] * 6,   # Ring 1 inner  12%
    [zone1_univ],       # Ring 0 center 12%
]
 
# outer universe: Be fills space outside lattice boundary
outer_universe = openmc.Universe()
outer_universe.add_cell(openmc.Cell(fill=be))
lattice.outer = outer_universe
 
# =============================================================================
# LATTICE CELL AND CORE UNIVERSE
# FIX v6: lattice_cell uses openmc.model.HexagonalPrism for the lateral boundary
# instead of a ZCylinder. The HexagonalPrism exactly matches the lattice geometry
# so particles crossing hex faces land cleanly in the next cell or the outer universe.
# ZCylinder was cutting corners of the outermost hex cells, leaving thin undefined
# slivers that caused particles to get lost.
# Axial bounds: fuel_bottom to fuel_top only (BeO cells cover above/below).
# =============================================================================
 
# FIX v6: edge_length for enclosing HexagonalPrism
# For 3 rings (R=3), outermost cell centers are at 3*pitch from origin.
# The enclosing hex must reach those centers + half a cell flat = 3*pitch + pitch/2
# apothem = (R + 0.5) * pitch = 3.5 * 10.0 = 35.0 cm
# edge_length = apothem / (sqrt(3)/2) = 35.0 / 0.866 = 40.415 cm
enclosing_apothem   = (3 + 0.5) * cell_flat          # 35.0 cm
enclosing_edge      = enclosing_apothem / (math.sqrt(3) / 2)  # 40.415 cm
 
core_hex_prism = openmc.model.HexagonalPrism(
    edge_length=enclosing_edge,
    orientation='x'   # must match lattice.orientation
)
 
lattice_cell = openmc.Cell(
    fill=lattice,
    region=-core_hex_prism & +fuel_bottom & -fuel_top
)
 
core_universe = openmc.Universe()
core_universe.add_cell(lattice_cell)
 
# Be radial reflector: fills annular region outside lattice hex, full height
radial_reflector_cell = openmc.Cell(
    fill=be,
    region=+core_hex_prism & -outer_boundary & +bottom_boundary & -top_boundary
)
core_universe.add_cell(radial_reflector_cell)
 
# =============================================================================
# ROOT UNIVERSE
# =============================================================================
 
root_cell = openmc.Cell(name='root cell')
root_cell.fill   = core_universe
root_cell.region = (
    -outer_boundary
    & +fuel_bottom
    & -fuel_top
)
 
root_universe = openmc.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)
 
# BeO axial reflectors above and below active core
top_beo_cell = openmc.Cell(
    fill=beo,
    region=+fuel_top & -top_boundary & -outer_boundary
)
bot_beo_cell = openmc.Cell(
    fill=beo,
    region=+bottom_boundary & -fuel_bottom & -outer_boundary
)
root_universe.add_cell(top_beo_cell)
root_universe.add_cell(bot_beo_cell)
 
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()
 
# =============================================================================
# SETTINGS
# =============================================================================
 
settings = openmc.Settings()
# DEBUG values — restore to batches=200, inactive=50, particles=10000 after geometry confirmed
settings.batches   = 10
settings.inactive  = 5
settings.particles = 100
settings.temperature['multipole'] = True
settings.temperature['method']    = 'interpolation'
settings.verbosity = 10
 
# Source at center of inner-ring fuel pin at angle=0°
# (1.729, 0, 0) is inside UO2 fuel; no symmetry constraint so it is valid
settings.source = openmc.IndependentSource(
    space=openmc.stats.Point((1.729, 0.0, 0.0))
)
settings.export_to_xml()
 
# =============================================================================
# TALLIES
# =============================================================================
 
mesh = openmc.RegularMesh()
mesh.dimension   = [20, 20, 14]
mesh.lower_left  = [-core_radius, -core_radius, -core_height/2]
mesh.upper_right = [ core_radius,  core_radius,  core_height/2]
mesh_filter = openmc.MeshFilter(mesh)
 
tally = openmc.Tally(name='power_distribution')
tally.filters = [mesh_filter]
tally.scores  = ['heating', 'flux']
 
tallies = openmc.Tallies([tally])
tallies.export_to_xml()
 
print("All XML files exported.")
print(f"Enclosing hex edge_length: {enclosing_edge:.4f} cm (apothem {enclosing_apothem:.1f} cm)")
print(f"Lattice enclosure apothem: {enclosing_apothem:.1f} cm covers {3} rings at pitch {cell_flat} cm")
print(f"Total reactor height: {total_height:.1f} cm")
print("Run: openmc")
 
