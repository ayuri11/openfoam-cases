import numpy as np
import openmc

# =============================================================================
# MATERIALS
# =============================================================================
# REFERENCE CODE USED: U-10Mo fuel (HEU 93%), single enrichment
# OUR HPR USES: UO2 fuel, three HALEU enrichment zones
# CHANGE: entire fuel material definition

# ---- KEEP FROM REFERENCE: structural materials (just update names/values) ----

# KEEP - same cladding material as reference (Haynes 230)
haynes = openmc.Material(name='Haynes230')
haynes.set_density('g/cm3', 8.97)
haynes.add_element('Ni', 0.57, 'wo')
haynes.add_element('Cr', 0.22, 'wo')
haynes.add_element('W',  0.14, 'wo')
haynes.add_element('Mo', 0.02, 'wo')
haynes.add_element('Fe', 0.01875, 'wo')
haynes.add_element('Co', 0.03125, 'wo')

# KEEP - same B4C control rod as reference
# CHANGE - add B-10 enrichment (96%) per your specs
b4c = openmc.Material(name='B4C')
b4c.set_density('g/cm3', 2.52)
b4c.add_nuclide('B10', 0.96, 'wo')  # CHANGE: was add_nuclide('B10',4,'ao') in reference
b4c.add_nuclide('B11', 0.04, 'wo')  # CHANGE: added B11 remainder
b4c.add_element('C', 1.0, 'ao')

# KEEP - same BeO reflector as reference
beo = openmc.Material(name='BeO')
beo.set_density('g/cm3', 3.025)
beo.add_element('Be', 1.0, 'ao')
beo.add_element('O',  1.0, 'ao')

# KEEP - same sodium coolant as reference
sodium = openmc.Material(name='Na')
sodium.set_density('g/cm3', 0.76)
sodium.add_element('Na', 1.0, 'ao')

# ---- CHANGE: remove U-10Mo, replace with UO2 three-zone HALEU ----

# CHANGE: new material - graphite moderator monolith (not in reference)
graphite = openmc.Material(name='Graphite')
graphite.set_density('g/cm3', 1.7)
graphite.add_element('C', 1.0, 'ao')
graphite.add_s_alpha_beta('c_Graphite')  # thermal scattering law

# CHANGE: Zone 1 - central region, 12% enrichment (replaces U-10Mo entirely)
fuel_zone1 = openmc.Material(name='UO2_12pct')
fuel_zone1.set_density('g/cm3', 10.4)
fuel_zone1.add_nuclide('U235', 0.12,  'wo')
fuel_zone1.add_nuclide('U238', 0.88,  'wo')
fuel_zone1.add_element('O',    2.0,   'ao')

# CHANGE: Zone 2 - middle region, 15% enrichment
fuel_zone2 = openmc.Material(name='UO2_15pct')
fuel_zone2.set_density('g/cm3', 10.4)
fuel_zone2.add_nuclide('U235', 0.15,  'wo')
fuel_zone2.add_nuclide('U238', 0.85,  'wo')
fuel_zone2.add_element('O',    2.0,   'ao')

# CHANGE: Zone 3 - outer region, 19.75% HALEU enrichment
fuel_zone3 = openmc.Material(name='UO2_1975pct')
fuel_zone3.set_density('g/cm3', 10.4)
fuel_zone3.add_nuclide('U235', 0.1975, 'wo')
fuel_zone3.add_nuclide('U238', 0.8025, 'wo')
fuel_zone3.add_element('O',    2.0,    'ao')

# CHANGE: updated materials list (removed U-10Mo, added graphite + 3 fuel zones)
materials = openmc.Materials([
    haynes, b4c, beo, sodium,
    graphite, fuel_zone1, fuel_zone2, fuel_zone3
])
materials.export_to_xml()


# =============================================================================
# GEOMETRY - PARAMETERS
# =============================================================================
# REFERENCE: annular cylindrical geometry, 1/8 symmetry
# AYURI HPR: hexagonal lattice of unit cells, 1/12 symmetry
# CHANGE: all geometry below is new - reference geometry does not apply

# KEEP concept: define dimensions as variables first, not hardcoded
# CHANGE: all values updated to your HPR specs (units: cm)

core_height    = 160.0   # CHANGE: was fuel_h in reference, now 160cm per your specs
hp_radius      = 0.795   # KEEP concept, CHANGE value: 7.95mm OD/2 = 0.795cm
hp_wall_thick  = 0.089   # NEW: Haynes 230 wall thickness ~0.9mm
fuel_pin_r     = 0.635   # NEW: fuel pin radius (to be confirmed from your refs)
ctrl_rod_r     = 0.795   # NEW: same OD as heat pipe per your specs
cell_flat      = 5.5     # CHANGE: unit cell flat-to-flat 55mm = 5.5cm

# Axial reflector thickness per your methodology outline
axial_ref_top    = 1.25  # NEW: 12.5cm top BeO reflector
axial_ref_bottom = 1.25  # NEW: 12.5cm bottom BeO reflector

# Total height including reflectors
total_height = core_height + axial_ref_top + axial_ref_bottom

# Radial reflector per your specs (~45cm active core radius)
core_radius      = 45.0  # CHANGE value
reflector_radius = 65.0  # NEW: outer reflector boundary


# =============================================================================
# GEOMETRY - SURFACES
# =============================================================================
# REFERENCE: ZCylinder rings + ZPlanes for annular geometry
# AYURI HPR: hexagonal prism surfaces for unit cells
# CHANGE: replace ZCylinder rings with hexagonal surfaces

# KEEP: ZPlanes for axial boundaries (same concept as reference)
top_boundary    = openmc.ZPlane(z0=+total_height/2, boundary_type='vacuum')
bottom_boundary = openmc.ZPlane(z0=-total_height/2, boundary_type='vacuum')
fuel_top        = openmc.ZPlane(z0=+core_height/2)
fuel_bottom     = openmc.ZPlane(z0=-core_height/2)
top_ref_plane   = openmc.ZPlane(z0=+core_height/2)   # BeO axial reflector start
bot_ref_plane   = openmc.ZPlane(z0=-core_height/2)   # BeO axial reflector start

# KEEP: outer boundary cylinder (same concept as reference reflector_OD)
outer_boundary  = openmc.ZCylinder(r=reflector_radius, boundary_type='vacuum')

# CHANGE: symmetry planes for 1/12 hex (reference used 1/8 with 2 planes)
# 1/12 of hexagon = 30 degree slice
# Two planes at 0 deg and 30 deg from x-axis
import math
angle1 = 0.0                  # 0 degrees
angle2 = math.radians(30.0)   # 30 degrees = 1/12 of 360
sym_plane_1 = openmc.Plane(
    a=math.sin(angle1), b=-math.cos(angle1), c=0, d=0,
    boundary_type='reflective'
)
sym_plane_2 = openmc.Plane(
    a=math.sin(angle2), b=-math.cos(angle2), c=0, d=0,
    boundary_type='reflective'
)


# =============================================================================
# GEOMETRY - UNIT CELL UNIVERSE
# =============================================================================
# REFERENCE: pin_cell_universe with annular fuel rings
# AYURI HPR: hexagonal unit cell with 12 fuel pins + 6 HPs + 1 central rod
# CHANGE: completely new unit cell definition for each zone

# --- Single heat pipe universe (KEEP concept from reference, same structure) ---
hp_inner = openmc.ZCylinder(r=hp_radius - hp_wall_thick)  # sodium vapor core
hp_outer = openmc.ZCylinder(r=hp_radius)                  # Haynes 230 wall

sodium_cell = openmc.Cell(fill=sodium,  region=-hp_inner)
wall_cell   = openmc.Cell(fill=haynes,  region=+hp_inner & -hp_outer)
hp_universe = openmc.Universe(cells=[sodium_cell, wall_cell])

# --- Single fuel pin universe (CHANGE: UO2 replaces U-10Mo) ---
fuel_pin_surf = openmc.ZCylinder(r=fuel_pin_r)

# Zone 1 pin (central region)
fp1_fuel = openmc.Cell(fill=fuel_zone1, region=-fuel_pin_surf)
fp1_mod  = openmc.Cell(fill=graphite,   region=+fuel_pin_surf)
fp1_universe = openmc.Universe(cells=[fp1_fuel, fp1_mod])

# Zone 2 pin (middle region)
fp2_fuel = openmc.Cell(fill=fuel_zone2, region=-fuel_pin_surf)
fp2_mod  = openmc.Cell(fill=graphite,   region=+fuel_pin_surf)
fp2_universe = openmc.Universe(cells=[fp2_fuel, fp2_mod])

# Zone 3 pin (outer region)
fp3_fuel = openmc.Cell(fill=fuel_zone3, region=-fuel_pin_surf)
fp3_mod  = openmc.Cell(fill=graphite,   region=+fuel_pin_surf)
fp3_universe = openmc.Universe(cells=[fp3_fuel, fp3_mod])

# --- Control rod universe (KEEP concept, same material) ---
cr_surf = openmc.ZCylinder(r=ctrl_rod_r)
cr_cell = openmc.Cell(fill=b4c,      region=-cr_surf)
cr_mod  = openmc.Cell(fill=graphite, region=+cr_surf)
cr_universe = openmc.Universe(cells=[cr_cell, cr_mod])


# =============================================================================
# GEOMETRY - HEX LATTICE
# =============================================================================
# REFERENCE: no lattice - single pin with angular symmetry
# AYURI HPR: HexLattice of 37 unit cells
# CHANGE: entirely new section - does not exist in reference code

# ADD: define the hex lattice for the full core
# 37 cells = 3 rings (R=3: 1 + 6 + 12 + 18 = 37)
# Ring 0 (center, 1 cell): Zone 1 - control rod center, 12% fuel
# Ring 1 (6 cells): Zone 1 - control rod center, 12% fuel
# Ring 2 (12 cells): Zone 2 - extra HP center, 15% fuel
# Ring 3 (18 cells): Zone 3 - extra HP center, 19.75% fuel

# NOTE: each universe here represents one full unit cell
# For now, pin universes are used as placeholders
# Full unit cell universes with all 12 pins + 6 HPs to be built next

# Placeholder: zone universes (replace with full unit cell universes later)
zone1_cell = openmc.Cell(fill=fuel_zone1)
zone1_univ = openmc.Universe(cells=[zone1_cell])

zone2_cell = openmc.Cell(fill=fuel_zone2)
zone2_univ = openmc.Universe(cells=[zone2_cell])

zone3_cell = openmc.Cell(fill=fuel_zone3)
zone3_univ = openmc.Universe(cells=[zone3_cell])

# ADD: HexLattice definition
lattice = openmc.HexLattice()
lattice.center = (0.0, 0.0)
lattice.pitch  = (cell_flat,)           # flat-to-flat pitch in cm
lattice.orientation = 'x'              # flat side faces x-axis

# Ring arrangement: outermost ring first in OpenMC HexLattice
# Ring 3 (18 cells) = zone3, Ring 2 (12 cells) = zone2,
# Ring 1 (6 cells) = zone1, Ring 0 (1 cell) = zone1
lattice.universes = [
    [zone3_univ] * 18,   # CHANGE: outer ring - 19.75% HALEU
    [zone2_univ] * 12,   # CHANGE: middle ring - 15%
    [zone1_univ] * 6,    # CHANGE: inner ring - 12%
    [zone1_univ],        # CHANGE: center cell - 12%
]

# ADD: fill the lattice into a containing cell
lattice_cell = openmc.Cell(fill=lattice)
core_universe = openmc.Universe(cells=[lattice_cell])


# =============================================================================
# GEOMETRY - ROOT CELL AND GEOMETRY EXPORT
# =============================================================================
# KEEP: root cell with boundary conditions
# CHANGE: use 1/12 symmetry planes instead of 1/8

root_cell = openmc.Cell(name='root cell')
root_cell.fill   = core_universe
root_cell.region = (
    -outer_boundary
    & +bottom_boundary
    & -top_boundary
    & +sym_plane_1      # CHANGE: 1/12 symmetry 
    & -sym_plane_2
)

root_universe = openmc.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)

geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()


# =============================================================================
# SETTINGS
# =============================================================================
# KEEP: same settings structure
# CHANGE: source point moved to center of hex core

settings = openmc.Settings()
settings.batches   = 100    # CHANGE: increase later for production run
settings.inactive  = 20     # KEEP: same as reference
settings.particles = 1000   # CHANGE: increase to 10000+ for production run
settings.temperature['multipole'] = True
settings.temperature['method']    = 'interpolation'

# CHANGE: source point at center of hex core (was at fuel_r offset in reference)
settings.source = openmc.IndependentSource(
    space=openmc.stats.Point((0, 0, 0))
)
settings.export_to_xml()


# =============================================================================
# TALLIES
# =============================================================================
# REFERENCE: DistribcellFilter per fuel cell
# AYURI HPR: mesh tally for spatial power distribution
# CHANGE: replace per-cell tallies with regular mesh tally

# ADD: cylindrical mesh tally covering active core
# NA >= 14 axial slices per your methodology requirement
mesh = openmc.RegularMesh()
mesh.dimension = [1, 1, 14]                      # CHANGE: 14 axial slices minimum
mesh.lower_left  = [-core_radius, -core_radius, -core_height/2]
mesh.upper_right = [ core_radius,  core_radius,  core_height/2]

mesh_filter = openmc.MeshFilter(mesh)

# KEEP concept: heating and flux tallies (same scores as reference)
tally = openmc.Tally(name='power_distribution')
tally.filters = [mesh_filter]
tally.scores  = ['heating', 'flux']   # KEEP: same as reference

tallies = openmc.Tallies([tally])
tallies.export_to_xml()

print("All XML files exported. Run: openmc")
