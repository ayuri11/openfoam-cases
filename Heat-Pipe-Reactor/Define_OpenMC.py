import numpy as np #imports NumPy as np (NumPy is the standard Python math library for arrays)
import openmc #imports OpenMC Monte Carlo neutronics library
import math
# =============================================================================
# DEFINING MATERIALS SELECTED (line 5 -77)
# =============================================================================
# REFERENCE CODE USED: U-10Mo fuel (HEU 93%), single enrichment
# OUR HPR USES: UO2 fuel, three HALEU enrichment zones
# CHANGE: entire fuel material definition

# KEEP - same cladding material as reference (Haynes 230)
haynes = openmc.Material(name='Haynes230') # creates a new material object; name= is used to label 
haynes.set_density('g/cm3', 8.97) # sets the mass density; openmc alw is CSG units
haynes.add_element('Ni', 0.57, 'wo')  # add_element: adds a natural element; wo: weight fraction
haynes.add_element('Cr', 0.22, 'wo')
haynes.add_element('W',  0.14, 'wo')
haynes.add_element('Mo', 0.02, 'wo')
haynes.add_element('Fe', 0.01875, 'wo')
haynes.add_element('Co', 0.03125, 'wo')
# weight fraction must sum to 1; Haynes230 is made up of the stated elements in %

# KEEP - same B4C control rod as reference
# CHANGE - add B-10 enrichment 
b4c = openmc.Material(name='B4C')
b4c.set_density('g/cm3', 2.52)
b4c.add_nuclide('B10', 3.84, 'ao')  # B10 is the neutron-absorbing isotope in boron; 3.84/4 boron atoms are B10 (96%)
b4c.add_nuclide('B11', 0.16, 'ao')  # 0.16/4 of boron atoms are B11 (4%)
b4c.add_element('C',   1.0,  'ao')  # 1 carbon atom
#add_nucleotide: adds a specific isotope; ao: atomic fraction

# BeO material for axial reflector
beo = openmc.Material(name='BeO')
beo.set_density('g/cm3', 3.025)
beo.add_element('Be', 1.0, 'ao') #BeO stoich 1:1
beo.add_element('O',  1.0, 'ao')

# Be material for radial
be = openmc.Material(name='Be')
be.set_density('g/cm3', 1.85)
be.add_element('Be', 1.0, 'ao')

# KEEP - same sodium coolant as reference
sodium = openmc.Material(name='Na')
sodium.set_density('g/cm3', 0.76)
sodium.add_element('Na', 1.0, 'ao')

# CHANGE: remove U-10Mo, replace with UO2 three-zone HALEU
# CHANGE: new material - graphite moderator monolith 
graphite = openmc.Material(name='Graphite')
graphite.set_density('g/cm3', 1.7)
graphite.add_element('C', 1.0, 'ao')
graphite.add_s_alpha_beta('c_Graphite')  # Thermal scattering law: The S(α,β) table captures the crystal lattice vibration physics 

# CHANGE: Zone 1 - central region, 12% enrichment (replaces U-10Mo entirely)
fuel_zone1 = openmc.Material(name='UO2_12pct')
fuel_zone1.set_density('g/cm3', 10.4)
fuel_zone1.add_nuclide('U235', 0.12, 'ao') # 12% of uranium atoms are U235
fuel_zone1.add_nuclide('U238', 0.88, 'ao') # remaining 88% are U238
fuel_zone1.add_nuclide('O16',  2.0,  'ao') # O at ratio 2:1

# CHANGE: Zone 2 - middle region, 15% enrichment
fuel_zone2 = openmc.Material(name='UO2_15pct')
fuel_zone2.set_density('g/cm3', 10.4)
fuel_zone2.add_nuclide('U235', 0.15, 'ao') 
fuel_zone2.add_nuclide('U238', 0.85, 'ao')
fuel_zone2.add_nuclide('O16',  2.0,  'ao') #UO2 stoich 

# CHANGE: Zone 3 - outer region, 19.75% HALEU enrichment
fuel_zone3 = openmc.Material(name='UO2_1975pct')
fuel_zone3.set_density('g/cm3', 10.4)
fuel_zone3.add_nuclide('U235', 0.1975, 'ao')
fuel_zone3.add_nuclide('U238', 0.8025, 'ao')
fuel_zone3.add_nuclide('O16',  2.0,   'ao')

# 3 ZONES: The three fuel zones use increasing enrichment outward to compensate for neutron leakage at the core edge.

# create collection object of the materials
materials = openmc.Materials([
    haynes, b4c, beo, be, sodium,
    graphite, fuel_zone1, fuel_zone2, fuel_zone3
])
materials.export_to_xml() # OpenMC can't read python directly; can read xml files


# =============================================================================
# DEFINING GEOMETRY - PARAMETERS (line 86 - 110)
# =============================================================================
# REFERENCE: annular cylindrical geometry, 1/8 symmetry
# AYURI HPR: hexagonal lattice of unit cells, full core (symmetry removed for debug)

# KEEP concept: define dimensions as variables first, not hardcoded
# CHANGE: all values updated to AYURI HPR specs (units: cm)

core_height    = 160.0   # cm - active fuel
hp_radius      = 0.795   # 7.95mm OD/2 = 0.795cm
hp_wall_thick  = 0.089   # Haynes 230 wall thickness 
fuel_pin_r     = 0.635   # fuel pin radius 
ctrl_rod_r     = 0.795   # same OD as heat pipe 
# cell_flat = 10.0cm verified to fit 12 fuel pins + 6 HPs + 1 central rod
# hex vertex radius = 10.0/sqrt(3) = 5.774cm
# hp_ring_r=4.879  → HP outer edge=5.674 < 5.774 ✓
# pin_ring2_r=3.299 → pin2 outer=3.934 < HP inner edge=4.084 ✓
# pin_ring1_r=1.729 → pin1 outer=2.364 < pin2 inner edge=2.664 ✓
# ctrl_rod outer=0.795 < pin1 inner=1.094 ✓
cell_flat = 10.0  # cm - unit cell flat-to-flat width; matches lattice pitch

# Axial reflector thickness 
axial_ref_top    = 12.5  # cm top BeO reflector
axial_ref_bottom = 12.5  # cm bottom BeO reflector

# Total height including reflectors: 185 cm
total_height = core_height + axial_ref_top + axial_ref_bottom 

# Radial reflector (~45cm active core radius)
core_radius      = 45.0  # defines where the graphite core ends and the Be radial reflector begins
reflector_radius = 65.0  # outer reflector boundary

# lattice enclosure radius — circumscribes all 4 rings of the hex lattice
# outermost ring centers sit 3*cell_flat from origin; add cell circumradius + margin
lattice_enclosure_r = 3.0 * cell_flat + cell_flat / math.sqrt(3) + 0.5  # = 36.274cm

# =============================================================================
# GEOMETRY - SURFACES (line 114-150)
# in OpenMC, geometry is built by defining mathematical surfaces (planes, cylinders, spheres) 
# and then combining them with boolean operators to create regions (cells)
# =============================================================================

# ZPlanes for axial boundaries 
top_boundary    = openmc.ZPlane(z0=+total_height/2, boundary_type='vacuum')
bottom_boundary = openmc.ZPlane(z0=-total_height/2, boundary_type='vacuum')
fuel_top        = openmc.ZPlane(z0=+core_height/2)
fuel_bottom     = openmc.ZPlane(z0=-core_height/2)
# vaccuum: neutrons that reach this surface escape; outermost boundaries
# reflective: neutrons bounce back; for symmetry planes

# outer boundary cylinder 
outer_boundary  = openmc.ZCylinder(r=reflector_radius, boundary_type='vacuum')

# FIX v5: symmetry planes REMOVED for full-core debug run
# The 1/12 symmetry wedge (0°-30°) had an inverted sign convention:
# all fuel pin centers (at 0°,60°,30°,90°...) sit on or outside the wedge
# boundaries, so source point rejection could not be resolved without
# flipping and rederiving the plane normals.
# Running full core (no symmetry) eliminates this issue entirely.
# Symmetry can be re-added once a clean keff is confirmed.
# NOTE: full core runs ~12x slower — use debug particle counts (100) until geometry confirmed

# lattice boundary as ZCylinder
lattice_boundary = openmc.ZCylinder(r=lattice_enclosure_r)

# =============================================================================
# GEOMETRY - UNIT CELL UNIVERSE
# =============================================================================

# single heat pipe universe 
hp_inner = openmc.ZCylinder(r=hp_radius - hp_wall_thick)  # sodium vapor core
hp_outer = openmc.ZCylinder(r=hp_radius)                  # Haynes 230 wall
sodium_cell = openmc.Cell(fill=sodium,  region=-hp_inner)
wall_cell   = openmc.Cell(fill=haynes,  region=+hp_inner & -hp_outer)
hp_universe = openmc.Universe(cells=[sodium_cell, wall_cell])
# bundles into a universe; represents 1 complete heat pipe cross-section

# single fuel pin universe
fuel_pin_surf = openmc.ZCylinder(r=fuel_pin_r)

# Zone 1 pin (central region, 12%)
fp1_fuel = openmc.Cell(fill=fuel_zone1, region=-fuel_pin_surf)
fp1_mod  = openmc.Cell(fill=graphite,   region=+fuel_pin_surf)
fp1_universe = openmc.Universe(cells=[fp1_fuel, fp1_mod])

# Zone 2 pin (middle region, 15%)
fp2_fuel = openmc.Cell(fill=fuel_zone2, region=-fuel_pin_surf)
fp2_mod  = openmc.Cell(fill=graphite,   region=+fuel_pin_surf)
fp2_universe = openmc.Universe(cells=[fp2_fuel, fp2_mod])

# Zone 3 pin (outer region, 19.75%)
fp3_fuel = openmc.Cell(fill=fuel_zone3, region=-fuel_pin_surf)
fp3_mod  = openmc.Cell(fill=graphite,   region=+fuel_pin_surf)
fp3_universe = openmc.Universe(cells=[fp3_fuel, fp3_mod])

# Control rod universe = B4C absorber + graphite
cr_surf = openmc.ZCylinder(r=ctrl_rod_r)
cr_cell = openmc.Cell(fill=b4c,      region=-cr_surf)
cr_mod  = openmc.Cell(fill=graphite, region=+cr_surf)
cr_universe = openmc.Universe(cells=[cr_cell, cr_mod])

# =============================================================================
# GEOMETRY - UNIT CELL UNIVERSE BUILDER
# =============================================================================
# Places 12 fuel pins + 6 heat pipes + 1 central rod inside a hexagonal graphite block.
# hex boundary built from 6 explicit half-space planes.
# all cylinder surfaces created once and stored; reused for graphite exclusion region.

def build_unit_cell(fp_universe, cr_universe, hp_universe, graphite):
    pin_ring1_r = 1.729  # cm - inner ring of 6 fuel pins
    pin_ring2_r = 3.299  # cm - outer ring of 6 fuel pins (offset 30° from inner)
    hp_ring_r   = 4.879  # cm - 6 heat pipes at hex corners

    cells     = []
    pin_surfs = []  # stored for graphite region reuse
    hp_surfs  = []

    # 6 fuel pins: inner ring
    for i in range(6):
        ang = math.radians(i * 60)
        x = pin_ring1_r * math.cos(ang)
        y = pin_ring1_r * math.sin(ang)
        pin_s = openmc.ZCylinder(x0=x, y0=y, r=fuel_pin_r)
        pin_surfs.append(pin_s)
        cells.append(openmc.Cell(fill=fp_universe, region=-pin_s))

    # 6 fuel pins: outer ring (30° offset)
    for i in range(6):
        ang = math.radians(i * 60 + 30)
        x = pin_ring2_r * math.cos(ang)
        y = pin_ring2_r * math.sin(ang)
        pin_s = openmc.ZCylinder(x0=x, y0=y, r=fuel_pin_r)
        pin_surfs.append(pin_s)
        cells.append(openmc.Cell(fill=fp_universe, region=-pin_s))

    # 6 heat pipes: at hex corners
    for i in range(6):
        ang = math.radians(i * 60)
        x = hp_ring_r * math.cos(ang)
        y = hp_ring_r * math.sin(ang)
        hp_s = openmc.ZCylinder(x0=x, y0=y, r=hp_radius)
        hp_surfs.append(hp_s)
        cells.append(openmc.Cell(fill=hp_universe, region=-hp_s))

    # 1 central rod
    cr_s = openmc.ZCylinder(x0=0, y0=0, r=ctrl_rod_r)
    cells.append(openmc.Cell(fill=cr_universe, region=-cr_s))

    # graphite: hex interior minus all pins/HPs/rod
    # hex boundary from 6 explicit half-space planes (apothem = cell_flat/2)
    apothem = cell_flat / 2.0
    hex_planes = []
    for k in range(6):
        theta = math.radians(k * 60)
        nx = math.cos(theta)
        ny = math.sin(theta)
        plane = openmc.Plane(a=nx, b=ny, c=0.0, d=apothem)
        hex_planes.append(plane)

    hex_interior = -hex_planes[0]
    for plane in hex_planes[1:]:
        hex_interior = hex_interior & -plane

    graphite_region = hex_interior
    for s in pin_surfs + hp_surfs:
        graphite_region = graphite_region & +s
    graphite_region = graphite_region & +cr_s

    cells.append(openmc.Cell(fill=graphite, region=graphite_region))
    return openmc.Universe(cells=cells)

# =============================================================================
# GEOMETRY - HEX LATTICE
# =============================================================================
# 37 cells = 3 rings (1 + 6 + 12 + 18 = 37)
# Ring 0/1: Zone 1 - 12% fuel + control rod center
# Ring 2:   Zone 2 - 15% fuel
# Ring 3:   Zone 3 - 19.75% fuel + extra HP center

zone1_univ = build_unit_cell(fp1_universe, cr_universe, hp_universe, graphite)
zone2_univ = build_unit_cell(fp2_universe, cr_universe, hp_universe, graphite)
zone3_univ = build_unit_cell(fp3_universe, hp_universe, hp_universe, graphite)
# zone3 central position uses hp_universe (extra heat pipe, not control rod)

lattice = openmc.HexLattice()
lattice.center      = (0.0, 0.0)
lattice.pitch       = (cell_flat,)  # center-to-center = flat-to-flat = 10.0cm
lattice.orientation = 'x'           # flat face perpendicular to x-axis

# OpenMC reads rings outermost-first
lattice.universes = [
    [zone3_univ] * 18,  # Ring 3 - outer,  19.75% HALEU
    [zone2_univ] * 12,  # Ring 2 - middle, 15%
    [zone1_univ] * 6,   # Ring 1 - inner,  12%
    [zone1_univ],       # Ring 0 - center, 12%
]

# outer universe: catches particles that cross lattice outer boundary; filled with Be
outer_universe = openmc.Universe()
outer_universe.add_cell(openmc.Cell(fill=be))
lattice.outer = outer_universe

# =============================================================================
# GEOMETRY - LATTICE CELL AND CORE UNIVERSE
# =============================================================================
# lattice_cell axially bounded by fuel_bottom/fuel_top (not vacuum boundaries)
# so it does NOT overlap with the BeO axial reflector cells below

lattice_cell = openmc.Cell(
    fill=lattice,
    region=-lattice_boundary & +fuel_bottom & -fuel_top
)

core_universe = openmc.Universe()
core_universe.add_cell(lattice_cell)

# Be radial reflector: annular region between lattice cylinder and outer vacuum
# spans full height including axial reflector zones (bottom_boundary to top_boundary)
radial_reflector_cell = openmc.Cell(
    fill=be,
    region=+lattice_boundary & -outer_boundary & +bottom_boundary & -top_boundary
)
core_universe.add_cell(radial_reflector_cell)

# =============================================================================
# GEOMETRY - ROOT CELL AND GEOMETRY EXPORT
# =============================================================================
# FIX v5: symmetry planes removed — root_cell uses full cylinder with no wedge restriction
# root_cell restricted to active fuel zone (fuel_bottom to fuel_top)
# BeO axial reflector cells sit cleanly above/below with no overlap

root_cell = openmc.Cell(name='root cell')
root_cell.fill   = core_universe
root_cell.region = (
    -outer_boundary   # inside outer vacuum cylinder (r < 65cm)
    & +fuel_bottom    # above bottom of active fuel zone
    & -fuel_top       # below top of active fuel zone
    # FIX v5: no symmetry plane conditions — full 360° core
)

root_universe = openmc.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)

# BeO axial reflector cells above and below active core
# fills 12.5cm gaps with BeO so neutrons reflect instead of streaming through void
# non-overlapping with root_cell (root_cell restricted to fuel zone only)
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

# FIX v5: source point at center of a zone1 fuel pin, guaranteed inside fuel cell
# Pin at inner ring (r=1.729cm), angle=0° → center at (1.729, 0, 0)
# FIX v5: no symmetry constraint → (1.729, 0, 0) is now valid (full 360° core)
# Previously this point was on the sym_plane_1 boundary and was rejected
settings.source = openmc.IndependentSource(
    space=openmc.stats.Point((1.729, 0.0, 0.0))
    # center of inner-ring fuel pin at 0° — inside UO2 fuel, guaranteed fissionable
    # fission source converges to true spatial distribution during inactive batches
)
settings.export_to_xml()

# =============================================================================
# TALLIES
# =============================================================================
# mesh tally covering active core; NA >= 14 axial slices
mesh = openmc.RegularMesh()
mesh.dimension = [20, 20, 14]  # 20 bins X, 20 bins Y, 14 axial slices
mesh.lower_left  = [-core_radius, -core_radius, -core_height/2]
mesh.upper_right = [ core_radius,  core_radius,  core_height/2]
mesh_filter = openmc.MeshFilter(mesh)

tally = openmc.Tally(name='power_distribution')
tally.filters = [mesh_filter]
tally.scores  = ['heating', 'flux']

tallies = openmc.Tallies([tally])
tallies.export_to_xml()

print("All XML files exported.")
print(f"Lattice enclosure radius: {lattice_enclosure_r:.3f} cm")
print(f"Total reactor height: {total_height:.1f} cm")
print("Run: openmc")  # output feeds into OpenFOAM as fvModels heat source


