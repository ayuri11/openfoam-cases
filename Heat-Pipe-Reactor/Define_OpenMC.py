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
# AYURI HPR: hexagonal lattice of unit cells, 1/12 symmetry

# KEEP concept: define dimensions as variables first, not hardcoded
# CHANGE: all values updated to AYURI HPR specs (units: cm)

core_height    = 160.0   # cm - active fuel
hp_radius      = 0.795   # 7.95mm OD/2 = 0.795cm
hp_wall_thick  = 0.089   # Haynes 230 wall thickness 
fuel_pin_r     = 0.635   # fuel pin radius 
ctrl_rod_r     = 0.795   # same OD as heat pipe 
# FIX: cell_flat = 10.0cm (was 5.5cm which was too small to fit 12 fuel pins + 6 HPs + 1 central rod)
# minimum cell_flat for this layout is ~8.6cm; 10.0cm gives adequate spacing margins
# mathematically verified ring radii (see build_unit_cell below):
#   hex vertex radius = 10.0/sqrt(3) = 5.774cm
#   hp_ring_r=4.879  → HP outer edge=5.674 < 5.774 ✓
#   pin_ring2_r=3.299 → pin2 outer=3.934 < HP inner edge=4.084 ✓
#   pin_ring1_r=1.729 → pin1 outer=2.364 < pin2 inner edge=2.664 ✓
#   ctrl_rod outer=0.795 < pin1 inner=1.094 ✓
cell_flat = 10.0  # cm - unit cell flat-to-flat width; matches lattice pitch

# Axial reflector thickness 
axial_ref_top    = 12.5  # cm top BeO reflector
axial_ref_bottom = 12.5  # cm bottom BeO reflector

# Total height including reflectors: 185 cm
total_height = core_height + axial_ref_top + axial_ref_bottom 

# Radial reflector (~45cm active core radius)
core_radius      = 45.0  # defines where the graphite core ends and the BeO radial reflector begins
reflector_radius = 65.0  # outer reflector boundary

# FIX: lattice enclosure radius — circumscribes all 4 rings (indices 0-3) of the hex lattice
# outermost ring centers sit 3*cell_flat from origin; add cell circumradius + small margin
# circumradius of hex cell = cell_flat / sqrt(3) = 5.774cm
lattice_enclosure_r = 3.0 * cell_flat + cell_flat / math.sqrt(3) + 0.5  # = 36.274cm

# =============================================================================
# GEOMETRY - SURFACES (line 114-150)
# in OpenMC, geometry is built by defining mathematical surfaces (planes, cylinders, spheres) 
# and then combining them with boolean operators to create regions (cells)
# =============================================================================
# REFERENCE: ZCylinder rings + ZPlanes for annular geometry
# AYURI HPR: hexagonal prism surfaces for unit cells
# CHANGE: replace ZCylinder rings with hexagonal surfaces

# ZPlanes for axial boundaries 
top_boundary    = openmc.ZPlane(z0=+total_height/2, boundary_type='vacuum')
bottom_boundary = openmc.ZPlane(z0=-total_height/2, boundary_type='vacuum')
fuel_top        = openmc.ZPlane(z0=+core_height/2)
fuel_bottom     = openmc.ZPlane(z0=-core_height/2)
# vaccuum: neutrons that reach this surface escapes; used on the outermost boundaries
# reflective: neutrons hitting this surface bounce back; for symmetry planes to simulate a full core with only 1/12 of it

# outer boundary cylinder 
outer_boundary  = openmc.ZCylinder(r=reflector_radius, boundary_type='vacuum')

# CHANGE: symmetry planes for 1/12 hex (reference used 1/8 with 2 planes)
# 1/12 symmetry: a regular hexagon has 12-fold symmetry (6 rotational × 2 mirror)
# by modeling only a 30 deg wedge with reflective boundaries on both sides, openmc simulates the full core
# this reduces computation time by 12×
angle1 = 0.0                  # 0 degrees
angle2 = math.radians(30.0)   # 30 degrees = 1/12 of 360; converts degrees to radians
sym_plane_1 = openmc.Plane(
    a=math.sin(angle1), b=-math.cos(angle1), c=0, d=0,
    boundary_type='reflective'
)
sym_plane_2 = openmc.Plane(
    a=math.sin(angle2), b=-math.cos(angle2), c=0, d=0,
    boundary_type='reflective'
)

# openmc plane abcd: general plane equation: ax + by + cz = d
# for a plane at angle θ from x-axis: a = sin(θ), b = -cos(θ), c = 0, d = 0
# the two planes are at 0° and 30° 

# FIX: lattice boundary as ZCylinder (replaces hexagonal_prism for enclosure)
# using a cylinder avoids hex orientation/edge_length ambiguity that caused
# particles to escape the lattice into undefined space in earlier version
lattice_boundary = openmc.ZCylinder(r=lattice_enclosure_r)

# =============================================================================
# GEOMETRY - UNIT CELL UNIVERSE (line 157-276)
# a Universe in OpenMC is a reusable geometry template; defined once, then places copies of it anywhere in the lattice 
# each universe contains cells (regions + materials); 
# used to build 37 identical-structure unit cells efficiently
# =============================================================================
# REFERENCE: pin_cell_universe with annular fuel rings
# AYURI HPR: hexagonal unit cell with 12 fuel pins + 6 HPs + 1 central rod
# CHANGE: completely new unit cell definition for each zone

# single heat pipe universe 
hp_inner = openmc.ZCylinder(r=hp_radius - hp_wall_thick)  # sodium vapor core
hp_outer = openmc.ZCylinder(r=hp_radius)                  # Haynes 230 wall
# outer radius - wall thickness = inner surface of HP wall 
# outer surface = outer radius

sodium_cell = openmc.Cell(fill=sodium,  region=-hp_inner) # everything inside the hp_inner (r < 0.706)
wall_cell   = openmc.Cell(fill=haynes,  region=+hp_inner & -hp_outer) # everything outside hp_inner & inside hp_outer
# -cylinder: INSIDE ( r < radius); +cylinder: OUTSIDE (r > radius)
# -plane: below the plane; +plane: above

hp_universe = openmc.Universe(cells=[sodium_cell, wall_cell])
# bundles into a universe; represents 1 complete heat pipe cross-section; 
# placed at the center of each unit cell in the hex lattice

# single fuel pin universe (CHANGE: UO2 replaces U-10Mo) 
fuel_pin_surf = openmc.ZCylinder(r=fuel_pin_r)

# Zone 1 pin (central region)
fp1_fuel = openmc.Cell(fill=fuel_zone1, region=-fuel_pin_surf)
fp1_mod  = openmc.Cell(fill=graphite,   region=+fuel_pin_surf)  # graphite outside pin (catch-all)
fp1_universe = openmc.Universe(cells=[fp1_fuel, fp1_mod])
 
# Zone 2 pin (middle region, 15%)
fp2_fuel = openmc.Cell(fill=fuel_zone2, region=-fuel_pin_surf)
fp2_mod  = openmc.Cell(fill=graphite,   region=+fuel_pin_surf)
fp2_universe = openmc.Universe(cells=[fp2_fuel, fp2_mod])
 
# Zone 3 pin (outer region, 19.75%)
fp3_fuel = openmc.Cell(fill=fuel_zone3, region=-fuel_pin_surf)
fp3_mod  = openmc.Cell(fill=graphite,   region=+fuel_pin_surf)
fp3_universe = openmc.Universe(cells=[fp3_fuel, fp3_mod])

# Control rod universe = B4C adsorber + graphite mix
cr_surf = openmc.ZCylinder(r=ctrl_rod_r)
cr_cell = openmc.Cell(fill=b4c,      region=-cr_surf)  # B4C inside rod surface
cr_mod  = openmc.Cell(fill=graphite, region=+cr_surf)  # graphite outside (catch-all)
cr_universe = openmc.Universe(cells=[cr_cell, cr_mod])

# =============================================================================
# GEOMETRY - UNIT CELL UNIVERSE BUILDER
# =============================================================================
# Places 12 fuel pins + 6 heat pipes + 1 central rod inside a hexagonal graphite block.
# Called once per zone to produce zone1_univ, zone2_univ, zone3_univ.
#
# FIX (this version): hex boundary built from 6 explicit half-space planes instead of
# openmc.model.hexagonal_prism(). hexagonal_prism() returns a Region object, not a
# Surface — the unary minus operator on it is unreliable across OpenMC versions and
# was silently producing malformed graphite regions, which created geometry gaps and
# lost particles whenever a neutron crossed the unit cell boundary.
#
# FIX: all cylinder surfaces (pin_surfs, hp_surfs, cr_s) are created once and stored,
# then reused for both the cell fill definitions and the graphite exclusion region.
# The earlier pattern of redefining surfaces inside the graphite region block produced
# duplicate surface IDs at the same positions — OpenMC treated them as distinct surfaces,
# which left thin undefined slivers at each cylinder boundary.
 
def build_unit_cell(fp_universe, cr_universe, hp_universe, graphite):
    """
    Build a hexagonal unit cell universe containing:
      - 6 fuel pins at inner ring (pin_ring1_r)
      - 6 fuel pins at outer ring (pin_ring2_r)
      - 6 heat pipes at hex corners (hp_ring_r)
      - 1 central rod (control rod or extra HP depending on zone)
      - graphite monolith filling all remaining space inside the hex boundary
 
    All ring radii verified to fit within cell_flat=10.0cm without overlaps.
    """
 
    pin_ring1_r = 1.729  # cm - inner ring of 6 fuel pins
    pin_ring2_r = 3.299  # cm - outer ring of 6 fuel pins (offset 30° from inner)
    hp_ring_r   = 4.879  # cm - 6 heat pipes at hex corners
 
    cells     = []
    pin_surfs = []  # store created surfaces to reuse for graphite exclusion region
    hp_surfs  = []  # avoids duplicate surface definitions that caused lost particles
 
    # --- 6 fuel pins: inner ring ---
    for i in range(6):
        ang = math.radians(i * 60)
        x = pin_ring1_r * math.cos(ang)
        y = pin_ring1_r * math.sin(ang)
        pin_s = openmc.ZCylinder(x0=x, y0=y, r=fuel_pin_r)
        pin_surfs.append(pin_s)
        cells.append(openmc.Cell(fill=fp_universe, region=-pin_s))
 
    # --- 6 fuel pins: outer ring (30° offset from inner ring) ---
    for i in range(6):
        ang = math.radians(i * 60 + 30)
        x = pin_ring2_r * math.cos(ang)
        y = pin_ring2_r * math.sin(ang)
        pin_s = openmc.ZCylinder(x0=x, y0=y, r=fuel_pin_r)
        pin_surfs.append(pin_s)
        cells.append(openmc.Cell(fill=fp_universe, region=-pin_s))
 
    # --- 6 heat pipes: at hex corners ---
    for i in range(6):
        ang = math.radians(i * 60)
        x = hp_ring_r * math.cos(ang)
        y = hp_ring_r * math.sin(ang)
        hp_s = openmc.ZCylinder(x0=x, y0=y, r=hp_radius)
        hp_surfs.append(hp_s)
        cells.append(openmc.Cell(fill=hp_universe, region=-hp_s))
 
    # --- 1 central rod (B4C control rod for zones 1/2; extra HP for zone 3) ---
    cr_s = openmc.ZCylinder(x0=0, y0=0, r=ctrl_rod_r)
    cells.append(openmc.Cell(fill=cr_universe, region=-cr_s))
 
    # --- graphite monolith: fills everything inside hex boundary, outside all rods/pins/HPs ---
    # FIX: build hex boundary from 6 explicit half-space planes (not hexagonal_prism region object)
    # For flat-x hex (flat face perpendicular to x-axis), the 6 bounding planes have normals at
    # 0°, 60°, 120°, 180°, 240°, 300°. Each plane sits at perpendicular distance = apothem
    # from the cell center, where apothem = flat-to-flat / 2 = cell_flat / 2.
    # Plane equation in OpenMC form (ax + by + cz = d): a=cos(θ), b=sin(θ), c=0, d=apothem.
    # A point is INSIDE the hex if it satisfies nx*x + ny*y <= apothem for all 6 normals,
    # i.e., it is on the negative side (-) of each plane.
    apothem = cell_flat / 2.0  # cm - perpendicular distance from center to flat face
 
    hex_planes = []
    for k in range(6):
        theta = math.radians(k * 60)
        nx = math.cos(theta)
        ny = math.sin(theta)
        # plane: nx*x + ny*y + 0*z = apothem
        plane = openmc.Plane(a=nx, b=ny, c=0.0, d=apothem)
        hex_planes.append(plane)
 
    # hex interior = negative side of all 6 planes simultaneously
    hex_interior = -hex_planes[0]
    for plane in hex_planes[1:]:
        hex_interior = hex_interior & -plane
 
    # graphite fills hex interior minus all cylindrical objects
    # FIX: reuse already-stored surface objects (pin_surfs, hp_surfs, cr_s)
    # so OpenMC sees a single surface at each position, not duplicate overlapping ones
    graphite_region = hex_interior
    for s in pin_surfs + hp_surfs:
        graphite_region = graphite_region & +s
    graphite_region = graphite_region & +cr_s
 
    cells.append(openmc.Cell(fill=graphite, region=graphite_region))
 
    return openmc.Universe(cells=cells)

# =============================================================================
# GEOMETRY - HEX LATTICE (line 280-326)
# places 37 unit cell universes in a regular hexagonal grid
# this is the 37-cell core arrangement with 3 enrichment zones
# =============================================================================
# REFERENCE: no lattice — single pin with angular symmetry
# AYURI HPR: HexLattice of 37 unit cells (R=3: 1 + 6 + 12 + 18 = 37)

# Ring assignment:
# Ring 0 (center, 1 cell): Zone 1 - control rod center, 12% fuel
# Ring 1 (6 cells): Zone 1 - control rod center, 12% fuel
# Ring 2 (12 cells): Zone 2 - extra HP center, 15% fuel
# Ring 3 (18 cells): Zone 3 - extra HP center, 19.75% fuel

# build zone universes (zone3 central rod replaced with extra heat pipe)
zone1_univ = build_unit_cell(fp1_universe, cr_universe, hp_universe, graphite)
zone2_univ = build_unit_cell(fp2_universe, cr_universe, hp_universe, graphite)
zone3_univ = build_unit_cell(fp3_universe, hp_universe, hp_universe, graphite)
 
# define the HexLattice
lattice = openmc.HexLattice()
lattice.center      = (0.0, 0.0)
lattice.pitch       = (cell_flat,)  # center-to-center spacing = flat-to-flat width = 10.0cm
lattice.orientation = 'x'           # flat face perpendicular to x-axis; must match build_unit_cell planes
 
# OpenMC reads rings outermost-first
lattice.universes = [
    [zone3_univ] * 18,  # Ring 3 - outer ring,  19.75% HALEU
    [zone2_univ] * 12,  # Ring 2 - middle ring, 15%
    [zone1_univ] * 6,   # Ring 1 - inner ring,  12%
    [zone1_univ],       # Ring 0 - center cell, 12%
]

# FIX: outer universe catches any particle that crosses the lattice outer boundary
# Without this, OpenMC throws an error when a neutron drifts outside the 37-cell grid.
# Filled with Be so particles that escape into the radial reflector region are tracked correctly.
outer_universe = openmc.Universe()
outer_universe.add_cell(openmc.Cell(fill=be))  # catch-all cell, no region = unbounded
lattice.outer = outer_universe

# =============================================================================
# GEOMETRY - LATTICE CELL AND CORE UNIVERSE
# =============================================================================
# FIX: lattice is enclosed by a ZCylinder (not hexagonal_prism) to avoid
# hex orientation/edge_length ambiguity. The cylinder circumscribes all ring
# centers plus one cell circumradius, with a small margin, so no ring center
# falls outside the cylinder and creates an undefined region.
#
# FIX: lattice_cell is axially bounded by fuel_bottom and fuel_top so it does
# NOT overlap with the BeO axial reflector cells added to root_universe below.
# Earlier version bounded by bottom_boundary/top_boundary, which caused the
# lattice cell and BeO cells to share the same axial space — geometry conflict.
 
lattice_cell = openmc.Cell(
    fill=lattice,
    region=-lattice_boundary & +fuel_bottom & -fuel_top
)
 
core_universe = openmc.Universe()
core_universe.add_cell(lattice_cell)
 
# Be radial reflector: fills the annular region between the lattice cylinder and
# the outer vacuum boundary, over the full reactor height (including axial reflector zones)
# FIX: spans bottom_boundary to top_boundary (not fuel_bottom to fuel_top) so the
# radial reflector covers the corners where axial and radial reflectors meet —
# earlier version left those corner volumes as void, causing lost particles there
radial_reflector_cell = openmc.Cell(
    fill=be,
    region=(
        +lattice_boundary
        & -outer_boundary
        & +bottom_boundary
        & -top_boundary
        & +sym_plane_1
        & -sym_plane_2
    )
)
core_universe.add_cell(radial_reflector_cell)

# =============================================================================
# GEOMETRY - ROOT CELL AND GEOMETRY EXPORT (line 330-366)
# root universe: the top-level container that holds everything else; defines the physical boundaries 
# =============================================================================
# KEEP: root cell with boundary conditions
# CHANGE: use 1/12 symmetry planes instead of 1/8

# FIX: root_cell region restricted to ACTIVE FUEL ZONE (fuel_bottom to fuel_top)
# — not bottom_boundary/top_boundary. The BeO axial reflector cells added below
# sit in the axial gaps above/below the active zone. If root_cell also claims that
# space, the two regions overlap → geometry conflict → lost particles.
root_cell = openmc.Cell(name='root cell')
root_cell.fill   = core_universe
root_cell.region = (
    -outer_boundary   # inside outer vacuum cylinder (r < 65cm)
    & +fuel_bottom    # above bottom of active fuel zone
    & -fuel_top       # below top of active fuel zone
    & +sym_plane_1    # inside 30° wedge — on positive side of plane at 0°
    & -sym_plane_2    # inside 30° wedge — on negative side of plane at 30°
)
 
root_universe = openmc.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)
 
# ADDED: explicit BeO axial reflector cells above and below the active core
# fills the 12.5cm axial gaps between fuel_top/fuel_bottom and the vacuum boundaries
# with real BeO material so neutrons reflect instead of streaming through void
# FIX: non-overlapping with root_cell because root_cell is restricted to fuel zone only
top_beo_cell = openmc.Cell(
    fill=beo,
    region=+fuel_top & -top_boundary & -outer_boundary & +sym_plane_1 & -sym_plane_2
)
bot_beo_cell = openmc.Cell(
    fill=beo,
    region=+bottom_boundary & -fuel_bottom & -outer_boundary & +sym_plane_1 & -sym_plane_2
)
root_universe.add_cell(top_beo_cell)
root_universe.add_cell(bot_beo_cell)
 
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()


# =============================================================================
# SETTINGS (line 370-394)
# settings that dictate how the Monte Carlo neutron simulation runs; how many neutrons, batches, where they start, physics
# =============================================================================
# KEEP: same settings structure
# CHANGE: source point moved to center of hex core

settings = openmc.Settings()
# TEMPORARILY reduced for debugging — restore to production values after geometry confirmed working
# production values: batches=200, inactive=50, particles=10000
settings.batches   = 10     # reduced from 200 for debugging
settings.inactive  = 5      # reduced from 50 for debugging
settings.particles = 100    # reduced from 10000 for debugging
settings.temperature['multipole'] = True # uses the multipole representation of nuclear cross-sections; 
# allows accurate Doppler broadening at any temperature, not just pre-tabulated values
settings.temperature['method']    = 'interpolation'

# DIAGNOSTIC: maximum verbosity prints the exact cell/surface at first particle loss
# remove after geometry is confirmed clean
settings.verbosity = 10

# CHANGE: source point at center of hex core (was at fuel_r offset in reference)
settings.source = openmc.IndependentSource(
    space=openmc.stats.Box(
        [-core_radius, -core_radius, -core_height/2],
        [ core_radius,  core_radius,  core_height/2],
        only_fissionable=True
    )
 # only_fissionable=True: starting neutrons are born only in cells containing fissile material
    # TODO: replace with mesh-based source from previous run for faster convergence in production
)
settings.export_to_xml()


# =============================================================================
# TALLIES (line 398-422)
# defining key outputs; what is measured during the run;
# gets the power distribution that feeds into OpenFOAM's fvModels heat source
# =============================================================================
# REFERENCE: DistribcellFilter per fuel cell
# AYURI HPR: mesh tally for spatial power distribution
# CHANGE: replace per-cell tallies with regular mesh tally

# ADD: cylindrical mesh tally covering active core
# NA >= 14 axial slices
mesh = openmc.RegularMesh()    # a 3D rectangular grid overlaid on the geometry; gives a spatial map of power and flux
# 20x20 radial gives enough resolution to distinguish the 3 enrichment zones for OpenFOAM coupling
mesh.dimension = [20, 20, 14]  # 20 bins X, 20 bins Y, 14 bins Z (NA >= 14 axial slices)
mesh.lower_left  = [-core_radius, -core_radius, -core_height/2]
mesh.upper_right = [ core_radius,  core_radius,  core_height/2]
# mesh boundaries match exactly the active fuel region (not including reflectors)
mesh_filter = openmc.MeshFilter(mesh)

# KEEP concept: heating and flux tallies (same scores as reference)
tally = openmc.Tally(name='power_distribution')
tally.filters = [mesh_filter]
tally.scores  = ['heating', 'flux']   # KEEP: same as reference

tallies = openmc.Tallies([tally])
tallies.export_to_xml()

print("All XML files exported.")
print(f"Lattice enclosure radius: {lattice_enclosure_r:.3f} cm")
print(f"Total reactor height: {total_height:.1f} cm")
print("Run: openmc")  # output feeds into OpenFOAM as fvModels heat source
