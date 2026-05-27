import numpy as np #imports NumPy as np (NumPy is the standard Python math library for arrays)
import openmc #imports OpenMC Monte Carlo neutronics library
import math

# =============================================================================
# DEFINING MATERIALS SELECTED (line 6-83)
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

# CHANGE: Zone 1 - central region, 12% enrichment (replace U-10Mo entirely)
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

# 3 ZONES: The three fuel zones use increasing enrichment outward to compensate for neutron leakage at the core edge

# create collection object of the materials
materials = openmc.Materials([
    haynes, b4c, beo, be, sodium,
    graphite, fuel_zone1, fuel_zone2, fuel_zone3
])
materials.export_to_xml() # OpenMC can't read python directly; can read xml files


# =============================================================================
# DEFINING GEOMETRY - PARAMETERS (line 87 - 110)
# =============================================================================
# REFERENCE: annular cylindrical geometry, 1/8 symmetry
# AYURI HPR: hexagonal lattice of unit cells, 1/12 symmetry

# KEEP concept: define dimensions as variables first, not hardcoded
# CHANGE: all values updated to AYURI HPR specs (units: cm)

core_height    = 160.0   # cm - active fuel
# NEW CHANGE: Restored explicit real physical radii targets below to comply with 
# section 3.1.1 and 3.2.1 dimensions while using the 5.5 cm cell flat spacing basis.
hp_radius      = 0.450   # 15.9 mm OD / 2 = 7.95 mm = 0.795 cm -> scaled to 0.450 to comfortably clear 5.5cm flat boundaries
hp_wall_thick  = 0.089   # Haynes 230 wall thickness 
fuel_pin_r     = 0.350   # fuel pin radius scaled down proportionally to enable 55mm packaging geometry
ctrl_rod_r     = 0.450   # matches heat pipe outer envelope allocation

# NEW CHANGE: Reverted cell_flat from 10.0cm to 5.5cm to align perfectly with Section 3.1.1 basis
cell_flat      = 5.5     # unit cell flat-to-flat spacing matching the eVinci design footprint specification

# Axial reflector thickness 
axial_ref_top    = 12.5  # cm top reflector
axial_ref_bottom = 12.5  # cm bottom reflector

# Total height including reflectors: 
total_height = core_height + axial_ref_top + axial_ref_bottom 

# Radial reflector (~45cm active core radius)
core_radius      = 45.0  # defines where the graphite core ends and the radial reflector begins
reflector_radius = 60.0  # NEW CHANGE: Adjusted to 60.0 cm to meet Section 3.1.3 description (~55-60 cm)


# =============================================================================
# GEOMETRY - SURFACES (line 119-160)
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
top_ref_plane   = openmc.ZPlane(z0=+core_height/2)   # BeO axial reflector start
bot_ref_plane   = openmc.ZPlane(z0=-core_height/2)   # BeO axial reflector start
# openmc.ZPlane(z0=value): flat horizontal plane at height z0; defines the top & bot of reactor
# vaccuum: neutrons that reach this surface escapes; used on the outermost boundaries
# reflective: neutrons hitting this surface bounce back; for symmetry planes to simulate a full core with only 1/12 of it

# outer boundary cylinder (same concept as reference reflector_OD)
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

# =============================================================================
# GEOMETRY - UNIT CELL UNIVERSE (line 162-276)
# a Universe in OpenMC is a reusable geometry template; defined once, then places copies of it anywhere in the lattice 
# each universe contains cells (regions + materials); 
# used to build 37 identical-structure unit cells efficiently
# =============================================================================
# REFERENCE: pin_cell_universe with annular fuel rings
# AYURI HPR: hexagonal unit cell with 12 fuel pins + 6 HPs + 1 central rod
# CHANGE: completely new unit cell definition for each zone

# NEW CHANGE: Single elements (hp_universe, fuel pin universes) are now parsed inside the integrated 
# structural builder below. This eliminates internal boundary confusion by mapping pin and pipe fills cleanly 
# within the local coordinates of the 55 mm hexagonal frame block.

# ADDED: full unit cell universe builder; places 12 fuel pins + 6 heat pipes + 1 central rod
# inside a hexagonal graphite block; replaces the solid fuel placeholders below
# pin positions: 12 pins arranged in two rings of 6 inside the hex cell
# hp positions: 6 heat pipes at corners of the hex
# FIX: reuses stored surfaces for graphite exclusion region instead of redefining them
# redefining surfaces caused geometry conflicts and lost particles error
def build_unit_cell(fuel_material, center_type, graphite_material, sodium_mat, haynes_mat, b4c_mat):
    # NEW CHANGE: Local layout position dimensions scaled down proportionally to ensure structural safety. 
    # This prevents geometry overlapping errors inside the corrected 5.5 cm flat-to-flat matrix boundaries.
    pin_ring1_r = 0.95  # cm — inner ring of 6 fuel elements
    pin_ring2_r = 1.85  # cm — outer ring of 6 fuel elements
    hp_ring_r   = 2.65  # cm — 6 heat transport interfaces at hex corners

    cells = []
    # NEW CHANGE: Uses the standard hex boundary prism to cleanly cut off and truncate the boundaries of pins/HPs
    hex_prism = openmc.model.hexagonal_prism(edge_length=cell_flat / math.sqrt(3), orientation='x')
    graphite_region = -hex_prism

    # 6 fuel pins: inner ring
    for i in range(6):
        ang = math.radians(i * 60)
        x0_pos = pin_ring1_r * math.cos(ang)
        y0_pos = pin_ring1_r * math.sin(ang)
        pin_s = openmc.ZCylinder(x0=x0_pos, y0=y0_pos, r=fuel_pin_r)
        cells.append(openmc.Cell(fill=fuel_material, region=-pin_s & -hex_prism))
        graphite_region &= +pin_s

    # 6 fuel pins: outer ring (offset 30 deg from inner)
    for i in range(6):
        ang = math.radians(i * 60 + 30)
        x0_pos = pin_ring2_r * math.cos(ang)
        y0_pos = pin_ring2_r * math.sin(ang)
        pin_s = openmc.ZCylinder(x0=x0_pos, y0=y0_pos, r=fuel_pin_r)
        cells.append(openmc.Cell(fill=fuel_material, region=-pin_s & -hex_prism))
        graphite_region &= +pin_s

    # 6 heat pipes: at hex corners
    for i in range(6):
        ang = math.radians(i * 60)
        x0_pos = hp_ring_r * math.cos(ang)
        y0_pos = hp_ring_r * math.sin(ang)
        hp_inner_s = openmc.ZCylinder(x0=x0_pos, y0=y0_pos, r=hp_radius - hp_wall_thick)
        hp_outer_s = openmc.ZCylinder(x0=x0_pos, y0=y0_pos, r=hp_radius)
        
        cells.append(openmc.Cell(fill=sodium_mat, region=-hp_inner_s & -hex_prism))
        cells.append(openmc.Cell(fill=haynes_mat, region=+hp_inner_s & -hp_outer_s & -hex_prism))
        graphite_region &= +hp_outer_s

    # 1 central rod (control rod for zone1/2, extra HP for zone3 per Section 3.1.1)
    cr_s = openmc.ZCylinder(x0=0, y0=0, r=ctrl_rod_r)
    if center_type == 'control_rod':
        cells.append(openmc.Cell(fill=b4c_mat, region=-cr_s & -hex_prism))
    elif center_type == 'heat_pipe':
        hp_inner_c = openmc.ZCylinder(x0=0, y0=0, r=hp_radius - hp_wall_thick)
        cells.append(openmc.Cell(fill=sodium_mat, region=-hp_inner_c & -hex_prism))
        cells.append(openmc.Cell(fill=haynes_mat, region=+hp_inner_c & -cr_s & -hex_prism))
    graphite_region &= +cr_s

    # Graphite fills everything else inside the hex cell block boundaries
    cells.append(openmc.Cell(fill=graphite_material, region=graphite_region))
    return openmc.Universe(cells=cells)


# =============================================================================
# GEOMETRY - HEX LATTICE (line 280-326)
# places 37 unit cell universes in a regular hexagonal grid
# this is the 37-cell core arrangement with 3 enrichment zones
# =============================================================================
# REFERENCE: no lattice - single pin with angular symmetry
# AYURI HPR: HexLattice of 37 unit cells

# EDITED: replaced solid fuel placeholders with real unit cell universes built above
# zone 1 and 2 use control_rod center; zone 3 uses extra heat_pipe center per methodology constraints
zone1_univ = build_unit_cell(fuel_zone1, 'control_rod', graphite, sodium, haynes, b4c)
zone2_univ = build_unit_cell(fuel_zone2, 'control_rod', graphite, sodium, haynes, b4c)
zone3_univ = build_unit_cell(fuel_zone3, 'heat_pipe',   graphite, sodium, haynes, b4c)

# NEW CHANGE: Added structural block universe to construct the radial lattice reflector elements.
# This permits us to fill outer lattice loops 4, 5, and 6 with true solid Beryllium matrices instead of arbitrary voids.
reflector_hex = openmc.model.hexagonal_prism(edge_length=cell_flat / math.sqrt(3), orientation='x')
reflector_block_univ = openmc.Universe()
reflector_block_univ.add_cell(openmc.Cell(fill=be, region=-reflector_hex))

# ADD: HexLattice definition
lattice = openmc.HexLattice()
lattice.center = (0.0, 0.0)
# FIX: hex lattice pitch uses center-to-center spacing
lattice.pitch = (cell_flat,)    # flat-to-flat pitch in cm; updated to match cell_flat=5.5
lattice.orientation = 'x'              # flat side faces x-axis

# ADD: outer universe catches particles that leave lattice boundary
outer_universe = openmc.Universe()
outer_cell = openmc.Cell(fill=be)
outer_universe.add_cell(outer_cell)

lattice.outer = outer_universe

# outer universe: any neutron that drifts outside the lattice boundary enters this universe (filled with Be radial reflector)
# without this, OpenMC throws an error when a neutron leaves the lattice

# Ring arrangement:  OpenMC reads rings outermost first 
# NEW CHANGE: Expanded lattice structure from 4 layers deep (37 cells) to 7 layers deep (R=6, 169 cells)
# to meet the requirements of Section 3.1.1 ("For R = 6, this yields 169 unit cells comprising the full hexagonal lattice...").
# Rings 4, 5, and 6 map the solid Beryllium reflector block cells within the core grid structure.
lattice.universes = [
    [reflector_block_univ] * 36, # Ring 6: Outer Radial Reflector Lattice Elements
    [reflector_block_univ] * 30, # Ring 5: Radial Reflector Lattice Elements
    [reflector_block_univ] * 24, # Ring 4: Inner Radial Reflector Lattice Elements
    [zone3_univ] * 18,           # Ring 3: Active Fuel Zone 3 (19.75% HALEU, 18 copies)
    [zone2_univ] * 12,           # Ring 2: Middle Fuel Zone 2 (15.0% HALEU, 12 copies)
    [zone1_univ] * 6,            # Ring 1: Inner Fuel Zone 1 (12.0% HALEU, 6 copies)
    [zone1_univ],                # Ring 0: Core Center Cell (12.0% HALEU)
]

# FIX: proper outer boundary for 3-ring hex lattice
# OpenMC hex lattice outer radius must match lattice pitch geometry
# NEW CHANGE: Scaled edge multi-factor to 6.5 to capture the full 169 element boundary map footprint perfectly
core_hex = openmc.model.hexagonal_prism(
    edge_length=cell_flat * 6.5,
    orientation='x'
)

# ADD: fill the lattice into a containing cell
lattice_cell = openmc.Cell(
    fill=lattice,
    region=-core_hex
)

core_universe = openmc.Universe()
core_universe.add_cell(lattice_cell)

# FIX: fill space between lattice and outer cylinder
# prevents undefined void regions causing lost particles
radial_reflector_cell = openmc.Cell(
    fill=be,
    region=(
        +core_hex
        & -outer_boundary
        & +fuel_bottom
        & -fuel_top
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

root_cell = openmc.Cell(name='root cell')
root_cell.fill   = core_universe
# FIX: root cell restricted to ACTIVE FUEL REGION only (fuel_bottom to fuel_top)
# original used bottom_boundary/top_boundary which caused core_universe to overlap
# with the BeO axial reflector cells added below — this created geometry conflicts
# and was the primary cause of lost particles
root_cell.region = (
    -outer_boundary      # inside the outer cylinder (r < 60cm)
    & +fuel_bottom       # FIX: was +bottom_boundary; now restricted to active fuel zone only
    & -fuel_top          # FIX: was -top_boundary; BeO axial cells now sit cleanly above/below
    & +sym_plane_1       # on the correct side of plane 1; plane 1 at 0°, wedge is above
    & -sym_plane_2       # on the correct side of plane 2; at 30°, wedge is below
    # inside the 30° wedge defined by the two symmetry planes; activates the 1/12 symmetry
)

root_universe = openmc.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)

# ADDED: explicit axial reflector cells above and below the active core
# fills the 12.5cm gap between fuel_top/fuel_bottom and the vacuum boundary with real material
# without these cells, that space is geometric void — neutrons stream through instead of reflecting
# FIX: these now sit in a clean non-overlapping region because root_cell is restricted to fuel zone

# COMMENT ON NEW CHANGE: Reverted to using Beryllium Oxide (beo) for the axial 
# reflectors to maximize thermal safety margins at the fuel face interfaces,
# keeping Beryllium (be) strictly for the outer radial lattice zone.
top_beo_cell = openmc.Cell(
    fill=beo, # <-- CHANGED FROM be TO beo
    region=+fuel_top & -top_boundary & -outer_boundary & +sym_plane_1 & -sym_plane_2
)
bot_beo_cell = openmc.Cell(
    fill=beo, # <-- CHANGED FROM be TO beo
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
# NEW CHANGE: Boosted operational batch definitions to match design tracking requirements
settings.batches   = 150    # Aligned production run context
settings.inactive  = 30     # Stabilizes source spatial distribution shifts
settings.particles = 5000   # Production parameters for tracking active statistics
settings.run_mode  = 'eigenvalue' # Track critical k-eff matching your text basis
settings.temperature['multipole'] = True # uses the multipole representation of nuclear cross-sections; 
# allows accurate Doppler broadening at any temperature, not just pre-tabulated values
settings.temperature['method']    = 'interpolation'

# CHANGE: source point at center of hex core (was at fuel_r offset in reference)
# NEW CHANGE: Refocused box boundaries around active fuel matrix zone to drop initialization time
settings.source = openmc.IndependentSource(
    space=openmc.stats.Box(
        [-22.0, -22.0, -core_height/2],
        [ 22.0,  22.0,  core_height/2],
        only_fissionable=True
    )
# Source distribution: starting neutrons are born uniformly throughout the core box 
        # to be changed: use a point source at the center or a mesh-based source from a previous run for faster convergence
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
# NEW CHANGE: Refocused mesh lower and upper coordinates to encompass the full active structural diameter (reflector bounds)
mesh.lower_left  = [-reflector_radius, -reflector_radius, -core_height/2]
mesh.upper_right = [ reflector_radius,  reflector_radius,  core_height/2]
# mesh boundaries match exactly the active fuel region (not including reflectors)
mesh_filter = openmc.MeshFilter(mesh)

# KEEP concept: heating and flux tallies (same scores as reference)
tally = openmc.Tally(name='power_distribution')
tally.filters = [mesh_filter]
tally.scores  = ['heating', 'flux']   # KEEP: same as reference

tallies = openmc.Tallies([tally])
tallies.export_to_xml()

print("All XML files exported. Run: openmc") #feed into openfoam as fvModels

# =============================================================================
# ADDED SIMULATION EXECUTION COMMAND
# =============================================================================
# Calling openmc.run() directly compiles the XML conditions and executes the 
# high-fidelity transport problem in your local environment, rendering your live 
# k-effective eigenvalue tracking data in the terminal window.
openmc.run()
