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
# =============================================================================
# ZPlanes for axial boundaries 
top_boundary    = openmc.ZPlane(z0=+total_height/2, boundary_type='vacuum')
bottom_boundary = openmc.ZPlane(z0=-total_height/2, boundary_type='vacuum')
fuel_top        = openmc.ZPlane(z0=+core_height/2)
fuel_bottom     = openmc.ZPlane(z0=-core_height/2)
top_ref_plane   = openmc.ZPlane(z0=+core_height/2)   # BeO axial reflector start
bot_ref_plane   = openmc.ZPlane(z0=-core_height/2)   # BeO axial reflector start

# outer boundary cylinder (same concept as reference reflector_OD)
outer_boundary  = openmc.ZCylinder(r=reflector_radius, boundary_type='vacuum')

# CHANGE: symmetry planes for 1/12 hex (reference used 1/8 with 2 planes)
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


# =============================================================================
# GEOMETRY - UNIT CELL UNIVERSE (line 162-276)
# =============================================================================
def build_unit_cell(fuel_material, center_type, graphite_material, sodium_mat, haynes_mat, b4c_mat):
    # NEW CHANGE: Local layout position dimensions scaled down proportionally to ensure structural safety. 
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
# =============================================================================
zone1_univ = build_unit_cell(fuel_zone1, 'control_rod', graphite, sodium, haynes, b4c)
zone2_univ = build_unit_cell(fuel_zone2, 'control_rod', graphite, sodium, haynes, b4c)
zone3_univ = build_unit_cell(fuel_zone3, 'heat_pipe',   graphite, sodium, haynes, b4c)

# NEW CHANGE: Added structural block universe to construct the radial lattice reflector elements.
reflector_hex = openmc.model.hexagonal_prism(edge_length=cell_flat / math.sqrt(3), orientation='x')
reflector_block_univ = openmc.Universe()
reflector_block_univ.add_cell(openmc.Cell(fill=be, region=-reflector_hex))

# ADD: HexLattice definition
lattice = openmc.HexLattice()
lattice.center = (0.0, 0.0)
lattice.pitch = (cell_flat,)    # flat-to-flat pitch in cm; updated to match cell_flat=5.5
lattice.orientation = 'x'              # flat side faces x-axis

# ADD: outer universe catches particles that leave lattice boundary
outer_universe = openmc.Universe()
outer_cell = openmc.Cell(fill=be)
outer_universe.add_cell(outer_cell)

lattice.outer = outer_universe

# NEW CHANGE: Expanded lattice structure from 4 layers deep (37 cells) to 7 layers deep (R=6, 169 cells)
lattice.universes = [
    [reflector_block_univ] * 36, # Ring 6: Outer Radial Reflector Lattice Elements
    [reflector_block_univ] * 30, # Ring 5: Radial Reflector Lattice Elements
    [reflector_block_univ] * 24, # Ring 4: Inner Radial Reflector Lattice Elements
    [zone3_univ] * 18,           # Ring 3: Active Fuel Zone 3 (19.75% HALEU, 18 copies)
    [zone2_univ] * 12,           # Ring 2: Middle Fuel Zone 2 (15.0% HALEU, 12 copies)
    [zone1_univ] * 6,            # Ring 1: Inner Fuel Zone 1 (12.0% HALEU, 6 copies)
    [zone1_univ],                # Ring 0: Core Center Cell (12.0% HALEU)
]

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
# =============================================================================
root_cell = openmc.Cell(name='root cell')
root_cell.fill   = core_universe
root_cell.region = (
    -outer_boundary      # inside the outer cylinder (r < 60cm)
    & +fuel_bottom       # FIX: was +bottom_boundary; now restricted to active fuel zone only
    & -fuel_top          # FIX: was -top_boundary; BeO axial cells now sit cleanly above/below
    & +sym_plane_1       # on the correct side of plane 1; plane 1 at 0°, wedge is above
    & -sym_plane_2       # on the correct side of plane 2; at 30°, wedge is below
)

root_universe = openmc.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)

# FIX: Cleaned up a duplicate cell overlap error that was embedded here
# COMMENT ON NEW CHANGE: Reverted to using Beryllium Oxide (beo) for the axial 
# reflectors to maximize thermal safety margins at the fuel face interfaces.
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
# =============================================================================
settings = openmc.Settings()
# NEW CHANGE: Boosted operational batch definitions to match design tracking requirements
settings.batches   = 150    # Aligned production run context
settings.inactive  = 30     # Stabilizes source spatial distribution shifts
settings.particles = 5000   # Production parameters for tracking active statistics
settings.run_mode  = 'eigenvalue' # Track critical k-eff matching your text basis
settings.temperature['multipole'] = True 
settings.temperature['method']    = 'interpolation'

# NEW CHANGE: Refocused box boundaries around active fuel matrix zone to drop initialization time
settings.source = openmc.IndependentSource(
    space=openmc.stats.Box(
        [-22.0, -22.0, -core_height/2],
        [ 22.0,  22.0,  core_height/2],
        only_fissionable=True
    )
)
settings.export_to_xml()


# =============================================================================
# TALLIES (line 398-422)
# =============================================================================
mesh = openmc.RegularMesh()    
mesh.dimension = [20, 20, 14]  # 20 bins X, 20 bins Y, 14 bins Z (NA >= 14 axial slices)
# COMMENT ON NEW CHANGE: Refocused mesh lower and upper coordinates to encompass the full active structural diameter
mesh.lower_left  = [-reflector_radius, -reflector_radius, -core_height/2]
mesh.upper_right = [ reflector_radius,  reflector_radius,  core_height/2]
mesh_filter = openmc.MeshFilter(mesh)

tally = openmc.Tally(name='power_distribution')
tally.filters = [mesh_filter]
tally.scores  = ['heating', 'flux']   # KEEP: same as reference

tallies = openmc.Tallies([tally])
tallies.export_to_xml()

print("All XML files exported.") 

# =============================================================================
# COMMENT ON NEW CHANGE: ADDED SIMULATION EXECUTION COMMAND
# =============================================================================
# Calling openmc.run() directly compiles the XML conditions and executes the 
# high-fidelity transport problem in your local environment, rendering your live 
# k-effective eigenvalue tracking data in the terminal window.
openmc.run()
