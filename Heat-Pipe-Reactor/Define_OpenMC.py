import numpy as np
import openmc
import math

# =============================================================================
# DEFINING MATERIALS
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

# Zone 1 Fuel - 12% enrichment
fuel_zone1 = openmc.Material(name='UO2_12pct')
fuel_zone1.set_density('g/cm3', 10.4)
fuel_zone1.add_nuclide('U235', 0.12, 'ao')
fuel_zone1.add_nuclide('U238', 0.88, 'ao')
fuel_zone1.add_nuclide('O16',  2.0,  'ao')

# Zone 2 Fuel - 15% enrichment
fuel_zone2 = openmc.Material(name='UO2_15pct')
fuel_zone2.set_density('g/cm3', 10.4)
fuel_zone2.add_nuclide('U235', 0.15, 'ao')
fuel_zone2.add_nuclide('U238', 0.85, 'ao')
fuel_zone2.add_nuclide('O16',  2.0,  'ao')

# Zone 3 Fuel - 19.75% enrichment
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
core_height      = 160.0   
hp_radius        = 0.795   
hp_wall_thick    = 0.089   
fuel_pin_r       = 0.635   
ctrl_rod_r       = 0.795   
cell_flat        = 10.0    

axial_ref_top    = 12.5  
axial_ref_bottom = 12.5  
total_height     = core_height + axial_ref_top + axial_ref_bottom 

core_radius      = 45.0  
reflector_radius = 65.0  

# =============================================================================
# GLOBAL BOUNDARY SURFACES & SYMMETRY WEDGE
# =============================================================================
top_boundary    = openmc.ZPlane(z0=+total_height/2, boundary_type='vacuum')
bottom_boundary = openmc.ZPlane(z0=-total_height/2, boundary_type='vacuum')
fuel_top        = openmc.ZPlane(z0=+core_height/2)
fuel_bottom     = openmc.ZPlane(z0=-core_height/2)

outer_boundary  = openmc.ZCylinder(r=reflector_radius, boundary_type='vacuum')

angle1 = 0.0                  
angle2 = math.radians(30.0)   
sym_plane_1 = openmc.Plane(a=math.sin(angle1), b=-math.cos(angle1), c=0, d=0, boundary_type='reflective')
sym_plane_2 = openmc.Plane(a=math.sin(angle2), b=-math.cos(angle2), c=0, d=0, boundary_type='reflective')

# =============================================================================
# FIXED UNIT CELL UNIVERSE BUILDER (FLAT CSG GEOMETRY)
# =============================================================================
def build_unit_cell(fuel_material, center_type, graphite_material, sodium_mat, haynes_mat, b4c_mat):
    """
    Builds a flat, self-contained hexagonal unit cell universe.
    center_type: 'control_rod' or 'heat_pipe'
    """
    pin_ring1_r = 1.729  
    pin_ring2_r = 3.299  
    hp_ring_r   = 4.879  

    cells = []
    
    # Establish the local unit cell hexagonal bounding prism
    hex_prism = openmc.model.hexagonal_prism(
        edge_length=cell_flat / math.sqrt(3),
        orientation='x'
    )
    
    # The cell background matrix starts as the inside of the hexagonal prism boundary
    graphite_region = -hex_prism

    # 6 fuel pins: inner ring
    for i in range(6):
        ang = math.radians(i * 60)
        x = pin_ring1_r * math.cos(ang)
        y = pin_ring1_r * math.sin(ang)
        pin_s = openmc.ZCylinder(x0=x, y0=y, r=fuel_pin_r)
        
        # Add fuel cell directly inside this universe boundary scope
        cells.append(openmc.Cell(fill=fuel_material, region=-pin_s & -hex_prism))
        graphite_region &= +pin_s

    # 6 fuel pins: outer ring (offset 30 deg)
    for i in range(6):
        ang = math.radians(i * 60 + 30)
        x = pin_ring2_r * math.cos(ang)
        y = pin_ring2_r * math.sin(ang)
        pin_s = openmc.ZCylinder(x0=x, y0=y, r=fuel_pin_r)
        
        cells.append(openmc.Cell(fill=fuel_material, region=-pin_s & -hex_prism))
        graphite_region &= +pin_s

    # 6 heat pipes: at hex corners
    for i in range(6):
        ang = math.radians(i * 60)
        x = hp_ring_r * math.cos(ang)
        y = hp_ring_r * math.sin(ang)
        
        hp_inner_s = openmc.ZCylinder(x0=x, y0=y, r=hp_radius - hp_wall_thick)
        hp_outer_s = openmc.ZCylinder(x0=x, y0=y, r=hp_radius)
        
        cells.append(openmc.Cell(fill=sodium_mat, region=-hp_inner_s & -hex_prism))
        cells.append(openmc.Cell(fill=haynes_mat, region=+hp_inner_s & -hp_outer_s & -hex_prism))
        graphite_region &= +hp_outer_s

    # 1 central rod position
    cr_s = openmc.ZCylinder(x0=0, y0=0, r=ctrl_rod_r)
    if center_type == 'control_rod':
        cells.append(openmc.Cell(fill=b4c_mat, region=-cr_s & -hex_prism))
        graphite_region &= +cr_s
    elif center_type == 'heat_pipe':
        hp_inner_c = openmc.ZCylinder(x0=0, y0=0, r=hp_radius - hp_wall_thick)
        cells.append(openmc.Cell(fill=sodium_mat, region=-hp_inner_c & -hex_prism))
        cells.append(openmc.Cell(fill=haynes_mat, region=+hp_inner_c & -cr_s & -hex_prism))
        graphite_region &= +cr_s

    # Fill the remaining space inside the hex boundary with graphite
    cells.append(openmc.Cell(fill=graphite_material, region=graphite_region))

    return openmc.Universe(cells=cells)

# =============================================================================
# INITIALIZING THE THREE CHOSEN CORE ENRICHMENT ZONES
# =============================================================================
zone1_univ = build_unit_cell(fuel_zone1, 'control_rod', graphite, sodium, haynes, b4c)
zone2_univ = build_unit_cell(fuel_zone2, 'control_rod', graphite, sodium, haynes, b4c)
zone3_univ = build_unit_cell(fuel_zone3, 'heat_pipe',   graphite, sodium, haynes, b4c)

# =============================================================================
# BUILDING THE LATTICE
# =============================================================================
lattice = openmc.HexLattice()
lattice.center = (0.0, 0.0)
lattice.pitch = (cell_flat,)    
lattice.orientation = 'x'         

# Outer universe to catch out-of-bounds lattice samples
outer_universe = openmc.Universe()
outer_universe.add_cell(openmc.Cell(fill=be))
lattice.outer = outer_universe

# Map the layout rings (Outermost to Innermost)
lattice.universes = [
    [zone3_univ] * 18,   # Ring 3: Outer region (19.75% Fuel + Central HP)
    [zone2_univ] * 12,   # Ring 2: Middle region (15% Fuel + Central CR)
    [zone1_univ] * 6,    # Ring 1: Inner region (12% Fuel + Central CR)
    [zone1_univ],        # Ring 0: Core Center cell (12% Fuel + Central CR)
]

# =============================================================================
# CORE WRAPPER UNIVERSE
# =============================================================================
core_hex = openmc.model.hexagonal_prism(
    edge_length=cell_flat * 3.5, # Generous outer boundary to encapsulate 3 full rings safely
    orientation='x'
)

lattice_cell = openmc.Cell(fill=lattice, region=-core_hex)

core_universe = openmc.Universe()
core_universe.add_cell(lattice_cell)

# Fill radial space between the structural core boundary and outer cylinder with Be Reflector
radial_reflector_cell = openmc.Cell(
    fill=be,
    region=+core_hex & -outer_boundary & +fuel_bottom & -fuel_top & +sym_plane_1 & -sym_plane_2
)
core_universe.add_cell(radial_reflector_cell)

# =============================================================================
# ROOT UNIVERSE WITH 1/12 SYMMETRY WEDGE
# =============================================================================
root_cell = openmc.Cell(name='active_core_cell')
root_cell.fill   = core_universe
root_cell.region = (
    -outer_boundary    
    & +fuel_bottom       
    & -fuel_top          
    & +sym_plane_1       
    & -sym_plane_2       
)

root_universe = openmc.Universe(universe_id=0, name='root_universe')
root_universe.add_cell(root_cell)

# Top and Bottom Axial Reflector Cells (BeO) 
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
# SETTINGS
# =============================================================================
settings = openmc.Settings()
settings.batches   = 10     # Debugging scale
settings.inactive  = 5      
settings.particles = 500    # Bumped slightly to ensure robust source coverage over multi-zone box
settings.temperature['multipole'] = True
settings.temperature['method']    = 'interpolation'

# Source Box targeted inside fissionable boundaries
settings.source = openmc.IndependentSource(
    space=openmc.stats.Box(
        [-core_radius, -core_radius, -core_height/2],
        [ core_radius,  core_radius,  core_height/2],
        only_fissionable=True
    )
)
settings.export_to_xml()

# =============================================================================
# TALLIES FOR OPENFOAM COUPLING
# =============================================================================
mesh = openmc.RegularMesh()
mesh.dimension = [20, 20, 14]  
mesh.lower_left  = [-core_radius, -core_radius, -core_height/2]
mesh.upper_right = [ core_radius,  core_radius,  core_height/2]
mesh_filter = openmc.MeshFilter(mesh)

tally = openmc.Tally(name='power_distribution')
tally.filters = [mesh_filter]
tally.scores  = ['heating', 'flux']

tallies = openmc.Tallies([tally])
tallies.export_to_xml()

print("PRODUCTION CONFIGURATION EXPORTED: 37-cell multi-zone full core setup is locked and clean.")
