import numpy as np 
import openmc 
import math

# =============================================================================
# DEFINING MATERIALS SELECTED
# =============================================================================
# REFERENCE CODE USED: U-10Mo fuel (HEU 93%), single enrichment
# OUR HPR USES: UO2 fuel, three HALEU enrichment zones

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

# Fuel zones normalized perfectly to 1.0 stoichiometry
fuel_zone1 = openmc.Material(name='UO2_12pct')
fuel_zone1.set_density('g/cm3', 10.4)
fuel_zone1.add_nuclide('U235', 0.12 * (1.0 / 3.0), 'ao') 
fuel_zone1.add_nuclide('U238', 0.88 * (1.0 / 3.0), 'ao') 
fuel_zone1.add_nuclide('O16',  2.0 / 3.0,          'ao') 

fuel_zone2 = openmc.Material(name='UO2_15pct')
fuel_zone2.set_density('g/cm3', 10.4)
fuel_zone2.add_nuclide('U235', 0.15 * (1.0 / 3.0), 'ao') 
fuel_zone2.add_nuclide('U238', 0.85 * (1.0 / 3.0), 'ao')
fuel_zone2.add_nuclide('O16',  2.0 / 3.0,          'ao') 

fuel_zone3 = openmc.Material(name='UO2_1975pct')
fuel_zone3.set_density('g/cm3', 10.4)
fuel_zone3.add_nuclide('U235', 0.1975 * (1.0 / 3.0), 'ao')
fuel_zone3.add_nuclide('U238', 0.8025 * (1.0 / 3.0), 'ao')
fuel_zone3.add_nuclide('O16',  2.0 / 3.0,            'ao')

materials = openmc.Materials([
    haynes, b4c, beo, be, sodium,
    graphite, fuel_zone1, fuel_zone2, fuel_zone3
])
materials.export_to_xml()

# =============================================================================
# GEOMETRY - PARAMETERS
# =============================================================================
core_height    = 160.0   
hp_radius      = 0.450   
hp_wall_thick  = 0.089   
fuel_pin_r     = 0.350   
ctrl_rod_r     = 0.450   
cell_flat      = 5.5     

axial_ref_top    = 12.5  
axial_ref_bottom = 12.5  
total_height = core_height + axial_ref_top + axial_ref_bottom 

core_radius      = 45.0  
reflector_radius = 60.0  

# =============================================================================
# GEOMETRY - SURFACES
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
# GEOMETRY - UNIT CELL UNIVERSE BUILDER
# =============================================================================
def build_unit_cell(fuel_material, center_type, graphite_material, sodium_mat, haynes_mat, b4c_mat):
    pin_ring1_r = 0.95  
    pin_ring2_r = 1.85  
    hp_ring_r   = 2.65  

    cells = []
    d = cell_flat / 2.0
    s32 = math.sqrt(3.0) / 2.0

    px1 = openmc.XPlane(x0=-d)
    px2 = openmc.XPlane(x0=d)
    p1  = openmc.Plane(a=0.5,  b=s32,  c=0.0, d=d)
    p2  = openmc.Plane(a=0.5,  b=-s32, c=0.0, d=d)
    p3  = openmc.Plane(a=-0.5, b=s32,  c=0.0, d=d)
    p4  = openmc.Plane(a=-0.5, b=-s32, c=0.0, d=d)

    hex_region = (+px1 & -px2 & -p1 & -p2 & +p3 & +p4)
    graphite_region = hex_region

    for i in range(6):
        ang = math.radians(i * 60)
        x0_pos = pin_ring1_r * math.cos(ang)
        y0_pos = pin_ring1_r * math.sin(ang)
        pin_s = openmc.ZCylinder(x0=x0_pos, y0=y0_pos, r=fuel_pin_r)
        cells.append(openmc.Cell(fill=fuel_material, region=-pin_s & hex_region))
        graphite_region &= +pin_s

    for i in range(6):
        ang = math.radians(i * 60 + 30)
        x0_pos = pin_ring2_r * math.cos(ang)
        y0_pos = pin_ring2_r * math.sin(ang)
        pin_s = openmc.ZCylinder(x0=x0_pos, y0=y0_pos, r=fuel_pin_r)
        cells.append(openmc.Cell(fill=fuel_material, region=-pin_s & hex_region))
        graphite_region &= +pin_s

    for i in range(6):
        ang = math.radians(i * 60)
        x0_pos = hp_ring_r * math.cos(ang)
        y0_pos = hp_ring_r * math.sin(ang)
        hp_inner_s = openmc.ZCylinder(x0=x0_pos, y0=y0_pos, r=hp_radius - hp_wall_thick)
        hp_outer_s = openmc.ZCylinder(x0=x0_pos, y0=y0_pos, r=hp_radius)
        cells.append(openmc.Cell(fill=sodium_mat, region=-hp_inner_s & hex_region))
        cells.append(openmc.Cell(fill=haynes_mat, region=(+hp_inner_s & -hp_outer_s) & hex_region))
        graphite_region &= +hp_outer_s

    cr_s = openmc.ZCylinder(x0=0, y0=0, r=ctrl_rod_r)
    if center_type == 'control_rod':
        cells.append(openmc.Cell(fill=b4c_mat, region=-cr_s & hex_region))
    elif center_type == 'heat_pipe':
        hp_inner_c = openmc.ZCylinder(x0=0, y0=0, r=hp_radius - hp_wall_thick)
        cells.append(openmc.Cell(fill=sodium_mat, region=-hp_inner_c & hex_region))
        cells.append(openmc.Cell(fill=haynes_mat, region=(+hp_inner_c & -cr_s) & hex_region))
    graphite_region &= +cr_s

    cells.append(openmc.Cell(fill=graphite_material, region=graphite_region))
    return openmc.Universe(cells=cells)

# =============================================================================
# GEOMETRY - HEX LATTICE
# =============================================================================
zone1_univ = build_unit_cell(fuel_zone1, 'control_rod', graphite, sodium, haynes, b4c)
zone2_univ = build_unit_cell(fuel_zone2, 'control_rod', graphite, sodium, haynes, b4c)
zone3_univ = build_unit_cell(fuel_zone3, 'heat_pipe',   graphite, sodium, haynes, b4c)

ref_d = cell_flat / 2.0
ref_s32 = math.sqrt(3.0) / 2.0
r_px1 = openmc.XPlane(x0=-ref_d)
r_px2 = openmc.XPlane(x0=ref_d)
r_p1  = openmc.Plane(a=0.5,  b=ref_s32,  c=0.0, d=ref_d)
r_p2  = openmc.Plane(a=0.5,  b=-ref_s32, c=0.0, d=ref_d)
r_p3  = openmc.Plane(a=-0.5, b=ref_s32,  c=0.0, d=ref_d)
r_p4  = openmc.Plane(a=-0.5, b=-ref_s32, c=0.0, d=ref_d)

reflector_hex_region = (+r_px1 & -r_px2 & -r_p1 & -r_p2 & +r_p3 & +r_p4)
reflector_block_univ = openmc.Universe()
reflector_block_univ.add_cell(openmc.Cell(fill=be, region=reflector_hex_region))

lattice = openmc.HexLattice()
lattice.center = (0.0, 0.0)
lattice.pitch = (cell_flat,)    
lattice.orientation = 'x'              

outer_universe = openmc.Universe()
outer_cell = openmc.Cell(fill=be)
outer_universe.add_cell(outer_cell)
lattice.outer = outer_universe

lattice.universes = [
    [reflector_block_univ] * 36, 
    [reflector_block_univ] * 30, 
    [reflector_block_univ] * 24, 
    [zone3_univ] * 18,           
    [zone2_univ] * 12,           
    [zone1_univ] * 6,            
    [zone1_univ],                
]

D_outer = (cell_flat * 6.5) * (math.sqrt(3.0) / 2.0)
s32 = math.sqrt(3.0) / 2.0
c_px1 = openmc.XPlane(x0=-D_outer)
c_px2 = openmc.XPlane(x0=D_outer)
c_p1  = openmc.Plane(a=0.5,  b=s32,  c=0.0, d=D_outer)
c_p2  = openmc.Plane(a=0.5,  b=-s32, c=0.0, d=D_outer)
c_p3  = openmc.Plane(a=-0.5, b=s32,  c=0.0, d=D_outer)
c_p4  = openmc.Plane(a=-0.5, b=-s32, c=0.0, d=D_outer)

core_hex_region = (+c_px1 & -c_px2 & -c_p1 & -c_p2 & +c_p3 & +c_p4)

lattice_cell = openmc.Cell(fill=lattice, region=core_hex_region)
core_universe = openmc.Universe()
core_universe.add_cell(lattice_cell)

radial_reflector_cell = openmc.Cell(
    fill=be,
    region=(~core_hex_region & -outer_boundary & +fuel_bottom & -fuel_top & +sym_plane_1 & -sym_plane_2)
)
core_universe.add_cell(radial_reflector_cell)

# =============================================================================
# GEOMETRY - ROOT CELL
# =============================================================================
root_cell = openmc.Cell(name='root cell')
root_cell.fill   = core_universe
root_cell.region = (-outer_boundary & +fuel_bottom & -fuel_top & +sym_plane_1 & -sym_plane_2)

root_universe = openmc.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)

top_beo_cell = openmc.Cell(fill=beo, region=+fuel_top & -top_boundary & -outer_boundary & +sym_plane_1 & -sym_plane_2)
bot_beo_cell = openmc.Cell(fill=beo, region=+bottom_boundary & -fuel_bottom & -outer_boundary & +sym_plane_1 & -sym_plane_2)
root_universe.add_cell(top_beo_cell)
root_universe.add_cell(bot_beo_cell)

geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

# =============================================================================
# SETTINGS (Reverted back to unconstrained source fallback configuration)
# =============================================================================
settings = openmc.Settings()
settings.batches   = 150    
settings.inactive  = 30     
settings.particles = 5000   
settings.run_mode  = 'eigenvalue' 
settings.temperature['multipole'] = True 
settings.temperature['method']    = 'interpolation'

# REVERTED: Working, broad spatial distribution box that doesn't trigger boundary rejection
# Added a dummy initial guess for k to anchor the generation loop
settings.keff_guess = 1.0
settings.source = openmc.IndependentSource(
    space=openmc.stats.Box(
        [-20.0, -20.0, -core_height/2],
        [ 20.0,  20.0,  core_height/2],
        only_fissionable=False
    )
)
settings.export_to_xml()

# =============================================================================
# TALLIES
# =============================================================================
mesh = openmc.RegularMesh()    
mesh.dimension = [20, 20, 14]  
mesh.lower_left  = [-reflector_radius, -reflector_radius, -core_height/2]
mesh.upper_right = [ reflector_radius,  reflector_radius,  core_height/2]
mesh_filter = openmc.MeshFilter(mesh)

tally = openmc.Tally(name='power_distribution')
tally.filters = [mesh_filter]
tally.scores  = ['heating', 'flux']   

tallies = openmc.Tallies([tally])
tallies.export_to_xml()

print("All XML configurations updated. Running K-Eigenvalue track...")
openmc.run()
