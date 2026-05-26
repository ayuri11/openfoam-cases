import openmc
import math

# --- Materials ---
graphite = openmc.Material(name='Graphite')
graphite.set_density('g/cm3', 1.7)
graphite.add_element('C', 1.0, 'ao')

fuel = openmc.Material(name='UO2_12pct')
fuel.set_density('g/cm3', 10.4)
fuel.add_nuclide('U235', 0.12, 'ao')
fuel.add_nuclide('U238', 0.88, 'ao')
fuel.add_nuclide('O16',  2.0,  'ao')

sodium = openmc.Material(name='Na')
sodium.set_density('g/cm3', 0.76)
sodium.add_element('Na', 1.0, 'ao')

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

openmc.Materials([graphite, fuel, sodium, haynes, b4c]).export_to_xml()

# --- Geometry Parameters ---
cell_flat   = 10.0   
fuel_pin_r  = 0.635  
hp_radius   = 0.795  
hp_wall_t   = 0.089  
ctrl_rod_r  = 0.795  
core_height = 160.0  
pin_ring1_r = 1.729  
pin_ring2_r = 3.299  
hp_ring_r   = 4.879  

# --- Pin and HP Universes (FIXED: No infinite background catch-alls) ---
hp_inner_s = openmc.ZCylinder(r=hp_radius - hp_wall_t)
hp_outer_s = openmc.ZCylinder(r=hp_radius)
hp_universe = openmc.Universe(cells=[
    openmc.Cell(fill=sodium,  region=-hp_inner_s),
    openmc.Cell(fill=haynes,  region=+hp_inner_s & -hp_outer_s),
])

fp_surf = openmc.ZCylinder(r=fuel_pin_r)
fp_universe = openmc.Universe(cells=[
    openmc.Cell(fill=fuel,     region=-fp_surf)
])

cr_surf = openmc.ZCylinder(r=ctrl_rod_r)
cr_universe = openmc.Universe(cells=[
    openmc.Cell(fill=b4c,      region=-cr_surf)
])

# --- Outer Boundaries ---
unit_cell_r  = (cell_flat / 2.0) - 0.01   
top_plane    = openmc.ZPlane(z0=+core_height/2, boundary_type='vacuum')
bottom_plane = openmc.ZPlane(z0=-core_height/2, boundary_type='vacuum')
outer_cyl    = openmc.ZCylinder(r=unit_cell_r, boundary_type='vacuum')

axial_region = +bottom_plane & -top_plane

# --- Build Root Universe ---
root_universe = openmc.Universe(universe_id=0)
cells = []

# Keep track of everything cutting into the graphite background
graphite_region = -outer_cyl & axial_region

# 6 fuel pins: inner ring
for i in range(6):
    ang = math.radians(i * 60)
    x = pin_ring1_r * math.cos(ang)
    y = pin_ring1_r * math.sin(ang)
    s = openmc.ZCylinder(x0=x, y0=y, r=fuel_pin_r)
    
    cells.append(openmc.Cell(fill=fp_universe, region=-s & axial_region))
    graphite_region &= +s

# 6 fuel pins: outer ring
for i in range(6):
    ang = math.radians(i * 60 + 30)
    x = pin_ring2_r * math.cos(ang)
    y = pin_ring2_r * math.sin(ang)
    s = openmc.ZCylinder(x0=x, y0=y, r=fuel_pin_r)
    
    cells.append(openmc.Cell(fill=fp_universe, region=-s & axial_region))
    graphite_region &= +s

# 6 heat pipes
for i in range(6):
    ang = math.radians(i * 60)
    x = hp_ring_r * math.cos(ang)
    y = hp_ring_r * math.sin(ang)
    s = openmc.ZCylinder(x0=x, y0=y, r=hp_radius)
    
    cells.append(openmc.Cell(fill=hp_universe, region=-s & axial_region))
    graphite_region &= +s

# 1 central control rod
cr_s = openmc.ZCylinder(x0=0, y0=0, r=ctrl_rod_r)
cells.append(openmc.Cell(fill=cr_universe, region=-cr_s & axial_region))
graphite_region &= +cr_s

# Background graphite cell explicitly bounded
graphite_cell = openmc.Cell(fill=graphite, region=graphite_region)
cells.append(graphite_cell)

root_universe.add_cells(cells)
openmc.Geometry(root_universe).export_to_xml()

# --- Settings ---
settings = openmc.Settings()
settings.batches   = 20
settings.inactive  = 5
settings.particles = 500
settings.run_mode  = 'fixed source' # Forcing fixed source since k-eff might be < 1 in an open cyl

settings.source = openmc.IndependentSource(
    space=openmc.stats.Point((pin_ring1_r, 0.0, 0.0))
)
settings.export_to_xml()

print("STEP 1 (Fixed) — Geometry exported successfully.")
