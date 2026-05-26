import openmc
import math

# =============================================================================
# STEP 1 (BOUNDED FLAT GEOMETRY) — No leaks outside outer_cyl
# =============================================================================

# --- Materials ---
graphite = openmc.Material(name='Graphite')
graphite.set_density('g/cm3', 1.7)
graphite.add_element('C', 1.0, 'ao')
graphite.add_s_alpha_beta('c_Graphite')

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

# --- Global Boundaries ---
# CRITICAL FIX: Make the bounding cylinder large enough to contain the heat pipes for Step 1 testing!
# Outermost edge of heat pipe is 4.879 + 0.795 = 5.674 cm. Let's make the boundary 5.80 cm.
unit_cell_r  = 5.80 
top_plane    = openmc.ZPlane(z0=+core_height/2, boundary_type='vacuum')
bottom_plane = openmc.ZPlane(z0=-core_height/2, boundary_type='vacuum')
outer_cyl    = openmc.ZCylinder(r=unit_cell_r, boundary_type='vacuum')

# Every cell in the problem MUST be bounded by this common master region
master_region = -outer_cyl & +bottom_plane & -top_plane

# --- Build Root Universe Directly ---
root_universe = openmc.Universe(universe_id=0)
cells = []

# Initialize graphite background as the full master region
graphite_region = master_region

# 6 fuel pins: inner ring
for i in range(6):
    ang = math.radians(i * 60)
    x = pin_ring1_r * math.cos(ang)
    y = pin_ring1_r * math.sin(ang)
    s = openmc.ZCylinder(x0=x, y0=y, r=fuel_pin_r)
    
    # Clip the pin cell to the master region
    cells.append(openmc.Cell(fill=fuel, region=-s & master_region))
    graphite_region &= +s

# 6 fuel pins: outer ring
for i in range(6):
    ang = math.radians(i * 60 + 30)
    x = pin_ring2_r * math.cos(ang)
    y = pin_ring2_r * math.sin(ang)
    s = openmc.ZCylinder(x0=x, y0=y, r=fuel_pin_r)
    
    # Clip the pin cell to the master region
    cells.append(openmc.Cell(fill=fuel, region=-s & master_region))
    graphite_region &= +s

# 6 heat pipes
for i in range(6):
    ang = math.radians(i * 60)
    x = hp_ring_r * math.cos(ang)
    y = hp_ring_r * math.sin(ang)
    
    s_inner = openmc.ZCylinder(x0=x, y0=y, r=hp_radius - hp_wall_t)
    s_outer = openmc.ZCylinder(x0=x, y0=y, r=hp_radius)
    
    # Clip heat pipe layers to the master region
    cells.append(openmc.Cell(fill=sodium, region=-s_inner & master_region))
    cells.append(openmc.Cell(fill=haynes, region=+s_inner & -s_outer & master_region))
    graphite_region &= +s_outer

# 1 central control rod
cr_s = openmc.ZCylinder(x0=0, y0=0, r=ctrl_rod_r)
cells.append(openmc.Cell(fill=b4c, region=-cr_s & master_region))
graphite_region &= +cr_s

# Background graphite cell
graphite_cell = openmc.Cell(fill=graphite, region=graphite_region)
cells.append(graphite_cell)

root_universe.add_cells(cells)
openmc.Geometry(root_universe).export_to_xml()

# --- Settings ---
settings = openmc.Settings()
settings.batches   = 20
settings.inactive  = 5
settings.particles = 500
settings.run_mode  = 'fixed source'

settings.source = openmc.IndependentSource(
    space=openmc.stats.Point((0.0, 0.0, 0.0))
)
settings.export_to_xml()

print(f"STEP 1 (Fixed Boundaries) — Exported with master boundary radius: {unit_cell_r} cm")
