import openmc
import math

# =============================================================================
# STEP 1 — single unit cell, no lattice, no symmetry
# goal: confirm unit cell geometry is clean before adding lattice complexity
# if this runs without lost particles, the unit cell build_unit_cell() is correct
# =============================================================================

# --- materials (minimal set needed for one unit cell) ---
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

# --- geometry parameters ---
cell_flat   = 10.0   # cm flat-to-flat
fuel_pin_r  = 0.635  # cm fuel pin radius
hp_radius   = 0.795  # cm heat pipe outer radius
hp_wall_t   = 0.089  # cm Haynes wall thickness
ctrl_rod_r  = 0.795  # cm control rod radius
core_height = 160.0  # cm — axial extent of unit cell

pin_ring1_r = 1.729  # cm inner ring 6 pins
pin_ring2_r = 3.299  # cm outer ring 6 pins (30° offset)
hp_ring_r   = 4.879  # cm 6 HPs at hex corners

# --- pin and HP universes ---
# heat pipe universe: sodium vapor + Haynes wall
hp_inner_s = openmc.ZCylinder(r=hp_radius - hp_wall_t)
hp_outer_s = openmc.ZCylinder(r=hp_radius)
hp_universe = openmc.Universe(cells=[
    openmc.Cell(fill=sodium,  region=-hp_inner_s),
    openmc.Cell(fill=haynes,  region=+hp_inner_s & -hp_outer_s),
])

# fuel pin universe: fuel inside, graphite outside (catch-all)
fp_surf = openmc.ZCylinder(r=fuel_pin_r)
fp_universe = openmc.Universe(cells=[
    openmc.Cell(fill=fuel,     region=-fp_surf),
    openmc.Cell(fill=graphite, region=+fp_surf),
])

# control rod universe: B4C inside, graphite outside (catch-all)
cr_surf = openmc.ZCylinder(r=ctrl_rod_r)
cr_universe = openmc.Universe(cells=[
    openmc.Cell(fill=b4c,      region=-cr_surf),
    openmc.Cell(fill=graphite, region=+cr_surf),
])

# --- build the unit cell directly (no function, fully explicit for debugging) ---
cells = []
pin_surfs = []
hp_surfs  = []

# 6 fuel pins: inner ring at 0°, 60°, 120°, 180°, 240°, 300°
for i in range(6):
    ang = math.radians(i * 60)
    x = pin_ring1_r * math.cos(ang)
    y = pin_ring1_r * math.sin(ang)
    s = openmc.ZCylinder(x0=x, y0=y, r=fuel_pin_r)
    pin_surfs.append(s)
    cells.append(openmc.Cell(fill=fp_universe, region=-s))

# 6 fuel pins: outer ring at 30°, 90°, 150°, 210°, 270°, 330°
for i in range(6):
    ang = math.radians(i * 60 + 30)
    x = pin_ring2_r * math.cos(ang)
    y = pin_ring2_r * math.sin(ang)
    s = openmc.ZCylinder(x0=x, y0=y, r=fuel_pin_r)
    pin_surfs.append(s)
    cells.append(openmc.Cell(fill=fp_universe, region=-s))

# 6 heat pipes: at hex corners 0°, 60°, 120°, 180°, 240°, 300°
for i in range(6):
    ang = math.radians(i * 60)
    x = hp_ring_r * math.cos(ang)
    y = hp_ring_r * math.sin(ang)
    s = openmc.ZCylinder(x0=x, y0=y, r=hp_radius)
    hp_surfs.append(s)
    cells.append(openmc.Cell(fill=hp_universe, region=-s))

# 1 central control rod at origin
cr_s = openmc.ZCylinder(x0=0, y0=0, r=ctrl_rod_r)
cells.append(openmc.Cell(fill=cr_universe, region=-cr_s))

# graphite fills everything outside all pins/HPs/rod
graphite_region = +cr_s
for s in pin_surfs + hp_surfs:
    graphite_region = graphite_region & +s
cells.append(openmc.Cell(fill=graphite, region=graphite_region))

# --- STEP 1 KEY: wrap unit cell in a CYLINDER not a hex prism ---
# using a ZCylinder as the outer boundary avoids any hex face tracking issues
# inscribed circle radius of hex = apothem = cell_flat/2 = 5.0 cm
# circumscribed circle radius = cell_flat/sqrt(3) = 5.774 cm
# use apothem (5.0) so cylinder stays inside the hex — conservative but avoids leaks
unit_cell_r = cell_flat / 2.0 - 0.01   # 4.99 cm — just inside the hex inscribed circle

top_plane    = openmc.ZPlane(z0=+core_height/2, boundary_type='vacuum')
bottom_plane = openmc.ZPlane(z0=-core_height/2, boundary_type='vacuum')
outer_cyl    = openmc.ZCylinder(r=unit_cell_r, boundary_type='vacuum')
# boundary_type='vacuum' = outermost boundary; neutrons that cross it escape

# bound the graphite cell laterally too (inside the outer cylinder)
# this is the critical fix: graphite must be bounded or it extends infinitely
cells[-1].region = cells[-1].region & -outer_cyl

# add axial bounds to all cells
for cell in cells:
    cell.region = cell.region & +bottom_plane & -top_plane

root_universe = openmc.Universe(universe_id=0)
root_universe.add_cells(cells)
# also add a vacuum catch-all outside the cylinder (redundant with vacuum BC but safe)
outside_cell = openmc.Cell(region=+outer_cyl)   # no fill = void/vacuum
root_universe.add_cell(outside_cell)

openmc.Geometry(root_universe).export_to_xml()

# --- settings: point source at center of first inner fuel pin ---
# inner ring pin at angle=0°: center at (pin_ring1_r, 0, 0) = (1.729, 0, 0)
settings = openmc.Settings()
settings.batches   = 20
settings.inactive  = 5
settings.particles = 200
settings.verbosity = 7   # reduced verbosity — only print warnings, not every crossing
settings.source = openmc.IndependentSource(
    space=openmc.stats.Point((pin_ring1_r, 0.0, 0.0))
    # center of inner-ring fuel pin at 0° — guaranteed inside UO2 fuel
)
settings.export_to_xml()

# --- simple flux tally ---
mesh = openmc.RegularMesh()
mesh.dimension   = [10, 10, 5]
mesh.lower_left  = [-unit_cell_r, -unit_cell_r, -core_height/2]
mesh.upper_right = [ unit_cell_r,  unit_cell_r,  core_height/2]
tally = openmc.Tally(name='flux')
tally.filters = [openmc.MeshFilter(mesh)]
tally.scores  = ['flux', 'heating']
openmc.Tallies([tally]).export_to_xml()

print("STEP 1 — single unit cell geometry exported.")
print(f"  Unit cell cylinder radius: {unit_cell_r:.4f} cm")
print(f"  Core height: {core_height:.1f} cm")
print(f"  Source at: ({pin_ring1_r:.3f}, 0.0, 0.0) — center of inner fuel pin")
print(f"  Geometry: 12 fuel pins + 6 HPs + 1 control rod + graphite")
print()
print("Run:  openmc")
print("Look for: no WARNING lines → unit cell geometry is clean")
print("If clean: proceed to step 2 (7-cell mini lattice)")
