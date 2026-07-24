// =============================================================================
// hpr_core.geo — gmsh geometry script for HPR core thermal model
// PHASE 3B: OpenFOAM IVTBC thermal-hydraulic analysis
// Geometry: 1/12 symmetry wedge of hexagonal graphite monolith
// Includes: fuel pin zones, heat pipe holes, reflector boundary
// Units: meters (OpenFOAM standard)
// =============================================================================

// =============================================================================
// DESIGN PARAMETERS (from methodology Section 3.2)
// all dimensions in meters
// =============================================================================

// core dimensions
core_height   = 1.60;    // m — active fuel height (160 cm)
cell_flat     = 0.055;   // m — unit cell flat-to-flat (5.5 cm)

// heat pipe geometry (Price et al. 2023)
hp_r          = 0.00450; // m — heat pipe outer radius (scaled: 4.5 mm)
hp_wall       = 0.00089; // m — Haynes 230 wall thickness (0.89 mm)

// fuel pin geometry
fp_r          = 0.00350; // m — fuel pin radius (scaled: 3.5 mm)

// control rod
cr_r          = 0.00450; // m — control rod radius = hp_r

// pin ring positions within unit cell (from build_unit_cell)
pin_ring1_r   = 0.0095;  // m — inner ring of 6 fuel pins (0.95 cm)
pin_ring2_r   = 0.0185;  // m — outer ring of 6 fuel pins (1.85 cm)
hp_ring_r     = 0.0265;  // m — 6 heat pipes at hex corners (2.65 cm)

// mesh refinement levels
lc_coarse     = 0.004;   // m — background graphite mesh size
lc_fine       = 0.001;   // m — mesh size near pins and HPs
lc_hp         = 0.0008;  // m — mesh size at heat pipe surfaces (critical for BC)

// axial layers
n_axial       = 28;      // number of axial layers (NA=28 per methodology)

// =============================================================================
// GEOMETRY: 1/12 SYMMETRY WEDGE
// A regular hexagon has 12-fold symmetry (6 rotational × 2 mirror)
// We model a 30° wedge from 0° to 30°
// The wedge contains ONE representative unit cell
// =============================================================================

// hexagon apothem (flat-to-center distance) = cell_flat/2
apothem = cell_flat / 2.0;   // 0.0275 m

// hex vertex radius = cell_flat / sqrt(3)
vertex_r = cell_flat / Sqrt(3.0);   // 0.03175 m

// 1/12 wedge: 30 degree slice
// Point 1: origin (core center)
// Points 2,3: on the two symmetry planes at radius = outer reflector boundary
outer_r = 0.45;   // m — outer reflector radius (45 cm)

// =============================================================================
// POINTS — WEDGE BOUNDARY
// =============================================================================

// origin
Point(1) = {0, 0, 0, lc_coarse};

// symmetry plane 1: at 0° (along +x axis)
Point(2) = {outer_r, 0, 0, lc_coarse};

// symmetry plane 2: at 30° 
Point(3) = {outer_r * Cos(Pi/6), outer_r * Sin(Pi/6), 0, lc_coarse};

// =============================================================================
// POINTS — FUEL PINS (inner ring, 1/12 wedge has 1 pin from inner ring)
// inner ring: 6 pins at 0°, 60°, 120°, 180°, 240°, 300°
// in 30° wedge: partial pin at 0° (on symmetry plane)
// use center of first inner pin at 0°: x=pin_ring1_r, y=0
// =============================================================================

// inner ring pin center at 0° (on symmetry plane 1)
fp1_x = pin_ring1_r;
fp1_y = 0.0;

// inner ring pin at 60° (partially in wedge)
fp2_x = pin_ring1_r * Cos(Pi/3);
fp2_y = pin_ring1_r * Sin(Pi/3);

// outer ring pin at 30° (on symmetry plane 2, center of wedge)
fp3_x = pin_ring2_r * Cos(Pi/6);
fp3_y = pin_ring2_r * Sin(Pi/6);

// =============================================================================
// POINTS — HEAT PIPES (at hex corners)
// hp_ring_r = 2.65 cm from cell center
// in 30° wedge: HP at 0° (on sym plane 1) and HP at 30° (on sym plane 2)
// =============================================================================

// HP at 0° (on symmetry plane 1)
hp1_x = hp_ring_r;
hp1_y = 0.0;

// HP at 30° (on symmetry plane 2 — middle of wedge)
hp2_x = hp_ring_r * Cos(Pi/6);
hp2_y = hp_ring_r * Sin(Pi/6);

// =============================================================================
// SIMPLIFIED GEOMETRY FOR PHASE 3B
// Instead of modeling all pins/HPs explicitly (complex intersection issues),
// we use a simplified approach:
// - Graphite monolith as main solid region
// - Heat pipe holes as circular cutouts with IVTBC boundary condition
// - Fuel pins treated as volumetric heat sources within graphite
// This follows Price et al. (2023) approach exactly
// =============================================================================

// --- Central region circle (represents core active zone) ---
// radius = distance to outermost HP + hp_r
core_active_r = hp_ring_r + hp_r + 0.005;   // ~0.032 m

// --- Outer boundary of graphite block ---
// use inscribed circle of hex (apothem) for simplicity in 2D cross section
// this avoids complex hex corner geometry while preserving area
graphite_r = apothem;   // 0.0275 m

// =============================================================================
// 2D CROSS-SECTION: SIMPLIFIED WEDGE
// Build as: outer arc → symmetry lines → subtract HP holes and fuel zones
// =============================================================================

// Center point
Point(10) = {0, 0, 0, lc_coarse};

// Outer boundary points (at graphite_r on symmetry planes)
Point(11) = {graphite_r, 0, 0, lc_coarse};           // on sym plane 1
Point(12) = {graphite_r*Cos(Pi/6), graphite_r*Sin(Pi/6), 0, lc_coarse}; // on sym plane 2

// Arc from point 11 to 12 (outer graphite boundary)
Point(13) = {0, 0, 0, lc_coarse};   // center of arc = origin
Circle(1) = {11, 10, 12};

// Symmetry lines
Line(2) = {10, 11};   // sym plane 1 (0°)
Line(3) = {12, 10};   // sym plane 2 (30° → back to origin)

// Wedge surface (before subtracting holes)
Line Loop(1) = {2, 1, 3};
Plane Surface(1) = {1};

// =============================================================================
// HEAT PIPE HOLES — circular cutouts in graphite
// These surfaces get the IVTBC boundary condition in OpenFOAM
// =============================================================================

// HP 1: at 0° on symmetry plane (only half circle in wedge)
Point(20) = {hp1_x, hp1_y, 0, lc_hp};           // center
Point(21) = {hp1_x + hp_r, 0, 0, lc_hp};         // right point
Point(22) = {hp1_x, hp_r, 0, lc_hp};             // top point  
Point(23) = {hp1_x - hp_r, 0, 0, lc_hp};         // left point

Circle(10) = {21, 20, 22};
Circle(11) = {22, 20, 23};
// half circle (in wedge): arc from right to top to left
Line Loop(10) = {10, 11};   // half circle
// Note: HP on symmetry plane = half pipe in wedge

// HP 2: at 30° (fully inside wedge — whole circle)
Point(30) = {hp2_x, hp2_y, 0, lc_hp};
Point(31) = {hp2_x + hp_r, hp2_y, 0, lc_hp};
Point(32) = {hp2_x, hp2_y + hp_r, 0, lc_hp};
Point(33) = {hp2_x - hp_r, hp2_y, 0, lc_hp};
Point(34) = {hp2_x, hp2_y - hp_r, 0, lc_hp};

Circle(20) = {31, 30, 32};
Circle(21) = {32, 30, 33};
Circle(22) = {33, 30, 34};
Circle(23) = {34, 30, 31};

Line Loop(20) = {20, 21, 22, 23};   // full circle
Plane Surface(20) = {20};           // HP 2 disk (to be subtracted)

// Fuel pin 1: at 0° on symmetry plane (half pin)
Point(40) = {fp1_x, 0, 0, lc_fine};
Point(41) = {fp1_x + fp_r, 0, 0, lc_fine};
Point(42) = {fp1_x, fp_r, 0, lc_fine};
Point(43) = {fp1_x - fp_r, 0, 0, lc_fine};

Circle(30) = {41, 40, 42};
Circle(31) = {42, 40, 43};
Line Loop(30) = {30, 31};   // half circle fuel pin 1

// Fuel pin 2: at 30° (fully inside wedge)
Point(50) = {fp3_x, fp3_y, 0, lc_fine};
Point(51) = {fp3_x + fp_r, fp3_y, 0, lc_fine};
Point(52) = {fp3_x, fp3_y + fp_r, 0, lc_fine};
Point(53) = {fp3_x - fp_r, fp3_y, 0, lc_fine};
Point(54) = {fp3_x, fp3_y - fp_r, 0, lc_fine};

Circle(40) = {51, 50, 52};
Circle(41) = {52, 50, 53};
Circle(42) = {53, 50, 54};
Circle(43) = {54, 50, 51};

Line Loop(40) = {40, 41, 42, 43};
Plane Surface(40) = {40};   // fuel pin 2 disk

// =============================================================================
// MAIN GRAPHITE SURFACE (wedge minus HP holes and fuel pin holes)
// =============================================================================

// Subtract HP2 and fuel pin 2 from wedge
// HP1 and fuel pin 1 are on symmetry plane — handled as boundary segments
Plane Surface(100) = {1, -20, -40};   // wedge minus HP2 minus fp2

// =============================================================================
// PHYSICAL GROUPS — tells OpenFOAM what each surface/volume is
// =============================================================================

// Volume (after extrusion — defined after Extrude command)
// Surfaces for boundary conditions:
Physical Surface("symmetry_1") = {2};           // 0° symmetry plane
Physical Surface("symmetry_2") = {3};           // 30° symmetry plane  
Physical Surface("heatPipe_surfaces") = {20};   // HP hole surfaces → IVTBC BC
Physical Surface("fuelPin_surfaces") = {40};    // fuel pin surfaces → heat source zone
Physical Surface("outer_graphite") = {1};       // outer graphite boundary

// =============================================================================
// EXTRUDE: 2D cross-section → 3D core
// Extrude in z direction for core_height with n_axial layers
// =============================================================================

// Extrude the graphite surface
ext[] = Extrude {0, 0, core_height} {
    Surface{100};
    Layers{n_axial};
    Recombine;
};

// Extrude the HP2 hole (creates cylinder surface for IVTBC BC)
ext_hp2[] = Extrude {0, 0, core_height} {
    Surface{20};
    Layers{n_axial};
    Recombine;
};

// Extrude the fuel pin 2 zone
ext_fp2[] = Extrude {0, 0, core_height} {
    Surface{40};
    Layers{n_axial};
    Recombine;
};

// =============================================================================
// PHYSICAL VOLUMES — the solid regions OpenFOAM solves in
// =============================================================================

Physical Volume("graphite_monolith") = {ext[1]};
Physical Volume("heatPipe_zone") = {ext_hp2[1]};
Physical Volume("fuelPin_zone") = {ext_fp2[1]};

// =============================================================================
// MESH SETTINGS
// =============================================================================

// use structured hex mesh for better convergence
Mesh.RecombineAll = 1;
Mesh.Algorithm = 6;       // Frontal-Delaunay for 2D
Mesh.Algorithm3D = 4;     // Frontal for 3D

// mesh size fields for refinement near HPs (critical for IVTBC accuracy)
Field[1] = Distance;
Field[1].SurfacesList = {20};   // distance from HP2 surface
Field[1].NNodesByEdge = 40;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_hp;
Field[2].SizeMax = lc_coarse;
Field[2].DistMin = 0.002;
Field[2].DistMax = 0.015;

Background Field = 2;

Mesh.CharacteristicLengthMin = 0.0005;
Mesh.CharacteristicLengthMax = 0.005;
