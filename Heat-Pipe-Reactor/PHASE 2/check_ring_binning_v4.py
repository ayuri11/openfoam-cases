# Usage (in the PHASE 2 folder):  python3 check_ring_binning_v4.py        # compares half-widths 45 and 60 cm
#                                 python3 check_ring_binning_v4.py 60     # one half-width
import numpy as np, math, sys, json
q = np.load('power_density_3d.npy'); nx, ny, nz = q.shape
q2 = q.sum(axis=2); P = 5.5
def circ(r):
    if r == 0: return [(0.0, 0.0)]
    return [(P*r*math.cos(math.radians(i*60/r)), P*r*math.sin(math.radians(i*60/r))) for i in range(6*r)]
def hexring(r, rot=0.0):
    if r == 0: return [(0.0, 0.0)]
    cs = [(r*P*math.cos(math.radians(60*k+rot)), r*P*math.sin(math.radians(60*k+rot))) for k in range(6)]
    return [(cs[k][0]+(cs[(k+1)%6][0]-cs[k][0])*i/r, cs[k][1]+(cs[(k+1)%6][1]-cs[k][1])*i/r) for k in range(6) for i in range(r)]
Rs = [float(a) for a in sys.argv[1:]] or [45.0, 60.0]
out = {}
for R in Rs:
    x = np.linspace(-R+R/nx, R-R/nx, nx); xx, yy = np.meshgrid(x, x, indexing='ij')
    print(f"\n=== mesh half-width {R:g} cm ({2*R/nx:.2f} cm voxels) ===")
    for name, fn in (("B hex, neighbours at 0,60 (OpenMC orientation 'x')", lambda r: hexring(r, 0.0)), ("C hex, neighbours at 30,90", lambda r: hexring(r, 30.0))):
        cen, ring = [], []
        for r in range(4):
            for c in fn(r): cen.append(c); ring.append(r)
        cen = np.array(cen); ring = np.array(ring)
        owner = ring[((xx[...,None]-cen[:,0])**2 + (yy[...,None]-cen[:,1])**2).argmin(axis=2)]
        w = [2.5e6*float(q2[owner == r].sum()/q2.sum()) for r in range(4)]
        print(f"{name:55s} ring W: {[f'{v:,.0f}' for v in w]}")
        if name.startswith('B'): out[str(R)] = w
json.dump(out, open('ring_binning_v4.json', 'w'), indent=1)
print("\nUse the layout-B line for the half-width that matches the run (the XML says 60).")
