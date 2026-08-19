"""
flow.py — generates an animated SVG banner from an actual potential-flow solve.

Physics
-------
Incompressible, irrotational flow past a circular cylinder. The streamfunction
psi satisfies Laplace's equation

    d2psi/dx2 + d2psi/dy2 = 0

with psi = 0 on the cylinder surface (it is a streamline) and psi -> U_inf * y
in the far field. Discretised with a five-point stencil on a uniform Cartesian
grid and solved directly with a sparse LU factorisation — at this size that is
both faster and less fiddly than iterating.

A point vortex of strength Gamma is superposed analytically. Since
psi_vortex = -(Gamma / 2*pi) * ln(r) is constant on r = R, the superposition
still satisfies the wall boundary condition exactly. Gamma controls the
asymmetry (Magnus lift) and is what we vary from run to run.

Velocities follow from u = dpsi/dy, v = -dpsi/dx. Tracer streamlines are
integrated with RK4, which also gives the travel time along each streamline —
that is what sets the per-streamline dash speed in the SVG, so particles
visibly accelerate through the gap over the cylinder.

Output
------
An SVG whose "particles" are dashed strokes with an animated stroke-dashoffset.
No SMIL, no JavaScript, no per-frame raster data: the whole animation is a
handful of CSS keyframes over static path geometry. Loops seamlessly forever
at essentially zero cost to whoever is viewing your profile.
"""

import argparse
import math
from datetime import date

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

R = 1.0                     # cylinder radius (nondimensional)
U_INF = 1.0                 # freestream speed
X_MIN, X_MAX = -3.5, 11.5   # domain extent in x (radii)
Y_MIN, Y_MAX = -2.6, 2.6    # domain extent in y (radii)
NX, NY = 451, 157           # grid points

N_STREAMLINES = 26          # tracers seeded across the inlet
SEED_HALF_SPAN = 2.15       # leaves a clean margin for the caption
RK4_DT = 0.035              # tracer integration step
RK4_MAX_STEPS = 4000

# SVG appearance
SVG_W, SVG_H = 1200, 416
BG = "#0d1117"              # matches GitHub dark; reads fine on light too
DASH_ON, DASH_OFF = 5.0, 55.0
BASE_PERIOD = 2.2           # seconds for the median streamline to advance one dash


# ----------------------------------------------------------------------------
# Solver
# ----------------------------------------------------------------------------

def solve_streamfunction(gamma):
    """Sparse direct solve of Laplace's equation for psi, plus analytic vortex.

    Five-point stencil on the interior fluid nodes. Constrained nodes (walls,
    inlet, far field, cylinder interior) get identity rows. The outlet uses a
    zero-gradient condition, psi[:, -1] - psi[:, -2] = 0.
    """
    x = np.linspace(X_MIN, X_MAX, NX)
    y = np.linspace(Y_MIN, Y_MAX, NY)
    X, Y = np.meshgrid(x, y)
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    r = np.hypot(X, Y)
    solid = r <= R

    idx = np.arange(NY * NX).reshape(NY, NX)
    rows, cols, vals = [], [], []
    rhs = np.zeros(NY * NX)

    ix2 = 1.0 / dx**2
    iy2 = 1.0 / dy**2

    for j in range(NY):
        for i in range(NX):
            n = idx[j, i]

            # Dirichlet: cylinder interior/surface is the psi = 0 streamline.
            if solid[j, i]:
                rows.append(n); cols.append(n); vals.append(1.0)
                rhs[n] = 0.0
                continue

            # Dirichlet: inlet, top, bottom carry the undisturbed freestream.
            if i == 0 or j == 0 or j == NY - 1:
                rows.append(n); cols.append(n); vals.append(1.0)
                rhs[n] = U_INF * y[j]
                continue

            # Neumann: outlet lets the flow leave without reflecting.
            if i == NX - 1:
                rows.append(n); cols.append(n); vals.append(1.0)
                rows.append(n); cols.append(idx[j, i - 1]); vals.append(-1.0)
                rhs[n] = 0.0
                continue

            rows.append(n); cols.append(n); vals.append(-2.0 * (ix2 + iy2))
            rows.append(n); cols.append(idx[j, i + 1]); vals.append(ix2)
            rows.append(n); cols.append(idx[j, i - 1]); vals.append(ix2)
            rows.append(n); cols.append(idx[j + 1, i]); vals.append(iy2)
            rows.append(n); cols.append(idx[j - 1, i]); vals.append(iy2)
            rhs[n] = 0.0

    A = sp.csr_matrix((vals, (rows, cols)), shape=(NY * NX, NY * NX))
    psi = spla.spsolve(A, rhs).reshape(NY, NX)

    resid = float(np.max(np.abs(A @ psi.ravel() - rhs)))

    # Analytic vortex: constant on r = R, so the wall BC survives untouched.
    if abs(gamma) > 1e-12:
        rc = np.maximum(r, R)
        psi = psi - (gamma / (2.0 * math.pi)) * np.log(rc / R)

    return x, y, psi, solid, A.nnz, resid


def velocity_field(x, y, psi):
    """u = dpsi/dy, v = -dpsi/dx, by central differences."""
    dpsi_dy, dpsi_dx = np.gradient(psi, y, x)
    return dpsi_dy, -dpsi_dx


# ----------------------------------------------------------------------------
# Tracer integration
# ----------------------------------------------------------------------------

def make_sampler(x, y, u, v):
    """Bilinear interpolation of (u, v) at arbitrary physical points.

    Note the index order: the fields are shaped (NY, NX), so the row index is
    the y index and the column index is the x index.
    """
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    nx, ny = len(x), len(y)

    def sample(px, py):
        fx = (px - x[0]) / dx
        fy = (py - y[0]) / dy
        if fx < 0.0 or fy < 0.0 or fx > nx - 1.001 or fy > ny - 1.001:
            return 0.0, 0.0
        i = int(fx)          # column / x
        j = int(fy)          # row / y
        tx = fx - i
        ty = fy - j
        w00 = (1 - tx) * (1 - ty)
        w10 = tx * (1 - ty)
        w01 = (1 - tx) * ty
        w11 = tx * ty
        uu = (w00 * u[j, i] + w10 * u[j, i + 1]
              + w01 * u[j + 1, i] + w11 * u[j + 1, i + 1])
        vv = (w00 * v[j, i] + w10 * v[j, i + 1]
              + w01 * v[j + 1, i] + w11 * v[j + 1, i + 1])
        return float(uu), float(vv)

    return sample


def trace(sample, x0, y0, dt=RK4_DT):
    """RK4 streamline trace. Returns (points, cumulative_time)."""
    pts = [(x0, y0)]
    t = 0.0
    px, py = x0, y0
    for _ in range(RK4_MAX_STEPS):
        k1x, k1y = sample(px, py)
        k2x, k2y = sample(px + 0.5 * dt * k1x, py + 0.5 * dt * k1y)
        k3x, k3y = sample(px + 0.5 * dt * k2x, py + 0.5 * dt * k2y)
        k4x, k4y = sample(px + dt * k3x, py + dt * k3y)

        speed = math.hypot(k1x, k1y)
        if speed < 1e-4:
            break

        px += (dt / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x)
        py += (dt / 6.0) * (k1y + 2 * k2y + 2 * k3y + k4y)
        t += dt

        if not (X_MIN <= px <= X_MAX and Y_MIN <= py <= Y_MAX):
            break
        if math.hypot(px, py) < R * 0.985:
            break
        pts.append((px, py))

    return pts, t


# ----------------------------------------------------------------------------
# SVG emission
# ----------------------------------------------------------------------------

def to_screen(pts):
    """Map physical coords onto the SVG viewBox (y flips)."""
    sx = SVG_W / (X_MAX - X_MIN)
    sy = SVG_H / (Y_MAX - Y_MIN)
    s = min(sx, sy)
    ox = (SVG_W - s * (X_MAX - X_MIN)) / 2.0
    oy = (SVG_H - s * (Y_MAX - Y_MIN)) / 2.0
    out = []
    for px, py in pts:
        out.append((ox + (px - X_MIN) * s,
                    SVG_H - oy - (py - Y_MIN) * s))
    return out, s, ox, oy


def resample(pts, n=110):
    """Resample a polyline to n points evenly spaced by arclength.

    The RK4 trace produces thousands of vertices clustered wherever the flow
    is slow. None of that detail survives at 1200px wide, and all of it costs
    file size, so we redistribute to a fixed budget.
    """
    seg = [math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(pts, pts[1:])]
    cum = [0.0]
    for d in seg:
        cum.append(cum[-1] + d)
    total = cum[-1]
    if total <= 0:
        return pts[:2]

    out, k = [], 0
    for t in np.linspace(0.0, total, n):
        while k < len(seg) and cum[k + 1] < t:
            k += 1
        if k >= len(seg):
            out.append(pts[-1])
            break
        f = (t - cum[k]) / seg[k] if seg[k] > 0 else 0.0
        out.append((pts[k][0] + f * (pts[k + 1][0] - pts[k][0]),
                    pts[k][1] + f * (pts[k + 1][1] - pts[k][1])))
    return out


def path_d(screen_pts):
    pts = resample(screen_pts)
    head = f"M{pts[0][0]:.1f} {pts[0][1]:.1f}"
    body = "".join(f"L{px:.1f} {py:.1f}" for px, py in pts[1:])
    return head + body


def polyline_length(pts):
    return sum(math.hypot(b[0] - a[0], b[1] - a[1])
               for a, b in zip(pts, pts[1:]))



def lerp_colour(f):
    """Deep blue -> cyan ramp, keyed to streamline speed."""
    a = (56, 96, 190)
    b = (86, 211, 226)
    return "#%02x%02x%02x" % tuple(int(a[k] + f * (b[k] - a[k])) for k in range(3))


def build_svg(streamlines, cx, cy, r_px, gamma, label):
    period = DASH_ON + DASH_OFF
    rules, paths, defs = [], [], []

    # Reference speed: median of (screen length / travel time).
    speeds = [ln / t for _, ln, t in streamlines if t > 0]
    ref = sorted(speeds)[len(speeds) // 2]

    for idx, (d, length, ttime) in enumerate(streamlines):
        speed = (length / ttime) if ttime > 0 else ref
        dur = BASE_PERIOD * ref / max(speed, 1e-6)
        delay = -(idx % 7) * (dur / 7.0)   # stagger so it never looks like a comb
        rules.append(
            f".s{idx}{{animation:flow {dur:.2f}s linear infinite;"
            f"animation-delay:{delay:.2f}s}}"
        )
        defs.append(f'<path id="p{idx}" d="{d}"/>')
        # Warmer/brighter where the flow is faster, so the acceleration through
        # the gap over the body is legible even in a still frame.
        f = min(max((speed / ref - 0.75) / 0.85, 0.0), 1.0)
        col = lerp_colour(f)
        paths.append(f'<use class="tr s{idx}" href="#p{idx}" '
                     f'xlink:href="#p{idx}" stroke="{col}"/>')

    css = (
        f"@keyframes flow{{to{{stroke-dashoffset:{-period:.1f}}}}}"
        f".tr{{fill:none;stroke-width:2.4;stroke-linecap:round;"
        f"stroke-dasharray:{DASH_ON} {DASH_OFF};opacity:.85}}"
        f".gl{{fill:none;stroke:#1f6feb;stroke-width:.7;opacity:.13}}"
        + "".join(rules)
        + "@media(prefers-reduced-motion:reduce){.tr{animation:none}}"
    )

    ghost = "".join(f'<use class="gl" href="#p{i}" xlink:href="#p{i}"/>'
                    for i in range(len(streamlines)))
    tracers = "".join(paths)

    return f'''<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 {SVG_W} {SVG_H}" width="{SVG_W}" height="{SVG_H}" role="img" aria-label="Animated potential flow past a cylinder">
<title>Potential flow past a cylinder — regenerated daily</title>
<style>{css}</style>\n<defs>{"".join(defs)}</defs>
<rect width="{SVG_W}" height="{SVG_H}" rx="10" fill="{BG}"/>
<g>{ghost}</g>
<g>{tracers}</g>
<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r_px:.1f}" fill="#0d1117" stroke="#30363d" stroke-width="2"/>
<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r_px:.1f}" fill="#58a6ff" opacity=".07"/>
<text x="20" y="{SVG_H - 14}" font-family="ui-monospace,SFMono-Regular,Menlo,monospace" font-size="11.5" fill="#6e7681">{label}</text>
</svg>'''


# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gamma", type=float, default=None,
                    help="vortex strength; defaults to a value derived from today's date")
    ap.add_argument("--out", default="assets/flow.svg")
    args = ap.parse_args()

    today = date.today()
    if args.gamma is None:
        # Smooth seasonal wander through +/- 3.2 so consecutive days differ
        # visibly but nothing ever looks broken.
        gamma = 3.2 * math.sin(today.toordinal() * 0.45)
    else:
        gamma = args.gamma

    x, y, psi, solid, iters, resid = solve_streamfunction(gamma)
    print(f"solved: {iters} nonzeros, max residual {resid:.2e}")

    u, v = velocity_field(x, y, psi)
    u[solid] = 0.0
    v[solid] = 0.0
    sample = make_sampler(x, y, u, v)

    seeds = np.linspace(-SEED_HALF_SPAN, SEED_HALF_SPAN, N_STREAMLINES)
    streamlines = []
    for y0 in seeds:
        pts, ttime = trace(sample, X_MIN + 0.05, float(y0))
        if len(pts) < 12:
            continue
        screen, s, ox, oy = to_screen(pts)
        d = path_d(screen)
        streamlines.append((d, polyline_length(screen), ttime))

    _, s, ox, oy = to_screen([(0.0, 0.0)])
    cx = ox + (0.0 - X_MIN) * s
    cy = SVG_H - oy - (0.0 - Y_MIN) * s
    r_px = R * s

    label = (f"potential flow  ·  Laplace, sparse direct, {NX}x{NY}  ·  "
             f"Gamma/2\u03c0RU = {gamma / (2 * math.pi):+.2f}  ·  {today.isoformat()}")

    svg = build_svg(streamlines, cx, cy, r_px, gamma, label)
    with open(args.out, "w") as f:
        f.write(svg)
    print(f"wrote {args.out}  ({len(svg) / 1024:.1f} KB, "
          f"{len(streamlines)} streamlines)")


if __name__ == "__main__":
    main()
