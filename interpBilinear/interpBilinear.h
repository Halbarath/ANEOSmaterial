/*  This file is part of interpBilinear.
 *  Copyright (c) 2020-2026 Thomas Meier & Christian Reinhardt
 *
 *  interpBilinear is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  interpBilinear is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with interpBilinear.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * interpBilinear: bilinear interpolation on a rectangular grid, and its
 * inversion along one axis.
 *
 * Header-only (every function is static inline) so that several libraries
 * linked into one binary can each carry their own copy, via git subtree,
 * without duplicate symbols. The consumers are ANEOSmaterial, scvheos and
 * reos3; EOSlib never calls it directly.
 *
 * The core is agnostic about what the numbers mean. A library that
 * interpolates in (log10 rho, log10 T) passes logarithms in and
 * exponentiates the results; a library that interpolates linear values
 * passes them as they are. Likewise it never prints, never exits and never
 * returns a sentinel value: every function that can fail returns a status
 * code and leaves its outputs untouched on failure. The calling library
 * maps the codes to its own contract (-1e50, exit, fallback to the first
 * row, counting in nLookupFail, ...).
 *
 * Conventions
 *   Grid: nx nodes on the outer axis x, ny on the inner axis y, both
 *   strictly increasing with at least two nodes. A table is one contiguous
 *   block, x-major: z[ix*ny + iy] = z(x[ix], y[iy]), so the inner index iy
 *   is contiguous in memory. (In the EOS libraries x is rho and y is T.)
 *
 *   Cell (ix, iy) is the rectangle [x[ix], x[ix+1]] x [y[iy], y[iy+1]],
 *   0 <= ix <= nx-2, 0 <= iy <= ny-2, with local coordinates (u, v) in
 *   [0, 1]^2. Both intervals are CLOSED: x == x[nx-1] gives the last cell
 *   with u == 1. A roundoff tolerance accepts and clamps coordinates within
 *   a band of
 *
 *     tol * |bound| + tolAbs
 *
 *   outside the first and the last node, so that a state on the table edge
 *   survives a unit conversion or a 10^x <-> log10 round trip. The two parts
 *   model the two kinds of axis:
 *     - a linear axis (ANEOSmaterial: rho, T) carries a relative roundoff,
 *       tol relative to the bound and tolAbs = 0; an absolute part would be
 *       wrong there (a density axis may start at 1e-25);
 *     - a log axis (scvheos, reos3: log10 rho, log10 T) carries the roundoff
 *       of log10(10^x), which is absolute (~1e-16) plus relative to |x|, so
 *       those libraries pass tolAbs > 0 as well. With tol alone a log bound
 *       at or near 0, which an arbitrary code-unit shift of log rho can
 *       produce, got no tolerance at all (up to 1.0.3).
 *   The functions taking one tol use tolAbs = 0; the *Tol variants
 *   (interpLocateAxisTol, interpLocateTol, interpInvertYTol, interpInvertXTol)
 *   take both. tol = tolAbs = 0 disables the tolerance.
 *
 *   The interpolant is C0; on a cell boundary the cell to the right/above
 *   is used (the last cell at the upper edge), so derivatives there are
 *   those of that cell.
 *
 * Version history
 *   1.0.0  extracted from ANEOSmaterial 1.0.4 (interpBilinear.c) and the
 *          scvheos 1.1 / reos3 1.1 cell arithmetic; closed intervals with
 *          roundoff tolerance; status codes instead of sentinels; one
 *          inversion with a strict (monotonic) and a tolerant mode.
 *   1.0.1  interpLocateAxis: edge tolerance relative to the bound instead
 *          of tol*(1 + |x|), which was an absolute 1e-13 on a linear axis
 *          starting at 1e-25 and let far-out-of-range densities through.
 *   1.0.2  INTERP_INVERT_SCAN renamed INTERP_INVERT_TOLERANT and the flag
 *          bScan renamed bFallback: the mode bisects like the strict one and
 *          only scans as a fallback, which the old name suggested otherwise.
 *          No change in behaviour.
 *   1.0.3  INTERP_ABI_VERSION, the layout version of the public types, so that the
 *          libraries vendoring this header can detect a copy with another layout in
 *          the same translation unit (see below). No change in behaviour.
 *          Edge tolerance band tol*|bound| + tolAbs: new interpLocateAxisTol,
 *          interpLocateTol, interpInvertYTol and interpInvertXTol with an absolute
 *          part for log axes (a log bound at 0 got no tolerance). The functions
 *          with one tol are unchanged (tolAbs = 0); the types and
 *          INTERP_ABI_VERSION are unchanged. New API: when the version is bumped
 *          (1.1.0), scvheos and reos3, which call the *Tol variants, must require it.
 *          interpInvertLine: the check for a non-bracketing segment after the strict
 *          bisection is gone, and with it INTERP_NONMONOTONIC (-7): it could not fire. The
 *          strict bisection keeps f_lo <= target <= f_hi on any line, monotonic or not, so it
 *          always ends on a bracketing segment and returns a root; on a non-monotonic line
 *          that is the root the bisection lands on (the comments, and the scvheos and reos3
 *          documentation built on them, said such a line fails). Strict mode: results
 *          bit-identical. Tolerant mode: a target equal to the first node value returns that
 *          node at once (k = 0, t = 0). On a monotonic line that is what it gave before, after
 *          the sign bisection had moved away to n-2 and the O(n) scan found the node; on a
 *          non-monotonic line the bisection could also end on another root, which is now the
 *          first node (the lowest-index convention of the tolerant mode). Removing
 *          INTERP_NONMONOTONIC is an API change.
 */
#ifndef INTERPBILINEAR_HINCLUDED
#define INTERPBILINEAR_HINCLUDED

#include <math.h>
#include <stddef.h>

#define INTERP_VERSION_TEXT   "1.0.3"
#define INTERP_VERSION_MAJOR  1
#define INTERP_VERSION_MINOR  0
#define INTERP_VERSION_PATCH  3

/*
 * Layout version of the public types INTERPGRID, INTERPCELL and INTERPINVERTOPT. Bump it whenever
 * one of them changes (members, order or types), independently of the version above.
 *
 * Several libraries vendor this header (ANEOSmaterial, scvheos, reos3, each as a git subtree) and
 * embed INTERPGRID in their material structs. EOSlib includes all their headers in one translation
 * unit, where the include guard lets only the first copy in: every library header therefore
 * checks that the copy in effect has the minimum version it needs and exactly the layout it was
 * written for (its <LIB>_INTERP_ABI_VERSION), and stops the compilation otherwise. A copy with a
 * different layout would make EOSlib and the library disagree on the struct silently.
 */
#define INTERP_ABI_VERSION 1

/*
 * Status codes. Everything that is not INTERP_SUCCESS is negative.
 */
#define INTERP_SUCCESS        0
#define INTERP_OUTSIDE_XMIN  -1   /* interpLocate*: x below the axis */
#define INTERP_OUTSIDE_XMAX  -2
#define INTERP_OUTSIDE_YMIN  -3
#define INTERP_OUTSIDE_YMAX  -4
#define INTERP_BELOW         -5   /* interpLocateAxis: below the first node;
                                     interpInvertLine (strict): target below the first node value */
#define INTERP_ABOVE         -6   /* ... above the last one */
                                  /* -7 unused (was INTERP_NONMONOTONIC, never returned) */
#define INTERP_NOT_ATTAINED  -8   /* interpInvertLine (tolerant): no segment attains the target */
#define INTERP_INVALID       -9   /* a NaN input */

#define INTERP_ALONG_X 0   /* vary ix at fixed (iy, v) */
#define INTERP_ALONG_Y 1   /* vary iy at fixed (ix, u) */

typedef struct interpGrid {
    int nx;            /* nodes on the outer axis */
    int ny;            /* nodes on the inner (contiguous) axis */
    const double *x;   /* nx strictly increasing nodes */
    const double *y;   /* ny strictly increasing nodes */
} INTERPGRID;

typedef struct interpCell {
    int ix;            /* 0 <= ix <= nx-2 */
    int iy;            /* 0 <= iy <= ny-2 */
    double u;          /* local coordinate along x, in [0, 1] */
    double v;          /* local coordinate along y, in [0, 1] */
} INTERPCELL;

/*
 * Options of interpInvertLine.
 *
 * Both modes search by bisection on the node index, O(log n). A bisection
 * whose end nodes bracket the target keeps a bracketing segment at every
 * step, on any line, monotonic or not, so it always ends on a root (on a
 * non-monotonic line on the one it bisects onto, which need not be the
 * lowest). The modes differ in what happens when the end nodes do not
 * bracket the target, which on a non-monotonic line does not mean that the
 * target is not attained (a dip below f_0, a rise above f_n-1):
 *
 *   bFallback 0: strict mode, for lines that increase in the normal case.
 *                The end nodes must bracket the target, f_0 <= target <=
 *                f_n-1: a target below fails with INTERP_BELOW, above with
 *                INTERP_ABOVE, also where a non-monotonic line attains it
 *                inside. A target between them always gives a root.
 *             1: tolerant mode, for lines of either direction and for
 *                defective tables. If the end nodes bracket the target, a
 *                sign bisection finds a root, and that is the whole search
 *                in the normal case; a target equal to f_0 is that node.
 *                Only if the end nodes do not bracket the target (which
 *                needs a non-monotonic line, or a target off the line) does
 *                it fall back to a linear scan for the lowest-index
 *                bracketing segment, and if none brackets the target
 *                exactly, to the edge and node tolerances below, before
 *                INTERP_NOT_ATTAINED.
 *   valTol    relative roundoff tolerance on the target at the line ends: a
 *             target within valTol*(1 + |target|) outside [f_0, f_n-1] is
 *             clamped onto the nearer end instead of failing (both modes).
 *   edgeTol   tolerant mode only: a target overshooting the first (last)
 *             node by at most edgeTol times the span of the first (last)
 *             segment is clamped onto that node.
 *   nodeTol   tolerant mode only: a target within nodeTol*max(|target|, |f_k|)
 *             of some node value f_k (roundoff on a plateau) is located at
 *             that node, lowest index first.
 *
 * The two option sets the EOS libraries use:
 */
typedef struct interpInvertOpt {
    int bFallback;
    double valTol;
    double edgeTol;
    double nodeTol;
} INTERPINVERTOPT;

#define INTERP_INVERT_STRICT    { 0, 1e-13, 0.0,  0.0   }   /* scvheos, reos3 */
#define INTERP_INVERT_TOLERANT  { 1, 0.0,   1e-4, 1e-12 }   /* ANEOSmaterial */

/* ---------------------------------------------------------------- access */

static inline double interpZ(const INTERPGRID *g, const double *z, int ix, int iy)
{
    return z[(size_t)ix * g->ny + iy];
}

static inline int interpAxisLength(const INTERPGRID *g, int axis)
{
    return (axis == INTERP_ALONG_Y) ? g->ny : g->nx;
}

static inline const double *interpAxisNodes(const INTERPGRID *g, int axis)
{
    return (axis == INTERP_ALONG_Y) ? g->y : g->x;
}

/* ---------------------------------------------------------------- locate */

/*
 * Index i with axis[i] <= x < axis[i+1] for a strictly increasing axis of
 * n >= 2 nodes, clamped to [0, n-2]: x < axis[0] gives 0 and x >= axis[n-1]
 * gives n-2, so the returned cell is always valid (and the interpolant then
 * extrapolates linearly). No tolerance, no failure; callers that must reject
 * out-of-range input use interpLocateAxis. O(log n).
 */
static inline int interpFindIndex(double x, const double *axis, int n)
{
    int lo = 0;
    int hi = n - 1;

    if (x < axis[0]) return 0;
    if (x >= axis[n - 1]) return n - 2;

    /* Invariant: axis[lo] <= x < axis[hi]. */
    while (hi - lo > 1) {
        int mid = lo + (hi - lo) / 2;
        if (axis[mid] <= x) lo = mid;
        else                hi = mid;
    }
    return lo;
}

/*
 * Locate x on the closed interval [axis[0], axis[n-1]] with the roundoff
 * tolerance band tol*|bound| + tolAbs at each bound (see the header comment).
 * Returns INTERP_SUCCESS with the segment index *pi (0 <= i <= n-2) and the
 * local coordinate *pw in [0, 1]; INTERP_BELOW or INTERP_ABOVE if x is outside
 * by more than the tolerance; INTERP_INVALID for a NaN. Outputs are untouched
 * on failure.
 */
static inline int interpLocateAxisTol(const double *axis, int n, double x, double tol, double tolAbs,
                                      int *pi, double *pw)
{
    const double lo = axis[0], hi = axis[n - 1];
    int i;

    if (x != x) return INTERP_INVALID;
    if (!(x >= lo)) { if (x >= lo - (tol * fabs(lo) + tolAbs)) x = lo; else return INTERP_BELOW; }
    if (!(x <= hi)) { if (x <= hi + (tol * fabs(hi) + tolAbs)) x = hi; else return INTERP_ABOVE; }

    i = interpFindIndex(x, axis, n);
    *pi = i;
    *pw = (x - axis[i]) / (axis[i + 1] - axis[i]);
    return INTERP_SUCCESS;
}

/* Same with a tolerance relative to the bound only (linear axes). */
static inline int interpLocateAxis(const double *axis, int n, double x, double tol, int *pi, double *pw)
{
    return interpLocateAxisTol(axis, n, x, tol, 0.0, pi, pw);
}

/*
 * Locate (x, y) on the grid with the edge tolerance band tol*|bound| + tolAbs.
 * Returns INTERP_SUCCESS and fills *cell, or one of INTERP_OUTSIDE_XMIN/XMAX/
 * YMIN/YMAX (x is checked first), or INTERP_INVALID for a NaN.
 */
static inline int interpLocateTol(const INTERPGRID *g, double x, double y, double tol, double tolAbs, INTERPCELL *cell)
{
    int iret;
    int ix, iy;
    double u, v;

    iret = interpLocateAxisTol(g->x, g->nx, x, tol, tolAbs, &ix, &u);
    if (iret == INTERP_BELOW) return INTERP_OUTSIDE_XMIN;
    if (iret == INTERP_ABOVE) return INTERP_OUTSIDE_XMAX;
    if (iret != INTERP_SUCCESS) return iret;

    iret = interpLocateAxisTol(g->y, g->ny, y, tol, tolAbs, &iy, &v);
    if (iret == INTERP_BELOW) return INTERP_OUTSIDE_YMIN;
    if (iret == INTERP_ABOVE) return INTERP_OUTSIDE_YMAX;
    if (iret != INTERP_SUCCESS) return iret;

    cell->ix = ix;
    cell->iy = iy;
    cell->u = u;
    cell->v = v;
    return INTERP_SUCCESS;
}

/* Same with a tolerance relative to the bounds only (linear axes). */
static inline int interpLocate(const INTERPGRID *g, double x, double y, double tol, INTERPCELL *cell)
{
    return interpLocateTol(g, x, y, tol, 0.0, cell);
}

/* -------------------------------------------------------------- evaluate */

/*
 * Bilinear value of a table in a located cell:
 *
 *   z(u,v) = (1-u)(1-v) f00 + u(1-v) f10 + (1-u) v f01 + u v f11
 *
 * No bounds checking: the caller owns the cell.
 */
static inline double interpCellValue(const INTERPGRID *g, const double *z, const INTERPCELL *c)
{
    const double *p0 = z + (size_t)c->ix * g->ny + c->iy;
    const double *p1 = p0 + g->ny;
    return (1.0 - c->u) * ((1.0 - c->v) * p0[0] + c->v * p0[1])
         +        c->u  * ((1.0 - c->v) * p1[0] + c->v * p1[1]);
}

/*
 * Partial derivatives of the bilinear interpolant in a located cell, with
 * respect to the grid coordinates x and y:
 *
 *   dz/dx = [ (1-v)(f10 - f00) + v (f11 - f01) ] / (x[ix+1] - x[ix])
 *   dz/dy = [ (1-u)(f01 - f00) + u (f11 - f10) ] / (y[iy+1] - y[iy])
 *
 * Exact for the interpolant; each is linear in the other coordinate and
 * jumps at the cell boundaries.
 */
static inline void interpCellGradient(const INTERPGRID *g, const double *z, const INTERPCELL *c,
                                      double *dzdx, double *dzdy)
{
    const double *p0 = z + (size_t)c->ix * g->ny + c->iy;
    const double *p1 = p0 + g->ny;
    const double f00 = p0[0], f01 = p0[1], f10 = p1[0], f11 = p1[1];

    *dzdx = ((1.0 - c->v) * (f10 - f00) + c->v * (f11 - f01)) / (g->x[c->ix + 1] - g->x[c->ix]);
    *dzdy = ((1.0 - c->u) * (f01 - f00) + c->u * (f11 - f10)) / (g->y[c->iy + 1] - g->y[c->iy]);
}

/*
 * z(x, y): locate and evaluate. Returns the status of interpLocate.
 */
static inline int interpValue(const INTERPGRID *g, const double *z, double x, double y, double tol, double *pz)
{
    INTERPCELL c;
    int iret = interpLocate(g, x, y, tol, &c);
    if (iret != INTERP_SUCCESS) return iret;
    *pz = interpCellValue(g, z, &c);
    return INTERP_SUCCESS;
}

/*
 * dz/dx and dz/dy at (x, y): locate and differentiate.
 */
static inline int interpGradient(const INTERPGRID *g, const double *z, double x, double y, double tol,
                                 double *dzdx, double *dzdy)
{
    INTERPCELL c;
    int iret = interpLocate(g, x, y, tol, &c);
    if (iret != INTERP_SUCCESS) return iret;
    interpCellGradient(g, z, &c, dzdx, dzdy);
    return INTERP_SUCCESS;
}

/*
 * Nearest-node lookup in an integer table of the same layout (phase tables).
 * Returns the status of interpLocate and the value at the node nearest to
 * (x, y) in *pz. Exactly halfway, the lower node wins.
 */
static inline int interpNearest(const INTERPGRID *g, const int *z, double x, double y, double tol, int *pz)
{
    INTERPCELL c;
    int ix, iy;
    int iret = interpLocate(g, x, y, tol, &c);
    if (iret != INTERP_SUCCESS) return iret;
    ix = c.ix + (c.u > 0.5);
    iy = c.iy + (c.v > 0.5);
    *pz = z[(size_t)ix * g->ny + iy];
    return INTERP_SUCCESS;
}

/* --------------------------------------------------------------- lines */

/*
 * Fixing one coordinate reduces the bilinear interpolant to a piecewise
 * linear function of the other, with node values
 *
 *   along y (fixed ix, u):  f_k = (1-u) z[ix][k]   + u z[ix+1][k],  k = 0..ny-1
 *   along x (fixed iy, v):  f_k = (1-v) z[k][iy]   + v z[k][iy+1],  k = 0..nx-1
 *
 * interpLineValue returns f_k. `fixed` is the segment index of the fixed
 * axis (0 <= fixed <= n_fixed-2) and w its local coordinate.
 */
static inline double interpLineValue(const INTERPGRID *g, const double *z, int axis, int fixed, double w, int k)
{
    if (axis == INTERP_ALONG_Y)
        return (1.0 - w) * interpZ(g, z, fixed, k) + w * interpZ(g, z, fixed + 1, k);
    else
        return (1.0 - w) * interpZ(g, z, k, fixed) + w * interpZ(g, z, k, fixed + 1);
}

/*
 * Coordinate on the varying axis of the point at segment k, position t.
 */
static inline double interpLineCoord(const INTERPGRID *g, int axis, int k, double t)
{
    const double *a = interpAxisNodes(g, axis);
    return a[k] + t * (a[k + 1] - a[k]);
}

/*
 * Node line k of the TABLE (not of the interpolant): the raw values
 * z[fixed][k] (along y) or z[k][fixed] (along x), 0 <= fixed <= n_fixed-1.
 * Returns 0 if the line is strictly increasing in k, otherwise the lowest
 * index k >= 1 with f_k <= f_k-1. Bilinear interpolation between two
 * monotonic node lines is monotonic, so checking every node line of the
 * table covers every interpolated line. Used for the monotonicity report
 * at initialisation.
 */
static inline int interpNodeLineMonotonic(const INTERPGRID *g, const double *z, int axis, int fixed)
{
    const int n = interpAxisLength(g, axis);
    const size_t stride = (axis == INTERP_ALONG_Y) ? 1 : (size_t)g->ny;
    const double *p = (axis == INTERP_ALONG_Y) ? z + (size_t)fixed * g->ny : z + fixed;
    int k;

    for (k = 1; k < n; k++) {
        if (!(p[(size_t)k * stride] > p[(size_t)(k - 1) * stride])) return k;
    }
    return 0;
}

/* -------------------------------------------------------------- invert */

/*
 * Linear scan for the lowest-index segment [k, k+1] whose node values
 * bracket the target (inclusive). Returns 1 and fills k, flo, fhi, or 0.
 */
static inline int interpScanLine(const INTERPGRID *g, const double *z, int axis, int fixed, double w,
                                 int n, double target, int *pk, double *pflo, double *pfhi)
{
    double fhi = interpLineValue(g, z, axis, fixed, w, 0);
    int k;

    for (k = 0; k < n - 1; k++) {
        double flo = fhi;
        fhi = interpLineValue(g, z, axis, fixed, w, k + 1);
        if ((flo <= target && target <= fhi) || (fhi <= target && target <= flo)) {
            *pk = k;
            *pflo = flo;
            *pfhi = fhi;
            return 1;
        }
    }
    return 0;
}

/*
 * Relaxed fallback of the tolerant mode once no segment brackets the target:
 * the edge tolerance (clamp onto the first/last node) and the node
 * tolerance (roundoff on a plateau), see INTERPINVERTOPT. Returns 1 and
 * fills k, t on success, 0 otherwise.
 */
static inline int interpRelaxedLocate(const INTERPGRID *g, const double *z, int axis, int fixed, double w,
                                      int n, double target, const INTERPINVERTOPT *opt, int *pk, double *pt)
{
    const double g0 = interpLineValue(g, z, axis, fixed, w, 0);
    const double g1 = interpLineValue(g, z, axis, fixed, w, n - 1);
    const double gb = interpLineValue(g, z, axis, fixed, w, 1);
    const double gt = interpLineValue(g, z, axis, fixed, w, n - 2);
    const double zbelow = g0 - opt->edgeTol * (gb - g0);
    const double zabove = g1 + opt->edgeTol * (g1 - gt);
    int k;

    if ((zbelow <= target && target <= g0) || (g0 <= target && target <= zbelow)) { *pk = 0;     *pt = 0.0; return 1; }
    if ((g1 <= target && target <= zabove) || (zabove <= target && target <= g1)) { *pk = n - 2; *pt = 1.0; return 1; }
    if (opt->nodeTol > 0.0) {
        for (k = 0; k < n - 1; k++) {
            double gk = interpLineValue(g, z, axis, fixed, w, k);
            if (fabs(target - gk) <= opt->nodeTol * fmax(fabs(target), fabs(gk))) { *pk = k; *pt = 0.0; return 1; }
        }
    }
    return 0;
}

/*
 * Invert the piecewise linear line: find the segment k and t in [0, 1] with
 *
 *   target = (1-t) f_k + t f_k+1.
 *
 * See INTERPINVERTOPT for the two modes and their failure codes. On a flat
 * segment (f_k == f_k+1 == target) t = 0, the lower node. When several roots
 * exist the result depends on the path: the bisection returns the root it
 * bisects onto (normally the one in the monotonic bulk of the line), the
 * scan fallback of the tolerant mode the lowest-index one. Neither mode
 * detects a non-monotonic line (that would take O(n) per call); the libraries
 * report such lines of their tables once, at initialisation.
 *
 * Returns INTERP_SUCCESS with *pk, *pt; INTERP_BELOW/INTERP_ABOVE (strict) or
 * INTERP_NOT_ATTAINED (tolerant) if the target is not found; INTERP_INVALID
 * for a NaN target. The bracket invariant needs finite line values, which
 * the libraries guarantee by validating their tables at load.
 */
static inline int interpInvertLine(const INTERPGRID *g, const double *z, int axis, int fixed, double w,
                                   double target, const INTERPINVERTOPT *opt, int *pk, double *pt)
{
    const int n = interpAxisLength(g, axis);
    double g0 = interpLineValue(g, z, axis, fixed, w, 0);
    double g1 = interpLineValue(g, z, axis, fixed, w, n - 1);
    int lo = 0, hi = n - 1;
    double flo = g0, fhi = g1;
    double t;

    if (target != target) return INTERP_INVALID;

    /*
     * Roundoff tolerance at the line ends: a target that left the line by a
     * 10^x <-> log10 round trip is clamped back onto the nearer end.
     */
    if (opt->valTol > 0.0 && fabs(target) < HUGE_VAL) {
        const double eps = opt->valTol * (1.0 + fabs(target));
        const double lo_end = fmin(g0, g1), hi_end = fmax(g0, g1);
        if (target < lo_end && target >= lo_end - eps) target = lo_end;
        if (target > hi_end && target <= hi_end + eps) target = hi_end;
    }

    if (!opt->bFallback) {
        /*
         * Strict mode: the end nodes must bracket the target. The bisection keeps
         * flo <= target < fhi (or target == fhi == g1), so it ends on a root.
         */
        if (!(target >= g0)) return INTERP_BELOW;
        if (!(target <= g1)) return INTERP_ABOVE;
        while (hi - lo > 1) {
            int mid = lo + (hi - lo) / 2;
            double fm = interpLineValue(g, z, axis, fixed, w, mid);
            if (fm <= target) { lo = mid; flo = fm; }
            else              { hi = mid; fhi = fm; }
        }
    } else {
        /*
         * The first node itself: the lowest-index root. The sign bisection below
         * needs flo != target and would move away from it (to n-2, then the scan).
         */
        if (g0 == target) { *pk = 0; *pt = 0.0; return INTERP_SUCCESS; }
        if ((g0 < target) != (g1 < target) || g1 == target) {
            /*
             * End nodes bracket the target: sign bisection, valid on any line. flo
             * stays strictly on the side of g0, fhi on the other side or at the
             * target, so the segment [lo, hi] brackets the target at the end.
             */
            while (hi - lo > 1) {
                int mid = lo + (hi - lo) / 2;
                double fm = interpLineValue(g, z, axis, fixed, w, mid);
                if ((fm < target) == (flo < target) && fm != target) { lo = mid; flo = fm; }
                else                                                 { hi = mid; fhi = fm; }
            }
        } else {
            /* Not bracketed by the end nodes: lowest-index bracketing segment, then the tolerances. */
            if (!interpScanLine(g, z, axis, fixed, w, n, target, &lo, &flo, &fhi)) {
                if (interpRelaxedLocate(g, z, axis, fixed, w, n, target, opt, pk, pt)) return INTERP_SUCCESS;
                return INTERP_NOT_ATTAINED;
            }
            hi = lo + 1;
        }
    }

    /* Solve the linear segment [lo, hi]. */
    if (fhi == flo) t = 0.0;
    else            t = (target - flo) / (fhi - flo);
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    *pk = lo;
    *pt = t;
    return INTERP_SUCCESS;
}

/*
 * Solve z(x, y) = target for y at given x. Locates x (tolerance band
 * tol*|bound| + tolAbs), then
 * inverts the line along y. On success fills *cell (so that other tables
 * can be evaluated in the same cell without a second search) and *py.
 *
 * Returns INTERP_OUTSIDE_XMIN/XMAX, INTERP_INVALID, or the status of
 * interpInvertLine.
 */
static inline int interpInvertYTol(const INTERPGRID *g, const double *z, double x, double target, double tol,
                                   double tolAbs, const INTERPINVERTOPT *opt, INTERPCELL *cell, double *py)
{
    int ix, iy, iret;
    double u, v;

    iret = interpLocateAxisTol(g->x, g->nx, x, tol, tolAbs, &ix, &u);
    if (iret == INTERP_BELOW) return INTERP_OUTSIDE_XMIN;
    if (iret == INTERP_ABOVE) return INTERP_OUTSIDE_XMAX;
    if (iret != INTERP_SUCCESS) return iret;

    iret = interpInvertLine(g, z, INTERP_ALONG_Y, ix, u, target, opt, &iy, &v);
    if (iret != INTERP_SUCCESS) return iret;

    cell->ix = ix; cell->iy = iy; cell->u = u; cell->v = v;
    *py = interpLineCoord(g, INTERP_ALONG_Y, iy, v);
    return INTERP_SUCCESS;
}

/* Same with a tolerance relative to the bounds only (linear axes). */
static inline int interpInvertY(const INTERPGRID *g, const double *z, double x, double target, double tol,
                                const INTERPINVERTOPT *opt, INTERPCELL *cell, double *py)
{
    return interpInvertYTol(g, z, x, target, tol, 0.0, opt, cell, py);
}

/*
 * Solve z(x, y) = target for x at given y. Mirror image of interpInvertY;
 * returns INTERP_OUTSIDE_YMIN/YMAX if y is outside the grid.
 */
static inline int interpInvertXTol(const INTERPGRID *g, const double *z, double y, double target, double tol,
                                   double tolAbs, const INTERPINVERTOPT *opt, INTERPCELL *cell, double *px)
{
    int ix, iy, iret;
    double u, v;

    iret = interpLocateAxisTol(g->y, g->ny, y, tol, tolAbs, &iy, &v);
    if (iret == INTERP_BELOW) return INTERP_OUTSIDE_YMIN;
    if (iret == INTERP_ABOVE) return INTERP_OUTSIDE_YMAX;
    if (iret != INTERP_SUCCESS) return iret;

    iret = interpInvertLine(g, z, INTERP_ALONG_X, iy, v, target, opt, &ix, &u);
    if (iret != INTERP_SUCCESS) return iret;

    cell->ix = ix; cell->iy = iy; cell->u = u; cell->v = v;
    *px = interpLineCoord(g, INTERP_ALONG_X, ix, u);
    return INTERP_SUCCESS;
}

/* Same with a tolerance relative to the bounds only (linear axes). */
static inline int interpInvertX(const INTERPGRID *g, const double *z, double y, double target, double tol,
                                const INTERPINVERTOPT *opt, INTERPCELL *cell, double *px)
{
    return interpInvertXTol(g, z, y, target, tol, 0.0, opt, cell, px);
}

/* ------------------------------------------------------- node derivatives */

/*
 * Derivative of f at node k of a strictly increasing axis x[0..n-1] by the
 * second-order central difference on a (possibly non-uniform) grid,
 *
 *   f'(x_k) = [h0^2 f_k+1 - h1^2 f_k-1 + (h1^2 - h0^2) f_k] / (h0 h1 (h0 + h1)),
 *   h0 = x_k - x_k-1,  h1 = x_k+1 - x_k,
 *
 * exact for quadratics, and the one-sided difference at the two ends. f is
 * accessed with a stride (1 along y, ny along x). Used to build the
 * sound-speed tables of scvheos and reos3 from P and u.
 */
static inline double interpNodeDeriv(const double *x, int n, const double *f, size_t stride, int k)
{
    if (k == 0)     return (f[stride] - f[0]) / (x[1] - x[0]);
    if (k == n - 1) return (f[(size_t)k * stride] - f[(size_t)(k - 1) * stride]) / (x[k] - x[k - 1]);
    {
        const double h0 = x[k] - x[k - 1], h1 = x[k + 1] - x[k];
        const double fm = f[(size_t)(k - 1) * stride], f0 = f[(size_t)k * stride], fp = f[(size_t)(k + 1) * stride];
        return (h0 * h0 * fp - h1 * h1 * fm + (h1 * h1 - h0 * h0) * f0) / (h0 * h1 * (h0 + h1));
    }
}

#endif