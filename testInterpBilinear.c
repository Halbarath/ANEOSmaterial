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
 * Test suite for interpBilinear.
 *
 * Runs on synthetic tables only, so it needs no EOS table file and tests
 * the arithmetic, not any particular material. Every check is an invariant
 * of the interpolant: a globally bilinear function must be reproduced to
 * roundoff, inversions must round-trip, and constructed defects (a dip, a
 * plateau, an edge overshoot) must produce the documented status codes.
 *
 * Usage: testInterpBilinear. Exit status is the number of failed checks.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "interpBilinear.h"

/* ---------------------------------------------------------------- harness */

static int nChecks = 0;
static int nFails = 0;
static const char *currentTest = "";

#define CHECK(cond, ...) do { \
    nChecks++; \
    if (!(cond)) { \
        nFails++; \
        fprintf(stderr, "FAIL [%s] %s:%d: ", currentTest, __FILE__, __LINE__); \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, "\n"); \
    } \
} while (0)

#define CHECK_REL(a, b, tol, ...) do { \
    double _a = (a), _b = (b); \
    double _den = fmax(fabs(_a), fabs(_b)); \
    double _err = (_den > 0.0) ? fabs(_a - _b) / _den : fabs(_a - _b); \
    nChecks++; \
    if (!(_err <= (tol))) { \
        nFails++; \
        fprintf(stderr, "FAIL [%s] %s:%d: ", currentTest, __FILE__, __LINE__); \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, " (got %.15e vs %.15e, rel err %.3e > %.1e)\n", _a, _b, _err, (double)(tol)); \
    } \
} while (0)

#define TEST(name) do { currentTest = name; printf("  %-68s", name); fflush(stdout); } while (0)
#define END_TEST() do { printf("%s\n", nFails == failsBefore ? "ok" : "FAILED"); } while (0)

/* ------------------------------------------------------------ test grids */

#define NX 37
#define NY 23

static double xAxis[NX], yAxis[NY];
static double zBilin[NX * NY];   /* globally bilinear: exact for the interpolant */
static double zMono[NX * NY];    /* increasing in both x and y */
static double zDip[NX * NY];     /* zMono with a dip along y in one x column */
static double zFlat[NX * NY];    /* zMono with a plateau along y */
static int    zPhase[NX * NY];   /* integer table for the nearest-node lookup */
static INTERPGRID G;

/* A deterministic pseudo-random number in [0, 1). */
static unsigned long rngState = 12345;
static double frand(void)
{
    rngState = rngState * 6364136223846793005UL + 1442695040888963407UL;
    return (double)(rngState >> 11) / 9007199254740992.0;
}

static double fBilin(double x, double y) { return 1.0 + 2.0 * x + 0.5 * y + 0.25 * x * y; }   /* positive terms: no cancellation */
static double fMono(double x, double y)  { return pow(x, 0.3) * pow(y, 1.2) + 0.1 * x; }

static void setupGrids(void)
{
    int i, j;
    G.nx = NX; G.ny = NY; G.x = xAxis; G.y = yAxis;

    /* Log-spaced axes with a slight non-uniform jitter so that no step is special. */
    for (i = 0; i < NX; i++) xAxis[i] = pow(10.0, -2.0 + 4.0 * (i + 0.3 * sin(1.7 * i)) / (NX - 1));
    for (j = 0; j < NY; j++) yAxis[j] = pow(10.0, -1.0 + 3.0 * (j + 0.3 * cos(2.3 * j)) / (NY - 1));
    for (i = 1; i < NX; i++) if (!(xAxis[i] > xAxis[i-1])) { fprintf(stderr, "bad x axis\n"); exit(99); }
    for (j = 1; j < NY; j++) if (!(yAxis[j] > yAxis[j-1])) { fprintf(stderr, "bad y axis\n"); exit(99); }

    for (i = 0; i < NX; i++)
        for (j = 0; j < NY; j++) {
            zBilin[i * NY + j] = fBilin(xAxis[i], yAxis[j]);
            zMono[i * NY + j]  = fMono(xAxis[i], yAxis[j]);
            zDip[i * NY + j]   = zMono[i * NY + j];
            zFlat[i * NY + j]  = zMono[i * NY + j];
            zPhase[i * NY + j] = (i < NX / 2 ? 1 : 2) + (j < NY / 2 ? 0 : 10);
        }

    /* Dip: in column ix = 10, nodes 8..11 go down and up again (a melt-like defect). */
    zDip[10 * NY + 9]  = zMono[10 * NY + 8] - 0.5 * (zMono[10 * NY + 9]  - zMono[10 * NY + 8]);
    zDip[10 * NY + 10] = zMono[10 * NY + 8] - 0.2 * (zMono[10 * NY + 9]  - zMono[10 * NY + 8]);

    /* Plateau: in column ix = 20, nodes 5..8 equal. */
    for (j = 6; j <= 8; j++) zFlat[20 * NY + j] = zFlat[20 * NY + 5];
}

/* Random interior point of the grid, well inside. */
static void randomPoint(double *x, double *y)
{
    *x = exp(log(xAxis[0]) + (0.01 + 0.98 * frand()) * (log(xAxis[NX-1]) - log(xAxis[0])));
    *y = exp(log(yAxis[0]) + (0.01 + 0.98 * frand()) * (log(yAxis[NY-1]) - log(yAxis[0])));
}

/* ------------------------------------------------------------------ tests */

static void testLocate(void)
{
    int failsBefore = nFails;
    INTERPCELL c;
    int i; double w;
    TEST("interpLocateAxis / interpLocate: closed interval, tolerance");

    /* Interior and nodes. */
    for (i = 0; i < NX - 1; i++) {
        int k; double u;
        CHECK(interpLocateAxis(xAxis, NX, xAxis[i], 0.0, &k, &u) == INTERP_SUCCESS && k == i && u == 0.0, "node %d", i);
        CHECK(interpLocateAxis(xAxis, NX, 0.5 * (xAxis[i] + xAxis[i+1]), 0.0, &k, &u) == INTERP_SUCCESS && k == i && fabs(u - 0.5) < 1e-12, "midpoint %d", i);
    }
    /* Upper end is closed: last cell with w = 1. */
    CHECK(interpLocateAxis(xAxis, NX, xAxis[NX-1], 0.0, &i, &w) == INTERP_SUCCESS && i == NX - 2 && w == 1.0, "upper end closed");
    CHECK(interpLocateAxis(xAxis, NX, xAxis[0], 0.0, &i, &w) == INTERP_SUCCESS && i == 0 && w == 0.0, "lower end closed");
    /* Outside, with and without tolerance. */
    CHECK(interpLocateAxis(xAxis, NX, xAxis[0] * (1.0 - 1e-15), 0.0, &i, &w) == INTERP_BELOW, "just below, no tolerance");
    CHECK(interpLocateAxis(xAxis, NX, xAxis[0] * (1.0 - 1e-15), 1e-13, &i, &w) == INTERP_SUCCESS && i == 0 && w == 0.0, "just below, clamped");
    CHECK(interpLocateAxis(xAxis, NX, xAxis[NX-1] * (1.0 + 1e-15), 1e-13, &i, &w) == INTERP_SUCCESS && i == NX - 2 && w == 1.0, "just above, clamped");
    CHECK(interpLocateAxis(xAxis, NX, 0.5 * xAxis[0], 1e-13, &i, &w) == INTERP_BELOW, "far below");
    CHECK(interpLocateAxis(xAxis, NX, 2.0 * xAxis[NX-1], 1e-13, &i, &w) == INTERP_ABOVE, "far above");
    CHECK(interpLocateAxis(xAxis, NX, -INFINITY, 1e-13, &i, &w) == INTERP_BELOW, "-inf");
    CHECK(interpLocateAxis(xAxis, NX, INFINITY, 1e-13, &i, &w) == INTERP_ABOVE, "+inf");
    CHECK(interpLocateAxis(xAxis, NX, NAN, 1e-13, &i, &w) == INTERP_INVALID, "NaN");
    /* The tolerance is relative to the bound, never absolute: on an axis
     * starting at 1e-25 (the ANEOS tables) half the first node is far out. */
    {
        static const double tiny[3] = { 1e-25, 1e-20, 1e-15 };
        CHECK(interpLocateAxis(tiny, 3, 0.5e-25, 1e-13, &i, &w) == INTERP_BELOW, "tiny axis: half the first node is below");
        CHECK(interpLocateAxis(tiny, 3, 1e-25 * (1.0 - 1e-15), 1e-13, &i, &w) == INTERP_SUCCESS && i == 0 && w == 0.0, "tiny axis: roundoff below clamped");
        CHECK(interpLocateAxis(tiny, 3, 1e-15 * (1.0 + 1e-15), 1e-13, &i, &w) == INTERP_SUCCESS && i == 1 && w == 1.0, "tiny axis: roundoff above clamped");
        CHECK(interpLocateAxis(tiny, 3, 2e-15, 1e-13, &i, &w) == INTERP_ABOVE, "tiny axis: twice the last node is above");
    }
    /* Outputs untouched on failure. */
    i = -7; w = -7.0;
    interpLocateAxis(xAxis, NX, 0.5 * xAxis[0], 0.0, &i, &w);
    CHECK(i == -7 && w == -7.0, "outputs touched on failure");

    /* 2D codes: x is checked first. */
    CHECK(interpLocate(&G, 0.5 * xAxis[0], yAxis[3], 0.0, &c) == INTERP_OUTSIDE_XMIN, "XMIN");
    CHECK(interpLocate(&G, 2.0 * xAxis[NX-1], yAxis[3], 0.0, &c) == INTERP_OUTSIDE_XMAX, "XMAX");
    CHECK(interpLocate(&G, xAxis[3], 0.5 * yAxis[0], 0.0, &c) == INTERP_OUTSIDE_YMIN, "YMIN");
    CHECK(interpLocate(&G, xAxis[3], 2.0 * yAxis[NY-1], 0.0, &c) == INTERP_OUTSIDE_YMAX, "YMAX");
    CHECK(interpLocate(&G, 0.5 * xAxis[0], 2.0 * yAxis[NY-1], 0.0, &c) == INTERP_OUTSIDE_XMIN, "x before y");
    CHECK(interpLocate(&G, xAxis[3], NAN, 0.0, &c) == INTERP_INVALID, "NaN y");
    CHECK(interpLocate(&G, xAxis[NX-1], yAxis[NY-1], 0.0, &c) == INTERP_SUCCESS && c.ix == NX - 2 && c.iy == NY - 2 && c.u == 1.0 && c.v == 1.0, "top corner");
    END_TEST();
}

static void testValue(void)
{
    int failsBefore = nFails;
    int n;
    TEST("interpValue reproduces a bilinear function and the nodes");
    for (n = 0; n < 2000; n++) {
        double x, y, z = 0.0;
        randomPoint(&x, &y);
        CHECK(interpValue(&G, zBilin, x, y, 0.0, &z) == INTERP_SUCCESS, "status at (%g,%g)", x, y);
        CHECK_REL(z, fBilin(x, y), 1e-13, "bilinear at (%g,%g)", x, y);
    }
    /* Every node exactly, including all four edges. */
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++) {
            double z = 0.0;
            CHECK(interpValue(&G, zMono, xAxis[i], yAxis[j], 0.0, &z) == INTERP_SUCCESS, "node status (%d,%d)", i, j);
            CHECK(z == zMono[i * NY + j], "node (%d,%d): %.17g vs %.17g", i, j, z, zMono[i * NY + j]);
        }
    /* Interpolation is bounded by the cell's node values. */
    for (n = 0; n < 2000; n++) {
        double x, y, z;
        INTERPCELL c = { 0, 0, 0.0, 0.0 };
        randomPoint(&x, &y);
        interpLocate(&G, x, y, 0.0, &c);
        z = interpCellValue(&G, zMono, &c);
        double f00 = interpZ(&G, zMono, c.ix, c.iy), f11 = interpZ(&G, zMono, c.ix + 1, c.iy + 1);
        CHECK(z >= f00 * (1 - 1e-15) && z <= f11 * (1 + 1e-15), "value outside cell range at (%g,%g)", x, y);
    }
    END_TEST();
}

static void testGradient(void)
{
    int failsBefore = nFails;
    int n;
    TEST("interpGradient is exact on the interpolant");
    for (n = 0; n < 2000; n++) {
        double x, y, dzdx = 0.0, dzdy = 0.0;
        randomPoint(&x, &y);
        CHECK(interpGradient(&G, zBilin, x, y, 0.0, &dzdx, &dzdy) == INTERP_SUCCESS, "status");
        CHECK_REL(dzdx, 2.0 + 0.25 * y, 1e-12, "dz/dx at (%g,%g)", x, y);
        CHECK_REL(dzdy, 0.5 + 0.25 * x, 1e-12, "dz/dy at (%g,%g)", x, y);
    }
    /* On a general table: central differences inside one cell are exact. */
    for (n = 0; n < 500; n++) {
        double x, y, dzdx = 0.0, dzdy = 0.0, zp = 0.0, zm = 0.0;
        INTERPCELL c = { 0, 0, 0.0, 0.0 };
        randomPoint(&x, &y);
        interpLocate(&G, x, y, 0.0, &c);
        if (c.u < 0.05 || c.u > 0.95 || c.v < 0.05 || c.v > 0.95) continue;
        double hx = 0.01 * (xAxis[c.ix + 1] - xAxis[c.ix]), hy = 0.01 * (yAxis[c.iy + 1] - yAxis[c.iy]);
        interpGradient(&G, zMono, x, y, 0.0, &dzdx, &dzdy);
        interpValue(&G, zMono, x + hx, y, 0.0, &zp); interpValue(&G, zMono, x - hx, y, 0.0, &zm);
        CHECK_REL(dzdx, (zp - zm) / (2 * hx), 1e-8, "dz/dx vs FD at (%g,%g)", x, y);
        interpValue(&G, zMono, x, y + hy, 0.0, &zp); interpValue(&G, zMono, x, y - hy, 0.0, &zm);
        CHECK_REL(dzdy, (zp - zm) / (2 * hy), 1e-8, "dz/dy vs FD at (%g,%g)", x, y);
    }
    END_TEST();
}

static void testNearest(void)
{
    int failsBefore = nFails;
    int v;
    TEST("interpNearest picks the nearest node");
    CHECK(interpNearest(&G, zPhase, xAxis[2], yAxis[3], 0.0, &v) == INTERP_SUCCESS && v == 1, "low corner node");
    CHECK(interpNearest(&G, zPhase, xAxis[NX-1], yAxis[NY-1], 0.0, &v) == INTERP_SUCCESS && v == 12, "top corner node");
    /* 0.3 of the way: lower node; 0.7: upper. */
    {
        int i = NX / 2 - 1, j = NY / 2 - 1;
        double x3 = xAxis[i] + 0.3 * (xAxis[i+1] - xAxis[i]), x7 = xAxis[i] + 0.7 * (xAxis[i+1] - xAxis[i]);
        double y3 = yAxis[j] + 0.3 * (yAxis[j+1] - yAxis[j]), y7 = yAxis[j] + 0.7 * (yAxis[j+1] - yAxis[j]);
        CHECK(interpNearest(&G, zPhase, x3, y3, 0.0, &v) == INTERP_SUCCESS && v == 1, "0.3/0.3 -> %d", v);
        CHECK(interpNearest(&G, zPhase, x7, y3, 0.0, &v) == INTERP_SUCCESS && v == 2, "0.7/0.3 -> %d", v);
        CHECK(interpNearest(&G, zPhase, x3, y7, 0.0, &v) == INTERP_SUCCESS && v == 11, "0.3/0.7 -> %d", v);
        CHECK(interpNearest(&G, zPhase, x7, y7, 0.0, &v) == INTERP_SUCCESS && v == 12, "0.7/0.7 -> %d", v);
    }
    CHECK(interpNearest(&G, zPhase, 0.5 * xAxis[0], yAxis[3], 0.0, &v) == INTERP_OUTSIDE_XMIN, "out of range");
    END_TEST();
}

static void testInvertRoundTrip(void)
{
    int failsBefore = nFails;
    static const INTERPINVERTOPT strict = INTERP_INVERT_STRICT;
    static const INTERPINVERTOPT tolerant = INTERP_INVERT_TOLERANT;
    int n;
    TEST("interpInvertX/Y round-trip on a monotonic table, both modes");
    for (n = 0; n < 3000; n++) {
        double x, y, z = 0.0, y2 = 0.0, x2 = 0.0, z2 = 0.0;
        INTERPCELL c = { 0, 0, 0.0, 0.0 }, c2 = { 0, 0, 0.0, 0.0 };
        randomPoint(&x, &y);
        interpValue(&G, zMono, x, y, 0.0, &z);
        interpLocate(&G, x, y, 0.0, &c);

        /* Along y, strict. */
        CHECK(interpInvertY(&G, zMono, x, z, 0.0, &strict, &c2, &y2) == INTERP_SUCCESS, "strict Y status at (%g,%g)", x, y);
        CHECK_REL(y2, y, 1e-12, "strict Y round trip at (%g,%g)", x, y);
        CHECK(c2.ix == c.ix && c2.iy == c.iy, "strict Y cell (%d,%d) vs (%d,%d)", c2.ix, c2.iy, c.ix, c.iy);
        CHECK_REL(c2.v, c.v, 1e-10, "strict Y local coordinate");
        /* The returned cell evaluates the table back to the target. */
        CHECK_REL(interpCellValue(&G, zMono, &c2), z, 1e-13, "strict Y cell value");

        /* Along y, tolerant. */
        CHECK(interpInvertY(&G, zMono, x, z, 0.0, &tolerant, &c2, &y2) == INTERP_SUCCESS, "tolerant Y status");
        CHECK_REL(y2, y, 1e-12, "tolerant Y round trip at (%g,%g)", x, y);

        /* Along x, both modes. */
        CHECK(interpInvertX(&G, zMono, y, z, 0.0, &strict, &c2, &x2) == INTERP_SUCCESS, "strict X status");
        CHECK_REL(x2, x, 1e-12, "strict X round trip at (%g,%g)", x, y);
        CHECK(interpInvertX(&G, zMono, y, z, 0.0, &tolerant, &c2, &x2) == INTERP_SUCCESS, "tolerant X status");
        CHECK_REL(x2, x, 1e-12, "tolerant X round trip at (%g,%g)", x, y);
        /* A second table evaluated in the returned cell agrees with a fresh lookup. */
        interpValue(&G, zBilin, x2, y, 0.0, &z2);
        CHECK_REL(interpCellValue(&G, zBilin, &c2), z2, 1e-12, "second table in returned cell");
    }
    /* Exactly at the line ends, with the strict valTol. */
    {
        INTERPCELL c = { 0, 0, 0.0, 0.0 }; double y2 = 0.0;
        double x = xAxis[7] * 1.37;
        double zlo = 0.0, zhi = 0.0;
        interpValue(&G, zMono, x, yAxis[0], 0.0, &zlo);
        interpValue(&G, zMono, x, yAxis[NY-1], 0.0, &zhi);
        CHECK(interpInvertY(&G, zMono, x, zlo, 0.0, &strict, &c, &y2) == INTERP_SUCCESS && y2 == yAxis[0], "y_min end");
        CHECK(interpInvertY(&G, zMono, x, zhi, 0.0, &strict, &c, &y2) == INTERP_SUCCESS && y2 == yAxis[NY-1], "y_max end");
        CHECK(interpInvertY(&G, zMono, x, zlo * (1 - 1e-15), 0.0, &strict, &c, &y2) == INTERP_SUCCESS && y2 == yAxis[0], "y_min end, roundoff below");
        CHECK(interpInvertY(&G, zMono, x, zhi * (1 + 1e-15), 0.0, &strict, &c, &y2) == INTERP_SUCCESS && y2 == yAxis[NY-1], "y_max end, roundoff above");
        CHECK(interpInvertY(&G, zMono, x, zlo * 0.5, 0.0, &strict, &c, &y2) == INTERP_BELOW, "below the line");
        CHECK(interpInvertY(&G, zMono, x, zhi * 2.0, 0.0, &strict, &c, &y2) == INTERP_ABOVE, "above the line");
        CHECK(interpInvertY(&G, zMono, x, zlo * 0.5, 0.0, &tolerant, &c, &y2) == INTERP_NOT_ATTAINED, "below the line, tolerant");
        CHECK(interpInvertY(&G, zMono, x, zhi * 2.0, 0.0, &tolerant, &c, &y2) == INTERP_NOT_ATTAINED, "above the line, tolerant");
        CHECK(interpInvertY(&G, zMono, 0.5 * xAxis[0], zlo, 0.0, &strict, &c, &y2) == INTERP_OUTSIDE_XMIN, "x outside");
        CHECK(interpInvertY(&G, zMono, x, NAN, 0.0, &strict, &c, &y2) == INTERP_INVALID, "NaN target");
        /* Outputs untouched on failure. */
        c.ix = -7; y2 = -7.0;
        interpInvertY(&G, zMono, x, zlo * 0.5, 0.0, &strict, &c, &y2);
        CHECK(c.ix == -7 && y2 == -7.0, "outputs touched on failure");
    }
    END_TEST();
}

static void testInvertDefects(void)
{
    int failsBefore = nFails;
    static const INTERPINVERTOPT strict = INTERP_INVERT_STRICT;
    static const INTERPINVERTOPT tolerant = INTERP_INVERT_TOLERANT;
    int k = 0; double t = 0.0;
    TEST("inversion on a non-monotonic line and a plateau");

    /* Dip in node column ix = 10. interpNodeLineMonotonic must report it. */
    CHECK(interpNodeLineMonotonic(&G, zDip, INTERP_ALONG_Y, 10) == 9, "dip reported at node %d", interpNodeLineMonotonic(&G, zDip, INTERP_ALONG_Y, 10));
    CHECK(interpNodeLineMonotonic(&G, zDip, INTERP_ALONG_Y, 11) == 0, "neighbouring column monotonic");
    CHECK(interpNodeLineMonotonic(&G, zMono, INTERP_ALONG_X, 5) == 0, "x line monotonic");
    CHECK(interpNodeLineMonotonic(&G, zDip, INTERP_ALONG_X, 9) != 0, "dip visible along x at iy = 9");

    /* Targets in the dip region are attained three times on the node line
     * (u = 0 at ix = 10). Tolerant mode must return a segment that reproduces
     * the target; strict mode either does so or reports NONMONOTONIC. */
    {
        double zlo = zDip[10 * NY + 10], zhi = zDip[10 * NY + 8];   /* dip bottom .. dip top */
        int m;
        for (m = 1; m < 20; m++) {
            double target = zlo + (zhi - zlo) * m / 20.0;
            int r = interpInvertLine(&G, zDip, INTERP_ALONG_Y, 10, 0.0, target, &tolerant, &k, &t);
            CHECK(r == INTERP_SUCCESS, "tolerant on dip: status %d", r);
            if (r == INTERP_SUCCESS) {
                double f = (1 - t) * interpLineValue(&G, zDip, INTERP_ALONG_Y, 10, 0.0, k) + t * interpLineValue(&G, zDip, INTERP_ALONG_Y, 10, 0.0, k + 1);
                CHECK_REL(f, target, 1e-13, "tolerant on dip reproduces target");
                /* The fallback convention: lowest-index root unless the sign bisection landed first. */
            }
            r = interpInvertLine(&G, zDip, INTERP_ALONG_Y, 10, 0.0, target, &strict, &k, &t);
            CHECK(r == INTERP_SUCCESS || r == INTERP_NONMONOTONIC, "strict on dip: status %d", r);
            if (r == INTERP_SUCCESS) {
                double f = (1 - t) * interpLineValue(&G, zDip, INTERP_ALONG_Y, 10, 0.0, k) + t * interpLineValue(&G, zDip, INTERP_ALONG_Y, 10, 0.0, k + 1);
                CHECK_REL(f, target, 1e-13, "strict on dip reproduces target when it succeeds");
            }
        }
    }
    /* Away from the dip the dip column is still inverted correctly by both. */
    {
        double target = 0.5 * (zDip[10 * NY + 15] + zDip[10 * NY + 16]);
        CHECK(interpInvertLine(&G, zDip, INTERP_ALONG_Y, 10, 0.0, target, &strict, &k, &t) == INTERP_SUCCESS && k == 15 && fabs(t - 0.5) < 1e-12, "strict above the dip");
        CHECK(interpInvertLine(&G, zDip, INTERP_ALONG_Y, 10, 0.0, target, &tolerant, &k, &t) == INTERP_SUCCESS && k == 15 && fabs(t - 0.5) < 1e-12, "tolerant above the dip");
    }

    /* Plateau in column ix = 20, nodes 5..8. The plateau value is attained on
     * the whole stretch [y5, y8]; which point is returned depends on the path
     * (the sign bisection stops at the node it hits, k = 4, t = 1, which is
     * y5 again; the strict bisection walks to the top of the plateau), but
     * the result must lie on the plateau and reproduce the target. */
    {
        double zp = zFlat[20 * NY + 5];
        CHECK(interpNodeLineMonotonic(&G, zFlat, INTERP_ALONG_Y, 20) == 6, "plateau reported");
        CHECK(interpInvertLine(&G, zFlat, INTERP_ALONG_Y, 20, 0.0, zp, &tolerant, &k, &t) == INTERP_SUCCESS, "tolerant on plateau");
        CHECK(interpLineCoord(&G, INTERP_ALONG_Y, k, t) >= yAxis[5] * (1 - 1e-15) && interpLineCoord(&G, INTERP_ALONG_Y, k, t) <= yAxis[8] * (1 + 1e-15), "tolerant on plateau: k=%d t=%g off the plateau", k, t);
        CHECK_REL((1 - t) * zFlat[20 * NY + k] + t * zFlat[20 * NY + k + 1], zp, 1e-15, "tolerant on plateau reproduces target");
        CHECK(interpInvertLine(&G, zFlat, INTERP_ALONG_Y, 20, 0.0, zp, &strict, &k, &t) == INTERP_SUCCESS, "strict on plateau");
        CHECK(interpLineCoord(&G, INTERP_ALONG_Y, k, t) >= yAxis[5] * (1 - 1e-15) && interpLineCoord(&G, INTERP_ALONG_Y, k, t) <= yAxis[8] * (1 + 1e-15), "strict on plateau: k=%d t=%g off the plateau", k, t);
        CHECK_REL((1 - t) * zFlat[20 * NY + k] + t * zFlat[20 * NY + k + 1], zp, 1e-15, "strict on plateau reproduces target");
        /* Roundoff above the plateau value but below the next node: nodeTol locates it at the plateau. */
        CHECK(interpInvertLine(&G, zFlat, INTERP_ALONG_Y, 20, 0.0, zp * (1 + 1e-14), &tolerant, &k, &t) == INTERP_SUCCESS, "plateau + roundoff (tolerant)");
        CHECK(interpInvertLine(&G, zFlat, INTERP_ALONG_Y, 20, 0.0, zp * (1 + 1e-14), &strict, &k, &t) == INTERP_SUCCESS, "plateau + roundoff (strict)");
    }

    /* Edge tolerance of the tolerant mode. */
    {
        double g0 = zMono[3 * NY + 0], g1 = zMono[3 * NY + 1];
        double span = g1 - g0;
        CHECK(interpInvertLine(&G, zMono, INTERP_ALONG_Y, 3, 0.0, g0 - 0.5 * tolerant.edgeTol * span, &tolerant, &k, &t) == INTERP_SUCCESS && k == 0 && t == 0.0, "edge overshoot within tolerance");
        CHECK(interpInvertLine(&G, zMono, INTERP_ALONG_Y, 3, 0.0, g0 - 2.0 * tolerant.edgeTol * span, &tolerant, &k, &t) == INTERP_NOT_ATTAINED, "edge overshoot beyond tolerance");
        CHECK(interpInvertLine(&G, zMono, INTERP_ALONG_Y, 3, 0.0, g0 - 0.5 * tolerant.edgeTol * span, &strict, &k, &t) == INTERP_BELOW, "strict has no edge tolerance");
        double gn = zMono[3 * NY + NY - 1], gm = zMono[3 * NY + NY - 2];
        CHECK(interpInvertLine(&G, zMono, INTERP_ALONG_Y, 3, 0.0, gn + 0.5 * tolerant.edgeTol * (gn - gm), &tolerant, &k, &t) == INTERP_SUCCESS && k == NY - 2 && t == 1.0, "top edge overshoot within tolerance");
    }
    END_TEST();
}

static void testNodeDeriv(void)
{
    int failsBefore = nFails;
    double f[NX], fy[NY * 3];
    int i, j;
    TEST("interpNodeDeriv is exact for quadratics (interior) and lines (ends)");
    for (i = 0; i < NX; i++) f[i] = 3.0 - 2.0 * xAxis[i] + 0.7 * xAxis[i] * xAxis[i];
    for (i = 1; i < NX - 1; i++)
        CHECK_REL(interpNodeDeriv(xAxis, NX, f, 1, i), -2.0 + 1.4 * xAxis[i], 1e-10, "interior node %d", i);
    for (i = 0; i < NX; i++) f[i] = 3.0 - 2.0 * xAxis[i];
    CHECK_REL(interpNodeDeriv(xAxis, NX, f, 1, 0), -2.0, 1e-12, "first node");
    CHECK_REL(interpNodeDeriv(xAxis, NX, f, 1, NX - 1), -2.0, 1e-12, "last node");
    /* With a stride: every third entry. */
    for (j = 0; j < NY; j++) fy[3 * j] = 1.0 + yAxis[j] * yAxis[j];
    for (j = 1; j < NY - 1; j++)
        CHECK_REL(interpNodeDeriv(yAxis, NY, fy, 3, j), 2.0 * yAxis[j], 1e-10, "strided node %d", j);
    END_TEST();
}

/* ------------------------------------------------------------------- main */

int main(void)
{
    printf("interpBilinear %s test suite: synthetic %d x %d grid\n", INTERP_VERSION_TEXT, NX, NY);
    setupGrids();

    testLocate();
    testValue();
    testGradient();
    testNearest();
    testInvertRoundTrip();
    testInvertDefects();
    testNodeDeriv();

    printf("%d checks, %d failed\n", nChecks, nFails);
    return nFails;
}