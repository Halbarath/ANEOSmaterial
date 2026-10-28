/*  This file is part of ANEOSmaterial.
 *  Copyright (c) 2020-2026 Thomas Meier & Christian Reinhardt
 *
 *  ANEOSmaterial is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ANEOSmaterial is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ANEOSmaterial.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Test suite for the ANEOSmaterial library.
 *
 * Every test is an invariant that must hold for any valid table, so no
 * reference values are hard-coded. Usage:
 *
 *   testANEOSmaterial [iMat] [dKpcUnit dMsolUnit]
 *
 * Defaults: iMat = 54 (ANEOS dunite), Earth-like code units. Pass negative
 * units to run in cgs. The table file for iMat must be in the working
 * directory. Exit status is the number of failed checks.
 *
 * Known table defects (cells with non-positive heat capacity or decreasing
 * entropy along T) are reported in the first test and skipped by the tests
 * that depend on them: there the library's documented response (some root of
 * a non-monotonic line, or the -1e50 sentinel) is not an invariant.
 */
/* fork() and waitpid() for the fatal paths of the strength-file reader (POSIX). */
#define _POSIX_C_SOURCE 200809L

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <sys/wait.h>
#include "ANEOSmaterial.h"

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

/*
 * Relative tolerance with a floor: the error is measured against max(|a|, |b|, scale). For
 * pressures the scale is pressureScale() below: P is interpolated from table values of both
 * signs (cold tension and thermal pressure), so where it crosses zero its roundoff is relative to
 * the size of those values, not to P itself.
 */
#define CHECK_NEAR(a, b, tol, scale, ...) do { \
    double _a = (a), _b = (b); \
    double _den = fmax(fmax(fabs(_a), fabs(_b)), (scale)); \
    double _err = fabs(_a - _b) / _den; \
    nChecks++; \
    if (!(_err <= (tol))) { \
        nFails++; \
        fprintf(stderr, "FAIL [%s] %s:%d: ", currentTest, __FILE__, __LINE__); \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, " (got %.15e vs %.15e, err %.3e > %.1e relative to %.3e)\n", _a, _b, _err, (double)(tol), _den); \
    } \
} while (0)

#define TEST(name) do { currentTest = name; printf("  %-44s", name); fflush(stdout); } while (0)
#define END_TEST() do { printf("%s\n", nFails == failsBefore ? "ok" : "FAILED"); } while (0)
/* Informational line inside a running test; re-indents so END_TEST aligns. */
#define NOTE(...) do { printf("\n    note: "); printf(__VA_ARGS__); printf("\n  %-44s", ""); } while (0)

static int isSentinel(double v) { return ANEOSFailed(v); }

/*
 * Pressure scale at (rho, T), code units: the largest |P| at the four table nodes of the cell. An
 * interpolated or inverted P near zero is a difference of values of this size, so this is the
 * scale of its roundoff (all P failures with a pure relative tolerance had |P| below ~1e-3 of it).
 */
static double pressureScale(ANEOSMATERIAL *m, double rho, double T)
{
    int i = interpFindIndex(rho * m->CodeUnitstoCGSforRho, m->rhoAxis, m->nRho);
    int j = interpFindIndex(T, m->TAxis, m->nT);
    double s = fmax(fmax(fabs(m->pArray[i][j]), fabs(m->pArray[i+1][j])),
                    fmax(fabs(m->pArray[i][j+1]), fabs(m->pArray[i+1][j+1])));
    return s * m->CGStoCodeUnitsforP;
}

/* Sample points: interior of the table, log-spaced in rho and T. */
/* 100 for a quick run (seconds); 1000 for a full sweep (minutes). */
#define NSAMPLE_RHO 1000
#define NSAMPLE_T   1000
#define SAMPLE_OFFSET 0.0137   /* keeps log-spaced samples off the table nodes */

static double sampleRho(ANEOSMATERIAL *m, int k)
{
    double lo = log(m->rhoAxis[0]), hi = log(m->rhoAxis[m->nRho-1]);
    double f = 0.05 + 0.90 * (k + 0.5 + SAMPLE_OFFSET) / NSAMPLE_RHO;
    return exp(lo + f * (hi - lo)) * m->CGStoCodeUnitsforRho;
}

static double sampleT(ANEOSMATERIAL *m, int k)
{
    double lo = log(m->TAxis[0]), hi = log(m->TAxis[m->nT-1]);
    double f = 0.05 + 0.90 * (k + 0.5 + SAMPLE_OFFSET) / NSAMPLE_T;
    return exp(lo + f * (hi - lo));
}

/* Centre of cell (i, j), in code units / Kelvin. */
static double cellRho(ANEOSMATERIAL *m, int i) { return 0.5 * (m->rhoAxis[i] + m->rhoAxis[i+1]) * m->CGStoCodeUnitsforRho; }
static double cellT(ANEOSMATERIAL *m, int j)   { return 0.5 * (m->TAxis[j] + m->TAxis[j+1]); }

/* 1 if S(T) is strictly increasing along the whole columns i and i+1. */
static int entropyMonotonicInColumn(ANEOSMATERIAL *m, int i)
{
    for (int j = 0; j < m->nT - 1; j++) {
        if (!(m->sArray[i][j+1] > m->sArray[i][j])) return 0;
        if (!(m->sArray[i+1][j+1] > m->sArray[i+1][j])) return 0;
    }
    return 1;
}

/*
 * 1 if the energy u_cgs lies above every defective segment (u(T) not increasing) of the columns
 * i and i+1, so that the u(T) lines through the state's density cell attain u once and only in
 * their monotonic part. A requirement that the whole columns be monotonic excludes every M-ANEOS
 * column (their cold rows are defective at all densities, du/dT < 0 up to 18-58 K), and with it
 * every state; this one only excludes states whose energy reaches into a defect.
 */
static int energyAboveDefectsInColumn(ANEOSMATERIAL *m, int i, double u_cgs)
{
    for (int k = i; k <= i + 1; k++)
        for (int j = 0; j < m->nT - 1; j++)
            if (!(m->uArray[k][j+1] > m->uArray[k][j]) && fmax(m->uArray[k][j], m->uArray[k][j+1]) >= u_cgs) return 0;
    return 1;
}

/* 1 if T(rho,u) has a unique root at (rho, T), so that every search path must
 * return it: u(rho,T) lies strictly inside its T cell, and no other cell of
 * the interpolated line u(rho, .) contains u. */
static int temperatureInvertible(ANEOSMATERIAL *m, double rho, double T)
{
    double rho_cgs = rho * m->CodeUnitstoCGSforRho;
    int i = interpFindIndex(rho_cgs, m->rhoAxis, m->nRho);
    int jT = interpFindIndex(T, m->TAxis, m->nT);
    double x = (rho_cgs - m->rhoAxis[i]) / (m->rhoAxis[i+1] - m->rhoAxis[i]);
    double y = (T - m->TAxis[jT]) / (m->TAxis[jT+1] - m->TAxis[jT]);
    double hlo = (1.0 - x) * m->uArray[i][jT]   + x * m->uArray[i+1][jT];
    double hhi = (1.0 - x) * m->uArray[i][jT+1] + x * m->uArray[i+1][jT+1];
    double u = hlo + y * (hhi - hlo);
    if (!(hlo < u && u < hhi)) return 0;
    double hprev = (1.0 - x) * m->uArray[i][0] + x * m->uArray[i+1][0];
    for (int j = 0; j < m->nT - 1; j++) {
        double hnext = (1.0 - x) * m->uArray[i][j+1] + x * m->uArray[i+1][j+1];
        if (j != jT && fmin(hprev, hnext) <= u && u <= fmax(hprev, hnext)) return 0;
        hprev = hnext;
    }
    return 1;
}

/* Relative tolerance for S after a round trip through u at (rho, T). Solving
 * for the position inside the T cell amplifies roundoff in u by u/du_cell,
 * and that position error maps to S through dS_cell; in cold dense cells the
 * thermal part of u is a tiny fraction of the cold-curve energy, so this can
 * exceed 1e-10 by orders of magnitude without any defect in the library. */
static double entropyRoundTripTol(ANEOSMATERIAL *m, double rho, double T, double S)
{
    double rho_cgs = rho * m->CodeUnitstoCGSforRho;
    int i = interpFindIndex(rho_cgs, m->rhoAxis, m->nRho);
    int j = interpFindIndex(T, m->TAxis, m->nT);
    double x = (rho_cgs - m->rhoAxis[i]) / (m->rhoAxis[i+1] - m->rhoAxis[i]);
    double u0 = (1.0 - x) * m->uArray[i][j]   + x * m->uArray[i+1][j];
    double u1 = (1.0 - x) * m->uArray[i][j+1] + x * m->uArray[i+1][j+1];
    double s0 = (1.0 - x) * m->sArray[i][j]   + x * m->sArray[i+1][j];
    double s1 = (1.0 - x) * m->sArray[i][j+1] + x * m->sArray[i+1][j+1];
    double du = fabs(u1 - u0), ds = fabs(s1 - s0);
    if (!(du > 0.0)) return INFINITY;
    double dt = 2048.0 * DBL_EPSILON * fmax(fabs(u0), fabs(u1)) / du;   /* cell coordinate error */
    return 1e-10 + dt * ds / fmax(fabs(S), ds);
}

/* ------------------------------------------------------------------ tests */

static void testInit(ANEOSMATERIAL *m, int iMat)
{
    int failsBefore = nFails;
    TEST("initialisation and table sanity");
    CHECK(m != NULL, "material is NULL");
    CHECK(m->iMat == iMat, "iMat %d != %d", m->iMat, iMat);
    CHECK(m->nRho >= 2 && m->nT >= 2, "table too small: %d x %d", m->nRho, m->nT);
    CHECK(m->rho0 > 0.0, "rho0 = %g", m->rho0);
    CHECK(ANEOSgetRho0(m) > 0.0, "rho0 in code units = %g", ANEOSgetRho0(m));
    for (int i = 1; i < m->nRho; i++)
        CHECK(m->rhoAxis[i] > m->rhoAxis[i-1], "rhoAxis not increasing at %d", i);
    for (int j = 1; j < m->nT; j++)
        CHECK(m->TAxis[j] > m->TAxis[j-1], "TAxis not increasing at %d", j);
    /* Tables must be finite everywhere. */
    int bad = 0;
    for (int i = 0; i < m->nRho; i++)
        for (int j = 0; j < m->nT; j++)
            if (!isfinite(m->uArray[i][j]) || !isfinite(m->pArray[i][j]) ||
                !isfinite(m->cArray[i][j]) || !isfinite(m->sArray[i][j])) bad++;
    CHECK(bad == 0, "%d non-finite table entries", bad);
    /* Column pointers must index the contiguous block. */
    CHECK(m->uArray[0] == m->uData, "uArray[0] is not uData");
    CHECK(m->uArray[m->nRho-1] == m->uData + (size_t)(m->nRho-1) * m->nT, "uArray column stride wrong");
    CHECK(m->grid.nx == m->nRho && m->grid.ny == m->nT && m->grid.x == m->rhoAxis && m->grid.y == m->TAxis, "interpBilinear grid view inconsistent");
    CHECK(m->CodeUnitstoCGSforU * m->CGStoCodeUnitsforU > 0.0, "unit factors");
    CHECK(m->nLookupFail == 0, "nLookupFail = %ld after init", m->nLookupFail);

    /* Thermodynamic consistency of the nodes: c_v > 0 and dS/dT > 0. Violations
     * are table defects (reported, not failed); the tests below skip them. */
    int nCvFlat = 0, nCvNeg = 0, nSBad = 0;
    for (int i = 0; i < m->nRho; i++)
        for (int j = 0; j < m->nT - 1; j++) {
            double du = m->uArray[i][j+1] - m->uArray[i][j];
            if (du == 0.0) nCvFlat++;
            if (du < 0.0) nCvNeg++;
            if (!(m->sArray[i][j+1] > m->sArray[i][j])) nSBad++;
        }
    if (nCvFlat || nCvNeg || nSBad)
        NOTE("%d cells with du/dT = 0, %d with du/dT < 0, %d with dS/dT <= 0 (table defects)",
             nCvFlat, nCvNeg, nSBad);
    if (nCvNeg) {
        /* Envelope of the defective region: for each density decade, the
         * highest T at which some cell has du/dT < 0. */
        printf("\n    note: highest T with du/dT < 0, per density decade:");
        int iStart = 0;
        while (iStart < m->nRho) {
            double decadeLo = pow(10.0, floor(log10(m->rhoAxis[iStart])));
            double Tbad = 0.0;
            int iEnd = iStart;
            while (iEnd < m->nRho && m->rhoAxis[iEnd] < 10.0 * decadeLo) {
                for (int j = 0; j < m->nT - 1; j++)
                    if (m->uArray[iEnd][j+1] < m->uArray[iEnd][j] && m->TAxis[j+1] > Tbad) Tbad = m->TAxis[j+1];
                iEnd++;
            }
            if (Tbad > 0.0) printf("\n          rho in [%g, %g) g/cc: T <= %g K", decadeLo, 10.0 * decadeLo, Tbad);
            iStart = iEnd;
        }
        printf("\n  %-44s", "");
    }
    END_TEST();
}

static void testForwardLookups(ANEOSMATERIAL *m)
{
    int failsBefore = nFails;
    TEST("forward lookups are finite and consistent");
    for (int a = 0; a < NSAMPLE_RHO; a++) for (int b = 0; b < NSAMPLE_T; b++) {
        double rho = sampleRho(m, a), T = sampleT(m, b);
        double u = ANEOSUofRhoT(m, rho, T);
        double P = ANEOSPofRhoT(m, rho, T);
        double c = ANEOSCofRhoT(m, rho, T);
        double S = ANEOSSofRhoT(m, rho, T);
        CHECK(isfinite(u) && !isSentinel(u), "u(rho=%g,T=%g) = %g", rho, T, u);
        CHECK(isfinite(P) && !isSentinel(P), "P(rho=%g,T=%g) = %g", rho, T, P);
        CHECK(isfinite(c) && !isSentinel(c) && c >= 0.0, "c(rho=%g,T=%g) = %g", rho, T, c);
        CHECK(isfinite(S) && !isSentinel(S), "S(rho=%g,T=%g) = %g", rho, T, S);
    }
    /* Table nodes must be reproduced exactly (bilinear interpolation),
     * including the edge nodes: the intervals are closed and the roundoff
     * tolerance ANEOS_EDGE_TOL absorbs the code<->cgs round trip. */
    for (int i = 0; i < m->nRho; i += (m->nRho - 1) / 7 + 1)
        for (int j = 0; j < m->nT; j += (m->nT - 1) / 7 + 1) {
            double rho = m->rhoAxis[i] * m->CGStoCodeUnitsforRho;
            double T = m->TAxis[j];
            CHECK_REL(ANEOSUofRhoT(m, rho, T) * m->CodeUnitstoCGSforU, m->uArray[i][j], 1e-13, "u at node (%d,%d)", i, j);
            CHECK_REL(ANEOSPofRhoT(m, rho, T) * m->CodeUnitstoCGSforP, m->pArray[i][j], 1e-13, "P at node (%d,%d)", i, j);
        }
    CHECK(m->nLookupFail == 0, "nLookupFail = %ld after in-table lookups", m->nLookupFail);
    END_TEST();
}

static void testTemperatureInversion(ANEOSMATERIAL *m)
{
    int failsBefore = nFails;
    TEST("T(rho,u) inverts u(rho,T)");
    const int maxList = 20;
    int nSkipped = 0, nSkippedWarm = 0;
    for (int a = 0; a < NSAMPLE_RHO; a++) for (int b = 0; b < NSAMPLE_T; b++) {
        double rho = sampleRho(m, a), T = sampleT(m, b);
        /* Several roots (table defect): which one is returned depends on the
         * search path, so the identity is not tested there. Points above 300 K
         * are listed, since those are states a simulation visits. */
        if (!temperatureInvertible(m, rho, T)) {
            nSkipped++;
            if (T > 300.0) {
                nSkippedWarm++;
                if (nSkippedWarm <= maxList)
                    fprintf(stderr, "    non-invertible at rho = %g g/cc, T = %g K\n", rho * m->CodeUnitstoCGSforRho, T);
                else if (nSkippedWarm == maxList + 1)
                    fprintf(stderr, "    ... (further non-invertible points not listed)\n");
            }
            continue;
        }
        double u = ANEOSUofRhoT(m, rho, T);
        double T2 = ANEOSTofRhoU(m, rho, u);
        CHECK_REL(T2, T, 1e-10, "T roundtrip at rho=%g", rho);
        /* P(rho,u) must agree with P(rho,T) at the recovered T. */
        CHECK_NEAR(ANEOSPofRhoU(m, rho, u), ANEOSPofRhoT(m, rho, T), 1e-10, pressureScale(m, rho, T), "P(rho,u) vs P(rho,T) at rho=%g T=%g", rho, T);
        CHECK_REL(ANEOSCofRhoU(m, rho, u), ANEOSCofRhoT(m, rho, T), 1e-10, "c(rho,u) vs c(rho,T) at rho=%g T=%g", rho, T);
    }
    /* T(rho,P): P(T) at fixed rho may have plateaus, so check consistency of P, not T. */
    for (int a = 0; a < NSAMPLE_RHO; a++) for (int b = 0; b < NSAMPLE_T; b++) {
        double rho = sampleRho(m, a), T = sampleT(m, b);
        double P = ANEOSPofRhoT(m, rho, T);
        double T2 = ANEOSTofRhoP(m, rho, P);
        if (isSentinel(T2)) { CHECK(0, "TofRhoP failed at rho=%g T=%g", rho, T); continue; }
        CHECK_NEAR(ANEOSPofRhoT(m, rho, T2), P, 1e-9, pressureScale(m, rho, T), "P(rho,T(rho,P)) at rho=%g T=%g", rho, T);
    }
    CHECK(m->nLookupFail == 0, "nLookupFail = %ld after inversions", m->nLookupFail);
    if (nSkipped)
        NOTE("%d of %d T round trips skipped (u(T) not invertible, table defect), %d of them above 300 K",
             nSkipped, NSAMPLE_RHO * NSAMPLE_T, nSkippedWarm);
    END_TEST();
}

static void testDensityInversion(ANEOSMATERIAL *m)
{
    int failsBefore = nFails;
    TEST("rho(P,T) and rho(u,T) invert the forward lookups");
    double rhoEdge = m->rhoAxis[0] * m->CGStoCodeUnitsforRho;
    for (int a = 0; a < NSAMPLE_RHO; a++) for (int b = 0; b < NSAMPLE_T; b++) {
        double rho = sampleRho(m, a), T = sampleT(m, b);
        double P = ANEOSPofRhoT(m, rho, T);
        double u = ANEOSUofRhoT(m, rho, T);
        double rhoP = ANEOSRhoofPT(m, P, T);
        double rhoU = ANEOSRhoofUT(m, u, T);
        /* Several roots are possible; the returned one must reproduce the
         * input. A result on the lower table edge (plateau in the vapour
         * region, lowest-index convention) is accepted as is: the code<->cgs
         * round trip can land an ulp below rhoAxis[0], which the forward
         * lookup would correctly reject. */
        if (isSentinel(rhoP)) CHECK(0, "RhoofPT failed at rho=%g T=%g", rho, T);
        else if (rhoP > rhoEdge * (1.0 + 1e-13))
            CHECK_NEAR(ANEOSPofRhoT(m, rhoP, T), P, 1e-9, pressureScale(m, rho, T), "P(rho(P,T),T) at rho=%g T=%g", rho, T);
        if (isSentinel(rhoU)) CHECK(0, "RhoofUT failed at rho=%g T=%g", rho, T);
        else if (rhoU > rhoEdge * (1.0 + 1e-13))
            CHECK_REL(ANEOSUofRhoT(m, rhoU, T), u, 1e-9, "u(rho(u,T),T) at rho=%g T=%g", rho, T);
    }
    /*
     * Exactly at a node the inversion must return that node, if that is the only density with this
     * P at this T: on a row j where P(rho) is strictly increasing. On rows with a maximum or plateau
     * (vapour dome, tension) the lowest-density root is returned by convention, which need not be
     * the node; there only the pressure is required. The first monotonic row at or above nT/2 is
     * used, so the exact check runs on every table that has one.
     */
    {
        int i = m->nRho / 2, j, jm = -1;
        for (j = m->nT / 2; j < m->nT && jm < 0; j++) {
            int k, mono = 1;
            for (k = 0; k < m->nRho - 1 && mono; k++) mono = (m->pArray[k+1][j] > m->pArray[k][j]);
            if (mono) jm = j;
        }
        j = (jm >= 0) ? jm : m->nT / 2;
        {
            double T = m->TAxis[j];
            double P = m->pArray[i][j] * m->CGStoCodeUnitsforP;
            double rho = ANEOSRhoofPT(m, P, T);
            if (jm >= 0) {
                CHECK_REL(rho * m->CodeUnitstoCGSforRho, m->rhoAxis[i], 1e-12, "rho(P,T) at node (%d,%d)", i, j);
            } else {
                CHECK(!isSentinel(rho), "RhoofPT failed at node (%d,%d)", i, j);
                if (!isSentinel(rho))
                    CHECK_NEAR(ANEOSPofRhoT(m, rho, T), P, 1e-12, pressureScale(m, rho, T), "P(rho(P,T),T) at node (%d,%d)", i, j);
                NOTE("no row with P(rho) increasing at or above T index %d (table defect): rho(P,T) at node (%d,%d) only reproduces P", m->nT / 2, i, j);
            }
        }
    }
    CHECK(m->nLookupFail == 0, "nLookupFail = %ld after density inversions", m->nLookupFail);
    END_TEST();
}

static void testFusedLookup(ANEOSMATERIAL *m)
{
    int failsBefore = nFails;
    TEST("ANEOSPCTofRhoU matches the separate calls");
    /* The separate calls go through T and recompute the cell coordinate from
     * it, so the two paths differ by roundoff in y, amplified where c or P
     * change steeply across a cell. */
    for (int a = 0; a < NSAMPLE_RHO; a++) for (int b = 0; b < NSAMPLE_T; b++) {
        double rho = sampleRho(m, a), T = sampleT(m, b);
        double u = ANEOSUofRhoT(m, rho, T);
        double c, Tf;
        double P = ANEOSPCTofRhoU(m, rho, u, &c, &Tf);
        CHECK_NEAR(P, ANEOSPofRhoU(m, rho, u), 1e-12, pressureScale(m, rho, T), "fused P at rho=%g", rho);
        CHECK_REL(c, ANEOSCofRhoU(m, rho, u), 1e-12, "fused c at rho=%g", rho);
        CHECK_REL(Tf, ANEOSTofRhoU(m, rho, u), 1e-12, "fused T at rho=%g", rho);
    }
    END_TEST();
}

static void testDerivatives(ANEOSMATERIAL *m)
{
    int failsBefore = nFails;
    int nSkipped = 0;
    TEST("analytic derivatives match the interpolant");
    /* Inside a cell the interpolant is bilinear, so central differences with a
     * step that stays inside the cell are exact for the partials at fixed T
     * and fixed rho, and P is linear in u at fixed rho. */
    for (int i = 1; i < m->nRho - 2; i += (m->nRho - 3) / 8 + 1)
        for (int j = 1; j < m->nT - 2; j += (m->nT - 3) / 8 + 1) {
            double rho = cellRho(m, i), T = cellT(m, j);
            double hr = 0.01 * (m->rhoAxis[i+1] - m->rhoAxis[i]) * m->CGStoCodeUnitsforRho;
            double hT = 0.01 * (m->TAxis[j+1] - m->TAxis[j]);

            double dPdRhoT = ANEOSdPdRhoofRhoT(m, rho, T);
            double fd = (ANEOSPofRhoT(m, rho + hr, T) - ANEOSPofRhoT(m, rho - hr, T)) / (2 * hr);
            CHECK_REL(dPdRhoT, fd, 1e-8, "dP/drho|T in cell (%d,%d)", i, j);

            double dPdT = ANEOSdPdTofRhoT(m, rho, T);
            fd = (ANEOSPofRhoT(m, rho, T + hT) - ANEOSPofRhoT(m, rho, T - hT)) / (2 * hT);
            CHECK_REL(dPdT, fd, 1e-8, "dP/dT|rho in cell (%d,%d)", i, j);

            /* The (rho, u) derivatives go through T(rho,u). Where u(T) is not
             * invertible at this rho (several roots, table defect) the root
             * that comes back depends on roundoff and may lie in a cell with
             * du/dT <= 0, so those states are skipped like in the other tests.
             * Where it is invertible, the 1e-10 round-trip error of T is
             * amplified where du/dT is small (cold cells), so the tolerance is
             * looser than for the fixed-T partials. A cell with du/dT <= 0 on
             * an otherwise invertible line gives the documented sentinel. */
            if (!temperatureInvertible(m, rho, T)) { nSkipped++; continue; }
            double u = ANEOSUofRhoT(m, rho, T);
            double up = ANEOSUofRhoT(m, rho, T + hT), um = ANEOSUofRhoT(m, rho, T - hT);
            double dPdU = ANEOSdPdUofRhoU(m, rho, u);
            if (up > um) {
                fd = (ANEOSPofRhoU(m, rho, up) - ANEOSPofRhoU(m, rho, um)) / (up - um);
                CHECK_REL(dPdU, fd, 1e-6, "dP/du|rho in cell (%d,%d)", i, j);
            } else {
                CHECK(isSentinel(dPdU), "dP/du|rho with du/dT <= 0 in cell (%d,%d) = %g", i, j, dPdU);
            }

            /* At constant u the stencil moves in T as well. The analytic value
             * is the local cell derivative, so only compare when both stencil
             * points stay inside this T cell; in cold dense cells a 1% shift in
             * rho changes u_cold by more than the whole thermal range of the cell. */
            double dPdRhoU = ANEOSdPdRhoofRhoU(m, rho, u);
            double Tp = ANEOSTofRhoU(m, rho + hr, u), Tm = ANEOSTofRhoU(m, rho - hr, u);
            if (!isSentinel(dPdRhoU) &&
                Tp >= m->TAxis[j] && Tp < m->TAxis[j+1] && Tm >= m->TAxis[j] && Tm < m->TAxis[j+1]) {
                fd = (ANEOSPofRhoU(m, rho + hr, u) - ANEOSPofRhoU(m, rho - hr, u)) / (2 * hr);
                CHECK_REL(dPdRhoU, fd, 1e-6, "dP/drho|u in cell (%d,%d)", i, j);
            }

            CHECK_REL(ANEOSdUdRhoofRhoU(m, rho, u), ANEOSPofRhoU(m, rho, u) / (rho * rho), 1e-14, "du/drho in cell (%d,%d)", i, j);
        }
    if (nSkipped) NOTE("%d (rho,u) derivative checks skipped (u(T) not invertible, table defect)", nSkipped);
    END_TEST();
}

static void testIsentrope(ANEOSMATERIAL *m)
{
    int failsBefore = nFails;
    TEST("isentropic evolution conserves tabulated entropy");
    int nSkipped = 0;
    for (int a = 1; a < NSAMPLE_RHO - 1; a++) for (int b = 2; b < NSAMPLE_T - 2; b++) {
        double rho1 = sampleRho(m, a), T1 = sampleT(m, b);
        /* S is checked through S(rho, u), which needs a unique T(rho,u). */
        if (!temperatureInvertible(m, rho1, T1)) { nSkipped++; continue; }
        double u1 = ANEOSUofRhoT(m, rho1, T1);
        double S1 = ANEOSSofRhoU(m, rho1, u1);

        /* Identity: S is always reproduced; u only if S(T) is monotonic at
         * this rho (otherwise another root of S may be returned, a table defect). */
        double uSame = ANEOSisentropicU(m, rho1, u1, rho1);
        if (isSentinel(uSame)) { CHECK(0, "isentropicU failed rho=%g->%g", rho1, rho1); continue; }
        CHECK_REL(ANEOSSofRhoU(m, rho1, uSame), S1, entropyRoundTripTol(m, rho1, T1, S1),
                  "S after isentrope to same rho at rho=%g T=%g", rho1, T1);
        if (entropyMonotonicInColumn(m, interpFindIndex(rho1 * m->CodeUnitstoCGSforRho, m->rhoAxis, m->nRho)))
            CHECK_REL(uSame, u1, 1e-10, "isentrope to same rho at rho=%g T=%g", rho1, T1);

        /* Small compression and expansion: S must be conserved where the
         * isentrope stays in the table and T(rho2, u2) is unique. */
        for (int dir = -1; dir <= 1; dir += 2) {
            double rho2 = rho1 * (1.0 + dir * 0.05);
            INTERPCELL cell;
            double T2;
            static const INTERPINVERTOPT opt = INTERP_INVERT_TOLERANT;
            if (interpInvertY(&m->grid, m->sData, rho2 * m->CodeUnitstoCGSforRho, S1, ANEOS_EDGE_TOL, &opt, &cell, &T2) != INTERP_SUCCESS)
                continue;                                              /* leaves the table: cold-curve fallback by contract */
            if (!temperatureInvertible(m, rho2, T2)) { nSkipped++; continue; }
            double u2 = ANEOSisentropicU(m, rho1, u1, rho2);
            if (isSentinel(u2)) { CHECK(0, "isentropicU failed rho=%g->%g", rho1, rho2); continue; }
            CHECK_REL(ANEOSSofRhoU(m, rho2, u2), S1, entropyRoundTripTol(m, rho2, T2, S1),
                      "S conserved rho=%g->%g T=%g", rho1, rho2, T1);
        }
    }
    if (nSkipped) NOTE("%d isentrope checks skipped (u(T) not invertible, table defect)", nSkipped);

    /*
     * Isentropes that leave the table: a hot state compressed far, or a cold one expanded far, has
     * an entropy that the s(T) line at rho2 does not attain. The state is clamped to the edge it
     * leaves through: the top row u(rho2, T_max) if S lies above the whole line, the cold curve
     * u(rho2, T_min) if below, counted in nLookupFail. (Up to ANEOSmaterial 1.2.0 both went to the
     * cold curve, which dropped the thermal energy of a hot compressed state.) The line's range is
     * taken from the table here, independently of the library, with a margin for the T(rho, u)
     * round trip in S(rho1, u1).
     */
    {
        static const double fac[4] = {1.0 / 50.0, 1.0 / 5.0, 5.0, 50.0};
        const int bs[2] = {2, NSAMPLE_T - 3};
        int nAbove = 0, nBelow = 0;
        long nf0 = m->nLookupFail;
        for (int a = 1; a < NSAMPLE_RHO - 1; a += 7) for (int q = 0; q < 2; q++) for (int f = 0; f < 4; f++) {
            double rho1 = sampleRho(m, a), T1 = sampleT(m, bs[q]), rho2 = fac[f] * rho1;
            double rho2_cgs = rho2 * m->CodeUnitstoCGSforRho, w, smin = INFINITY, smax = -INFINITY;
            int i;
            if (!temperatureInvertible(m, rho1, T1)) continue;
            if (interpLocateAxis(m->rhoAxis, m->nRho, rho2_cgs, ANEOS_EDGE_TOL, &i, &w) != INTERP_SUCCESS) continue;
            for (int j = 0; j < m->nT; j++) {
                double s = (1.0 - w) * m->sArray[i][j] + w * m->sArray[i+1][j];
                smin = fmin(smin, s);
                smax = fmax(smax, s);
            }
            double u1 = ANEOSUofRhoT(m, rho1, T1), S1 = ANEOSSofRhoT(m, rho1, T1);
            double margin = 1e-8 * (fabs(smin) + fabs(smax));
            if (S1 > smax + margin) {
                long nf = m->nLookupFail;
                double u2 = ANEOSisentropicU(m, rho1, u1, rho2);
                CHECK(u2 == ANEOSUofRhoT(m, rho2, m->TAxis[m->nT-1]), "isentrope above the table at rho=%g->%g T=%g: u2=%g, top row %g, cold curve %g",
                      rho1, rho2, T1, u2, ANEOSUofRhoT(m, rho2, m->TAxis[m->nT-1]), ANEOSUCold(m, rho2));
                CHECK(m->nLookupFail == nf + 1, "clamp to the top row counted %ld times", m->nLookupFail - nf);
                nAbove++;
            } else if (S1 < smin - margin) {
                long nf = m->nLookupFail;
                double u2 = ANEOSisentropicU(m, rho1, u1, rho2);
                CHECK(u2 == ANEOSUCold(m, rho2), "isentrope below the table at rho=%g->%g T=%g: u2=%g, cold curve %g",
                      rho1, rho2, T1, u2, ANEOSUCold(m, rho2));
                CHECK(m->nLookupFail == nf + 1, "clamp to the cold curve counted %ld times", m->nLookupFail - nf);
                nBelow++;
            }
        }
        CHECK(nAbove > 0, "no isentrope leaving the table at the top was tested");
        m->nLookupFail = nf0;
        NOTE("%d isentropes clamped to the top row, %d to the cold curve", nAbove, nBelow);
    }
    END_TEST();
}

/* 1 if P(rho', u) is strictly increasing in rho' across the density cell
 * containing rho, sampled at 32 points where u is attainable and at rho itself.
 * Near the cold curve the tabulated P(rho,u) can have a maximum or minimum
 * (large effective Gamma at low T, the tension branch below rho0), and
 * rho(P,u) is then double-valued: not a testable state. rho itself is part of
 * the sequence because the state can sit at such a minimum between two sample
 * points (forsterite at 1.83 g/cc and 35 K: P(rho', u) - P is +1e-5 at the
 * lower node, 0 at rho and positive again above it, a double root that no
 * sign change brackets). */
static int pressureMonotoneInCell(ANEOSMATERIAL *m, double rho, double u)
{
    double rho_cgs = rho * m->CodeUnitstoCGSforRho;
    int i = interpFindIndex(rho_cgs, m->rhoAxis, m->nRho);
    double prev = -INFINITY;
    int bDone = 0;
    for (int k = 0; k <= 32; k++) {
        double r_cgs = m->rhoAxis[i] + (m->rhoAxis[i+1] - m->rhoAxis[i]) * k / 32.0;
        for (int pass = 0; pass < 2; pass++) {
            double r;
            if (pass == 0) {
                /* rho itself, in order, before the first sample above it. */
                if (bDone || r_cgs <= rho_cgs) continue;
                r = rho;
                bDone = 1;
            } else {
                r = r_cgs * m->CGStoCodeUnitsforRho;
            }
            double T = ANEOSTofRhoU(m, r, u);
            if (T == m->TAxis[0]) continue;              /* fallback: u not attainable here */
            double P = ANEOSPofRhoT(m, r, T);
            if (!(P > prev)) return 0;
            prev = P;
        }
    }
    return 1;
}

/*
 * One rho(P,u) / T(P,u) check at the state (rho, T). Returns 0 if checked, 1 if skipped because
 * u(T) is not invertible or u reaches into a defect of the cell's columns (table defects), 2 if
 * skipped because P(rho,u) is not monotone in rho
 * across the cell (double-valued rho(P,u) near the cold curve).
 */
static int checkRhoofPU(ANEOSMATERIAL *m, double rho, double T)
{
    /* The residual P(rho,u) - P needs a unique T(rho,u) in the state's density cell, so a state
     * whose energy reaches into a defective segment of its columns is untestable. */
    int i = interpFindIndex(rho * m->CodeUnitstoCGSforRho, m->rhoAxis, m->nRho);
    if (!temperatureInvertible(m, rho, T)) return 1;
    double P = ANEOSPofRhoT(m, rho, T), u = ANEOSUofRhoT(m, rho, T);
    if (!energyAboveDefectsInColumn(m, i, u * m->CodeUnitstoCGSforU)) return 1;
    if (!pressureMonotoneInCell(m, rho, u)) return 2;
    double rho2 = ANEOSRhoofPU(m, P, u);
    if (isSentinel(rho2)) { CHECK(0, "RhoofPU failed at rho=%g T=%g (density cell %d of %d)", rho, T, i, m->nRho - 1); return 0; }
    /* Possibly another root; the returned state must reproduce (P, u). The root is converged
     * to 1e-10 in density, so the P residual is up to 1e-10 rho dP/drho|u, which relative to P
     * is (K/P) 1e-10 and exceeds 1e-7 where P crosses zero; the bound includes that term. */
    {
        double h = 1e-6 * rho2, P2 = ANEOSPofRhoU(m, rho2, u);
        double slope = fmax(fabs(ANEOSPofRhoU(m, rho2 + h, u) - P2), fabs(P2 - ANEOSPofRhoU(m, rho2 - h, u))) / h;
        CHECK(fabs(P2 - P) <= 1e-7 * fabs(P) + 4.0 * 1e-10 * rho2 * slope, "P at rho(P,u), rho=%g T=%g: %.15e vs %.15e", rho, T, P2, P);
    }
    double T2 = ANEOSTofPU(m, P, u);
    CHECK(!isSentinel(T2) && T2 > 0.0, "TofPU failed at rho=%g T=%g", rho, T);
    return 0;
}

static void testRhoofPU(ANEOSMATERIAL *m)
{
    int failsBefore = nFails;
    TEST("rho(P,u) and T(P,u)");
    int nSkipped = 0, nSkippedBump = 0, nEdge = 0;
    long nf0 = m->nLookupFail;
    for (int a = 1; a < NSAMPLE_RHO - 1; a += 3) for (int b = 1; b < NSAMPLE_T - 1; b += 3) {
        int r = checkRhoofPU(m, sampleRho(m, a), sampleT(m, b));
        if (r == 1) nSkipped++;
        if (r == 2) nSkippedBump++;
    }
    /*
     * The first and the last density cell. The log-spaced samples above cover 5% to 95% of the
     * axis and never reach them, which hid that ANEOSRhoofPU did not scan the last node (a root in
     * the last cell was never bracketed, up to ANEOSmaterial 1.2.0).
     */
    for (int edge = 0; edge < 2; edge++) {
        /* Positions inside the cell off the 1/32 grid of pressureMonotoneInCell, which needs
         * rho to differ from its sample points (P must increase strictly along them). */
        static const double frac[3] = {0.31, 0.53, 0.77};
        int i = edge ? m->nRho - 2 : 0;
        for (int k = 0; k < 3; k++) {
            double rho = (m->rhoAxis[i] + frac[k] * (m->rhoAxis[i+1] - m->rhoAxis[i])) * m->CGStoCodeUnitsforRho;
            for (int b = 1; b < NSAMPLE_T - 1; b += 10) {
                int r = checkRhoofPU(m, rho, sampleT(m, b));
                if (r == 0) nEdge++;
                if (r == 1) nSkipped++;
                if (r == 2) nSkippedBump++;
            }
        }
    }
    CHECK(nEdge > 0, "no rho(P,u) state in the first or last density cell could be tested");
    m->nLookupFail = nf0;   /* the cold-curve probes above may have logged fallbacks */
    if (nSkipped || nSkippedBump)
        NOTE("%d rho(P,u) checks skipped (u(T) not invertible, or u reaches into a table defect), %d skipped (P(rho,u) not monotone in rho near the cold curve)",
             nSkipped, nSkippedBump);
    END_TEST();
}

static void testFailureContract(ANEOSMATERIAL *m)
{
    int failsBefore = nFails;
    TEST("out-of-table inputs return the documented values");
    long nf0 = m->nLookupFail;
    double rhoLo = 0.5 * m->rhoAxis[0] * m->CGStoCodeUnitsforRho;
    double rhoHi = 2.0 * m->rhoAxis[m->nRho-1] * m->CGStoCodeUnitsforRho;
    double rhoOk = sampleRho(m, NSAMPLE_RHO / 2);
    double TLo = 0.5 * m->TAxis[0];
    double THi = 2.0 * m->TAxis[m->nT-1];
    double uHuge = 1e6 * ANEOSUofRhoT(m, rhoOk, m->TAxis[m->nT-1] * 0.999);

    CHECK(ANEOSUofRhoT(m, rhoLo, 300.0) == ANEOS_FAIL_VALUE, "UofRhoT below rhoMin");
    CHECK(ANEOSUofRhoT(m, rhoHi, 300.0) == ANEOS_FAIL_VALUE, "UofRhoT above rhoMax");
    CHECK(ANEOSUofRhoT(m, rhoOk, TLo) == ANEOS_FAIL_VALUE, "UofRhoT below TMin");
    CHECK(ANEOSUofRhoT(m, rhoOk, THi) == ANEOS_FAIL_VALUE, "UofRhoT above TMax");
    CHECK(ANEOSPofRhoT(m, rhoHi, 300.0) == ANEOS_FAIL_VALUE, "PofRhoT above rhoMax");
    CHECK(ANEOSCofRhoT(m, rhoHi, 300.0) == ANEOS_FAIL_VALUE, "CofRhoT above rhoMax");
    CHECK(ANEOSRhoofPT(m, 1.0, THi) == ANEOS_FAIL_VALUE, "RhoofPT above TMax");
    CHECK(ANEOSRhoofUT(m, 1.0, TLo) == ANEOS_FAIL_VALUE, "RhoofUT below TMin");
    CHECK(ANEOSTofRhoP(m, rhoHi, 1.0) == ANEOS_FAIL_VALUE, "TofRhoP above rhoMax");

    /* TofRhoU falls back to TAxis[0], by contract. */
    CHECK(ANEOSTofRhoU(m, rhoHi, 1.0) == m->TAxis[0], "TofRhoU fallback above rhoMax");
    CHECK(ANEOSTofRhoU(m, rhoOk, uHuge) == m->TAxis[0], "TofRhoU fallback for huge u");

    /* Fused lookup out of rho range. */
    {
        double c, T;
        double P = ANEOSPCTofRhoU(m, rhoHi, 1.0, &c, &T);
        CHECK(P == ANEOS_FAIL_VALUE && c == ANEOS_FAIL_VALUE && T == m->TAxis[0], "PCTofRhoU out of range: P=%g c=%g T=%g", P, c, T);
    }
    /* Chained lookups: a failed intermediate gives the sentinel, counted exactly once. */
    {
        long nf = m->nLookupFail;
        CHECK(ANEOSUofPT(m, 1.0, THi) == ANEOS_FAIL_VALUE, "UofPT with rho(P, T) failing (T above TMax)");
        CHECK(m->nLookupFail == nf + 1, "UofPT failure counted once (%ld)", m->nLookupFail - nf);
        CHECK(ANEOSUofRhoP(m, rhoHi, 1.0) == ANEOS_FAIL_VALUE, "UofRhoP with T(rho, P) failing (rho above rhoMax)");
        CHECK(m->nLookupFail == nf + 2, "UofRhoP failure counted once (%ld)", m->nLookupFail - nf);
        CHECK(ANEOSPofUT(m, 1.0, TLo) == ANEOS_FAIL_VALUE, "PofUT with rho(u, T) failing (T below TMin)");
        CHECK(m->nLookupFail == nf + 3, "PofUT failure counted once (%ld)", m->nLookupFail - nf);
        /* On success they are the composition of the two lookups. */
        double rho = rhoOk, T = sampleT(m, NSAMPLE_T / 2), P = ANEOSPofRhoT(m, rho, T), u = ANEOSUofRhoT(m, rho, T);
        CHECK(ANEOSUofPT(m, P, T) == ANEOSUofRhoT(m, ANEOSRhoofPT(m, P, T), T), "UofPT = u(rho(P, T), T)");
        CHECK(ANEOSUofRhoP(m, rho, P) == ANEOSUofRhoT(m, rho, ANEOSTofRhoP(m, rho, P)), "UofRhoP = u(rho, T(rho, P))");
        CHECK(ANEOSPofUT(m, u, T) == ANEOSPofRhoT(m, ANEOSRhoofUT(m, u, T), T), "PofUT = P(rho(u, T), T)");
        CHECK(m->nLookupFail == nf + 3, "successful chained lookups counted (%ld)", m->nLookupFail - nf);
    }
    /* Gamma with failed sound speed. */
    if (!m->fixedGamma) CHECK(ANEOSGammaofRhoT(m, rhoHi, 300.0) == ANEOS_FAIL_VALUE, "Gamma with failed c");

    CHECK(m->nLookupFail > nf0, "nLookupFail did not increase (%ld -> %ld)", nf0, m->nLookupFail);
    END_TEST();
}

static void testOptionalTables(ANEOSMATERIAL *m)
{
    int failsBefore = nFails;
    TEST("phase, melt curve and strength (if present)");
    double rhoOk = sampleRho(m, NSAMPLE_RHO / 2), TOk = sampleT(m, NSAMPLE_T / 2);

    if (m->PhaseArray == NULL) {
        CHECK(ANEOSPhaseofRhoT(m, rhoOk, TOk) == -1, "phase without table must be -1");
    } else {
        int ph = ANEOSPhaseofRhoT(m, rhoOk, TOk);
        CHECK(ph >= 0, "phase in table = %d", ph);
        CHECK(ANEOSPhaseofRhoU(m, rhoOk, ANEOSUofRhoT(m, rhoOk, TOk)) == ph, "phase(rho,u) != phase(rho,T)");
        CHECK(ANEOSPhaseofRhoT(m, 2.0 * m->rhoAxis[m->nRho-1] * m->CGStoCodeUnitsforRho, TOk) == -1, "phase out of table");
    }

    if (m->T_melt == NULL) {
        CHECK(m->yieldStrengthModel < 0, "strength model %d without melt curve", m->yieldStrengthModel);
    } else {
        double rLo = 0.5 * m->rhoAxis[0] * m->CGStoCodeUnitsforRho;
        double rHi = 2.0 * m->rhoAxis[m->nRho-1] * m->CGStoCodeUnitsforRho;
        CHECK(ANEOSTmeltofRho(m, rLo) == m->T_melt[0], "Tmelt clamp below table");
        CHECK(ANEOSTmeltofRho(m, rHi) == m->T_melt[m->nRho-1], "Tmelt clamp above table");
        CHECK(ANEOSFailed(ANEOSTmeltofRho(m, -1.0)), "Tmelt for negative rho");
        double Tm = ANEOSTmeltofRho(m, rhoOk);
        CHECK(Tm >= 0.0 && isfinite(Tm), "Tmelt(rho=%g) = %g", rhoOk, Tm);
        /* Linear interpolation must lie between the neighbouring nodes. */
        int i = m->nRho / 2;
        double rmid = cellRho(m, i);
        double a = fmin(m->T_melt[i], m->T_melt[i+1]), b = fmax(m->T_melt[i], m->T_melt[i+1]);
        double Tmid = ANEOSTmeltofRho(m, rmid);
        CHECK(Tmid >= a && Tmid <= b, "Tmelt mid-cell %g outside [%g,%g]", Tmid, a, b);
    }

    double Y0, YM, mui, mud, xi;
    int model = ANEOSYieldParameters(m, &Y0, &YM, &mui, &mud, &xi);
    CHECK(model == m->yieldStrengthModel, "yield model mismatch");
    if (model > 0) {
        CHECK(Y0 >= 0.0 && YM > Y0 && mui > 0.0 && xi > 0.0, "strength parameters Y0=%g YM=%g mui=%g xi=%g", Y0, YM, mui, xi);
        if (model == 2) CHECK(mud >= 0.0, "model 2 needs mud, got %g", mud);
        double G = ANEOSGammaofRhoT(m, rhoOk, TOk);
        CHECK(G > 0.0 && isfinite(G), "shear modulus %g", G);
    }
    END_TEST();
}

/*
 * Parse a strength file in a child process. Returns 0 if the reader accepted it, 1 if it stopped
 * the program (assert(0), or exit(1) with -DNDEBUG), -1 for anything else.
 */
static int strengthFileRejected(ANEOSMATERIAL *m, const char *path)
{
    fflush(NULL);
    pid_t pid = fork();
    if (pid == 0) {
        if (freopen("/dev/null", "w", stderr) == NULL) _exit(3);   /* the expected diagnostic */
        int r = ANEOSReadStrengthParameters(m, (char *) path);
        _exit(r == 0 ? 0 : 2);
    }
    int status = 0;
    if (pid < 0 || waitpid(pid, &status, 0) != pid) return -1;
    if (WIFEXITED(status) && WEXITSTATUS(status) == 0) return 0;
    if ((WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT) || (WIFEXITED(status) && WEXITSTATUS(status) == 1)) return 1;
    return -1;
}

static void writeFile(const char *path, const char *text)
{
    FILE *fp = fopen(path, "w");
    if (fp == NULL) { perror(path); exit(1); }
    fputs(text, fp);
    fclose(fp);
}

/*
 * The strength parameters are validated when they are read, so a bad file stops the run at
 * initialisation and not at the first strength evaluation: YM > Y0 (the intact strength has a
 * pole at P = (Y0 - YM)/mui otherwise), mui > 0, xi > 0, a positive fixed Gamma, mud for model 2,
 * a known model and finite values ("nan" and "inf" are numbers for strtod).
 */
static void testStrengthFileValidation(ANEOSMATERIAL *m)
{
    int failsBefore = nFails;
    TEST("strength file validation at load");
    const char *path = "testANEOSmaterial_strength.tmp";
    static const struct { const char *text; int bRejected; const char *what; } cases[] = {
        {"yieldStrengthModel = 1\nY0 = 1e8\nYM = 3.5e10\nmui = 1.5\nxi = 1.2\nGamma = -1\n", 0, "model 1, no mud"},
        {"yieldStrengthModel = 2\nY0 = 1e8\nYM = 3.5e10\nmui = 1.5\nmud = 0.6\nxi = 1.2\nGamma = 3e11\n", 0, "model 2, fixed Gamma"},
        {"yieldStrengthModel = 2\nY0 = 0\nYM = 3.5e10\nmui = 1.5\nmud = 0\nxi = 1.2\nGamma = -1\n", 0, "model 2, Y0 = 0, mud = 0"},
        {"yieldStrengthModel = 0\nY0 = 1e8\nYM = 1e8\nmui = 0\nxi = 0\nGamma = 0\n", 0, "model 0 (no limiter): relations not needed"},
        {"yieldStrengthModel = 1\nY0 = 3.5e10\nYM = 3.5e10\nmui = 1.5\nxi = 1.2\nGamma = -1\n", 1, "YM == Y0"},
        {"yieldStrengthModel = 1\nY0 = 3.5e10\nYM = 1e8\nmui = 1.5\nxi = 1.2\nGamma = -1\n", 1, "YM < Y0"},
        {"yieldStrengthModel = 1\nY0 = 1e8\nYM = 3.5e10\nmui = 0\nxi = 1.2\nGamma = -1\n", 1, "mui = 0"},
        {"yieldStrengthModel = 1\nY0 = 1e8\nYM = 3.5e10\nmui = 1.5\nxi = 0\nGamma = -1\n", 1, "xi = 0"},
        {"yieldStrengthModel = 1\nY0 = 1e8\nYM = 3.5e10\nmui = 1.5\nxi = 1.2\nGamma = 0\n", 1, "fixed Gamma = 0"},
        {"yieldStrengthModel = 2\nY0 = 1e8\nYM = 3.5e10\nmui = 1.5\nxi = 1.2\nGamma = -1\n", 1, "model 2 without mud"},
        {"yieldStrengthModel = 2\nY0 = 1e8\nYM = 3.5e10\nmui = 1.5\nmud = -0.5\nxi = 1.2\nGamma = -1\n", 1, "model 2, mud < 0"},
        {"yieldStrengthModel = 3\nY0 = 1e8\nYM = 3.5e10\nmui = 1.5\nxi = 1.2\nGamma = -1\n", 1, "unknown model 3"},
        {"yieldStrengthModel = 1\nY0 = nan\nYM = 3.5e10\nmui = 1.5\nxi = 1.2\nGamma = -1\n", 1, "Y0 = nan"},
        {"yieldStrengthModel = 1\nY0 = 1e8\nYM = inf\nmui = 1.5\nxi = 1.2\nGamma = -1\n", 1, "YM = inf"},
        {"yieldStrengthModel = 2\nY0 = 1e8\nYM = 3.5e10\nmui = 1.5\nmud = nan\nxi = 1.2\nGamma = -1\n", 1, "mud = nan"},
        {"yieldStrengthModel = 1\nY0 = 1e8\nYM = 3.5e10\nmui = 1.5\nxi = 1.2\n", 1, "Gamma missing"},
    };
    for (size_t k = 0; k < sizeof(cases) / sizeof(cases[0]); k++) {
        writeFile(path, cases[k].text);
        int r = strengthFileRejected(m, path);
        CHECK(r == cases[k].bRejected, "%s: %s (got %d)", cases[k].what, cases[k].bRejected ? "must be rejected" : "must be accepted", r);
    }
    /* An accepted file is read into the material (parsed here in the parent, then restored). */
    {
        ANEOSMATERIAL save = *m;
        writeFile(path, cases[1].text);
        CHECK(ANEOSReadStrengthParameters(m, (char *) path) == 0, "valid file read");
        CHECK(m->yieldStrengthModel == 2 && m->Y0 == 1e8 && m->YM == 3.5e10 && m->mui == 1.5 && m->mud == 0.6 &&
              m->xi == 1.2 && m->fixedGamma == 1 && m->Gamma == 3e11, "values of the valid file");
        *m = save;
    }
    remove(path);
    END_TEST();
}

/*
 * Load a table file in a child process. Returns 0 if ANEOSinitMaterialFromFile accepted it, 1 if
 * it stopped the program (assert(0), or exit(1) with -DNDEBUG), -1 for anything else.
 */
static int tableFileRejected(const char *path)
{
    fflush(NULL);
    pid_t pid = fork();
    if (pid == 0) {
        if (freopen("/dev/null", "w", stderr) == NULL) _exit(3);   /* the expected diagnostic */
        ANEOSMATERIAL *t = ANEOSinitMaterialFromFile(54, (char *) path, -1.0, -1.0);
        ANEOSfinalizeMaterial(t);
        _exit(0);
    }
    int status = 0;
    if (pid < 0 || waitpid(pid, &status, 0) != pid) return -1;
    if (WIFEXITED(status) && WEXITSTATUS(status) == 0) return 0;
    if ((WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT) || (WIFEXITED(status) && WEXITSTATUS(status) == 1)) return 1;
    return -1;
}

/* A small table in the file format of ANEOSinitMaterialFromFile (P, u, s, c stored T-major). */
#define VT_NRHO 4
#define VT_NT   3
struct validTable {
    double rho0;
    double rho[VT_NRHO], T[VT_NT];
    double p[VT_NT][VT_NRHO], u[VT_NT][VT_NRHO], s[VT_NT][VT_NRHO], c[VT_NT][VT_NRHO];
};

static void validTableInit(struct validTable *t)
{
    t->rho0 = 2.0;
    for (int i = 0; i < VT_NRHO; i++) t->rho[i] = 1.0 + i;
    for (int j = 0; j < VT_NT; j++) t->T[j] = pow(10.0, 1.0 + j);
    for (int j = 0; j < VT_NT; j++) for (int i = 0; i < VT_NRHO; i++) {
        t->p[j][i] = 1e7 * t->rho[i] * t->T[j];
        t->u[j][i] = 1e7 * t->T[j] + 1e9 * t->rho[i];
        t->s[j][i] = 1e7 * log(t->T[j]) - 1e7 * log(t->rho[i]) + 1e8;
        t->c[j][i] = 1e5 * sqrt(t->T[j]);
    }
}

static void validTableWrite(const struct validTable *t, const char *path)
{
    FILE *fp = fopen(path, "wb");
    int nRho = VT_NRHO, nT = VT_NT;
    if (fp == NULL) { perror(path); exit(1); }
    fwrite(&t->rho0, sizeof(double), 1, fp);
    fwrite(&nRho, sizeof(int), 1, fp);
    fwrite(&nT, sizeof(int), 1, fp);
    fwrite(t->rho, sizeof(double), VT_NRHO, fp);
    fwrite(t->T, sizeof(double), VT_NT, fp);
    fwrite(t->p, sizeof(double), VT_NRHO * VT_NT, fp);
    fwrite(t->u, sizeof(double), VT_NRHO * VT_NT, fp);
    fwrite(t->s, sizeof(double), VT_NRHO * VT_NT, fp);
    fwrite(t->c, sizeof(double), VT_NRHO * VT_NT, fp);
    fclose(fp);
}

static void meltCurveWrite(const char *path, const double *Tm)
{
    FILE *fp = fopen(path, "wb");
    int nRho = VT_NRHO;
    if (fp == NULL) { perror(path); exit(1); }
    fwrite(&nRho, sizeof(int), 1, fp);
    fwrite(Tm, sizeof(double), VT_NRHO, fp);
    fclose(fp);
}

/*
 * The table is validated at load: rho0 > 0, axes positive and strictly increasing (the
 * interpBilinear core relies on it), every entry of P, u, s, c finite, the melt curve finite and
 * >= 0. A defective table stops the run at initialisation instead of producing NaNs later.
 */
static void testTableValidation(void)
{
    int failsBefore = nFails;
    TEST("table validation at load");
    const char *path = "testANEOSmaterial_table.tmp";
    const char *melt = "testANEOSmaterial_table.tmp_tmelt";
    struct validTable t;

    /* A valid table loads and evaluates. */
    validTableInit(&t);
    validTableWrite(&t, path);
    CHECK(tableFileRejected(path) == 0, "valid table must be accepted");
    {
        ANEOSMATERIAL *v = ANEOSinitMaterialFromFile(54, (char *) path, -1.0, -1.0);
        CHECK(v->nRho == VT_NRHO && v->nT == VT_NT, "valid table: sizes");
        CHECK_REL(ANEOSPofRhoT(v, 1.5, 50.0), 1e7 * 1.5 * 50.0, 1e-14, "valid table: bilinear P");
        ANEOSfinalizeMaterial(v);
    }

    /* One defect at a time. */
    for (int k = 0; k < 14; k++) {
        const char *what = "";
        validTableInit(&t);
        switch (k) {
            case 0:  t.rho0 = 0.0;              what = "rho0 = 0"; break;
            case 1:  t.rho0 = NAN;              what = "rho0 = NaN"; break;
            case 2:  t.rho[2] = t.rho[1];       what = "rho axis not strictly increasing"; break;
            case 3:  t.rho[3] = 0.5;            what = "rho axis decreasing at the end"; break;
            case 4:  t.rho[0] = 0.0;            what = "rho axis starts at 0"; break;
            case 5:  t.rho[3] = INFINITY;       what = "rho axis infinite"; break;
            case 6:  t.T[1] = t.T[2];           what = "T axis not strictly increasing"; break;
            case 7:  t.T[0] = -1.0;             what = "T axis negative"; break;
            case 8:  t.T[1] = NAN;              what = "T axis NaN"; break;
            case 9:  t.p[1][2] = NAN;           what = "P NaN"; break;
            case 10: t.u[2][3] = INFINITY;      what = "u infinite"; break;
            case 11: t.s[0][0] = NAN;           what = "s NaN"; break;
            case 12: t.c[2][0] = -INFINITY;     what = "c -inf"; break;
            case 13: t.p[2][3] = NAN;           what = "P NaN in the last entry"; break;
        }
        validTableWrite(&t, path);
        CHECK(tableFileRejected(path) == 1, "%s: table must be rejected", what);
    }

    /* Melt curve: 0 is allowed (no condensed phase), negative or non-finite is not. */
    validTableInit(&t);
    validTableWrite(&t, path);
    {
        double Tm[VT_NRHO] = {0.0, 1500.0, 2000.0, 2500.0};
        meltCurveWrite(melt, Tm);
        CHECK(tableFileRejected(path) == 0, "melt curve with T_melt = 0 must be accepted");
        Tm[2] = NAN;
        meltCurveWrite(melt, Tm);
        CHECK(tableFileRejected(path) == 1, "melt curve with NaN must be rejected");
        Tm[2] = -1.0;
        meltCurveWrite(melt, Tm);
        CHECK(tableFileRejected(path) == 1, "negative melt curve must be rejected");
    }
    remove(melt);
    remove(path);
    END_TEST();
}

/*
 * Initialise a material with the unit pair (dKpcUnit, dMsolUnit) in a child process. Returns 0 if
 * the library accepted it, 1 if it stopped the program (assert(0), or exit(1) with -DNDEBUG), -1
 * for anything else.
 */
static int initUnitsRejected(int iMat, double dKpcUnit, double dMsolUnit)
{
    fflush(NULL);
    pid_t pid = fork();
    if (pid == 0) {
        if (freopen("/dev/null", "w", stderr) == NULL) _exit(3);   /* the expected diagnostic */
        ANEOSMATERIAL *m = ANEOSinitMaterial(iMat, dKpcUnit, dMsolUnit);
        ANEOSfinalizeMaterial(m);
        _exit(0);
    }
    int status = 0;
    if (pid < 0 || waitpid(pid, &status, 0) != pid) return -1;
    if (WIFEXITED(status) && WEXITSTATUS(status) == 0) return 0;
    if ((WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT) || (WIFEXITED(status) && WEXITSTATUS(status) == 1)) return 1;
    return -1;
}

/*
 * One rule for the unit pair in EOSlib and every sublibrary: both units > 0 (and finite) are code
 * units, both <= 0 is cgs, anything else is fatal. (Before, exactly one unit <= 0 made igeos divide
 * by zero, Tillotson assert, and ANEOS, SCvH and REOS3 silently use cgs.)
 */
static void testUnitPair(int iMat)
{
    int failsBefore = nFails;
    TEST("unit pair: code units, cgs or fatal");
    static const struct { double dKpcUnit, dMsolUnit; int bRejected; const char *what; } cases[] = {
        {2.06701e-13, 4.80438e-08, 0, "both > 0: code units"},
        {0.0, 0.0, 0, "both 0: cgs"},
        {-1.0, -1.0, 0, "both < 0: cgs"},
        {0.0, -1.0, 0, "0 and < 0: cgs"},
        {2.06701e-13, 0.0, 1, "dMsolUnit = 0 only"},
        {0.0, 4.80438e-08, 1, "dKpcUnit = 0 only"},
        {-1.0, 4.80438e-08, 1, "dKpcUnit < 0 only"},
        {2.06701e-13, NAN, 1, "dMsolUnit NaN"},
        {INFINITY, 4.80438e-08, 1, "dKpcUnit inf"},
        {NAN, NAN, 1, "both NaN"},
    };
    for (size_t k = 0; k < sizeof(cases) / sizeof(cases[0]); k++) {
        int r = initUnitsRejected(iMat, cases[k].dKpcUnit, cases[k].dMsolUnit);
        CHECK(r == cases[k].bRejected, "%s: %s (got %d)", cases[k].what, cases[k].bRejected ? "must be fatal" : "must be accepted", r);
    }
    END_TEST();
}

/*
 * The accessors EOSlib uses instead of reading the axes: ANEOSTmin, ANEOSUCold and ANEOSIsInTable.
 * The status must agree with the lookups: inside means T(rho, u) inverts without the T_min
 * fallback (on an invertible column), outside in rho means the forward lookups fail.
 */
static void testIsInTable(ANEOSMATERIAL *m)
{
    int failsBefore = nFails;
    TEST("ANEOSIsInTable, ANEOSUCold, ANEOSTmin");
    const double rhoMin = m->rhoAxis[0] * m->CGStoCodeUnitsforRho, rhoMax = m->rhoAxis[m->nRho-1] * m->CGStoCodeUnitsforRho;
    const double TMin = m->TAxis[0], TMax = m->TAxis[m->nT-1];
    long nf0 = m->nLookupFail;

    CHECK(ANEOSTmin(m) == TMin, "ANEOSTmin");
    for (int a = 0; a < NSAMPLE_RHO; a += 7) {
        double rho = sampleRho(m, a);
        double uc = ANEOSUCold(m, rho), uTop = ANEOSUofRhoT(m, rho, TMax);
        CHECK(uc == ANEOSUofRhoT(m, rho, TMin), "UCold = u(rho, T_min) at rho=%g", rho);
        CHECK(ANEOSIsInTable(m, rho, uc) == ANEOS_SUCCESS, "on the cold curve at rho=%g", rho);
        CHECK(ANEOSIsInTable(m, rho, uTop) == ANEOS_SUCCESS, "on the top row at rho=%g", rho);
        CHECK(ANEOSIsInTable(m, rho, uc - 1e-6 * fabs(uc) - 1e-6 * fabs(uTop - uc)) == ANEOS_OUTSIDE_TMIN, "below the cold curve at rho=%g", rho);
        CHECK(ANEOSIsInTable(m, rho, uTop + 1e-6 * fabs(uTop) + 1e-6 * fabs(uTop - uc)) == ANEOS_OUTSIDE_TMAX, "above the top row at rho=%g", rho);
        CHECK(ANEOSIsInTable(m, rho, NAN) == ANEOS_OUTSIDE_TMIN, "NaN u at rho=%g", rho);
        /* Inside on an invertible column: T(rho, u) needs no fallback. */
        for (int b = 0; b < NSAMPLE_T; b += 97) {
            double T = sampleT(m, b), u = ANEOSUofRhoT(m, rho, T);
            if (!temperatureInvertible(m, rho, T)) continue;
            long nf = m->nLookupFail;
            CHECK(ANEOSIsInTable(m, rho, u) == ANEOS_SUCCESS, "in-table state at rho=%g T=%g", rho, T);
            ANEOSTofRhoU(m, rho, u);
            CHECK(m->nLookupFail == nf, "in-table state needed the T fallback at rho=%g T=%g", rho, T);
        }
    }
    /* Edges closed (with the roundoff tolerance), outside beyond. */
    CHECK(ANEOSIsInTable(m, rhoMin, ANEOSUCold(m, rhoMin)) == ANEOS_SUCCESS, "on rhoMin");
    CHECK(ANEOSIsInTable(m, rhoMax, ANEOSUCold(m, rhoMax)) == ANEOS_SUCCESS, "on rhoMax");
    CHECK(m->nLookupFail == nf0, "status queries and in-table lookups counted (%ld -> %ld)", nf0, m->nLookupFail);
    CHECK(ANEOSIsInTable(m, 0.5 * rhoMin, 1.0) == ANEOS_OUTSIDE_RHOMIN, "below rhoMin");
    CHECK(ANEOSIsInTable(m, 2.0 * rhoMax, 1.0) == ANEOS_OUTSIDE_RHOMAX, "above rhoMax");
    CHECK(ANEOSIsInTable(m, NAN, 1.0) == ANEOS_OUTSIDE_RHOMIN, "NaN rho");
    CHECK(ANEOSIsInTable(m, -1.0, 1.0) == ANEOS_OUTSIDE_RHOMIN, "negative rho");
    CHECK(m->nLookupFail == nf0, "a status query is not a failure (%ld -> %ld)", nf0, m->nLookupFail);
    /* UCold outside the table: the sentinel, counted. */
    CHECK(ANEOSUCold(m, 2.0 * rhoMax) == ANEOS_FAIL_VALUE, "UCold above rhoMax");
    CHECK(m->nLookupFail == nf0 + 1, "UCold failure counted (%ld -> %ld)", nf0, m->nLookupFail);
    m->nLookupFail = nf0;
    END_TEST();
}

/* ------------------------------------------------------------------- main */

int main(int argc, char *argv[])
{
    int iMat = MAT_ANEOS_DUNITE;
    double dKpcUnit = 2.06701e-13;   /* Earth radius */
    double dMsolUnit = 4.80438e-08;  /* Earth mass */

    if (argc > 1) iMat = atoi(argv[1]);
    if (argc > 3) { dKpcUnit = atof(argv[2]); dMsolUnit = atof(argv[3]); }

    printf("ANEOSmaterial %s (interpBilinear %s) test suite: iMat = %d, dKpcUnit = %g, dMsolUnit = %g\n",
           ANEOS_VERSION_TEXT, INTERP_VERSION_TEXT, iMat, dKpcUnit, dMsolUnit);

    ANEOSMATERIAL *m = ANEOSinitMaterial(iMat, dKpcUnit, dMsolUnit);
    printf("  table %s: nRho = %d, nT = %d\n", m->matName, m->nRho, m->nT);

    testInit(m, iMat);
    testForwardLookups(m);
    testTemperatureInversion(m);
    testDensityInversion(m);
    testFusedLookup(m);
    testDerivatives(m);
    testIsentrope(m);
    testRhoofPU(m);
    testOptionalTables(m);
    testIsInTable(m);
    testStrengthFileValidation(m);
    testTableValidation();
    testUnitPair(iMat);
    testFailureContract(m);   /* last: deliberately increments nLookupFail */

    ANEOSfinalizeMaterial(m);

    printf("%d checks, %d failed\n", nChecks, nFails);
    return nFails;
}