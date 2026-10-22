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
 * that depend on them, since the library's documented response to them
 * (lowest root, or the -1e50 sentinel) is the one being tested.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ANEOSmaterial.h"
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

#define TEST(name) do { currentTest = name; printf("  %-44s", name); fflush(stdout); } while (0)
#define END_TEST() do { printf("%s\n", nFails == failsBefore ? "ok" : "FAILED"); } while (0)
/* Informational line inside a running test; re-indents so END_TEST aligns. */
#define NOTE(...) do { printf("\n    note: "); printf(__VA_ARGS__); printf("\n  %-44s", ""); } while (0)

static int isSentinel(double v) { return v < -1e40; }

/* Sample points: interior of the table, log-spaced in rho and T. */
#define NSAMPLE_RHO 12
#define NSAMPLE_T   12
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

/* 1 if T(rho,u) is uniquely defined by the lowest-root convention at (rho, T):
 * u(rho,T) exceeds every tabulated value below T in this rho cell and u is
 * increasing across the cell containing T. */
static int temperatureInvertible(ANEOSMATERIAL *m, double rho, double T)
{
    double rho_cgs = rho * m->CodeUnitstoCGSforRho;
    int i = findIndex(rho_cgs, m->rhoAxis, m->nRho);
    int jT = findIndex(T, m->TAxis, m->nT);
    double x = (rho_cgs - m->rhoAxis[i]) / (m->rhoAxis[i+1] - m->rhoAxis[i]);
    double hlo = (1.0 - x) * m->uArray[i][jT]   + x * m->uArray[i+1][jT];
    double hhi = (1.0 - x) * m->uArray[i][jT+1] + x * m->uArray[i+1][jT+1];
    if (!(hhi > hlo)) return 0;
    for (int j = 0; j < jT; j++) {
        double h = (1.0 - x) * m->uArray[i][j] + x * m->uArray[i+1][j];
        if (!(h < hlo)) return 0;
    }
    return 1;
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
    /* Interior table nodes must be reproduced exactly (bilinear interpolation).
     * The edge nodes are excluded: the code<->cgs round trip can land one ulp
     * outside the table, which is correctly reported as out of range. */
    for (int i = 1; i < m->nRho - 1; i += (m->nRho - 2) / 7 + 1)
        for (int j = 1; j < m->nT - 1; j += (m->nT - 2) / 7 + 1) {
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
                fprintf(stderr, "    non-invertible at rho = %g g/cc, T = %g K\n", rho * m->CodeUnitstoCGSforRho, T);
            }
            continue;
        }
        double u = ANEOSUofRhoT(m, rho, T);
        double T2 = ANEOSTofRhoU(m, rho, u);
        CHECK_REL(T2, T, 1e-10, "T roundtrip at rho=%g", rho);
        /* P(rho,u) must agree with P(rho,T) at the recovered T. */
        CHECK_REL(ANEOSPofRhoU(m, rho, u), ANEOSPofRhoT(m, rho, T), 1e-10, "P(rho,u) vs P(rho,T) at rho=%g T=%g", rho, T);
        CHECK_REL(ANEOSCofRhoU(m, rho, u), ANEOSCofRhoT(m, rho, T), 1e-10, "c(rho,u) vs c(rho,T) at rho=%g T=%g", rho, T);
    }
    /* T(rho,P): P(T) at fixed rho may have plateaus, so check consistency of P, not T. */
    for (int a = 0; a < NSAMPLE_RHO; a++) for (int b = 0; b < NSAMPLE_T; b++) {
        double rho = sampleRho(m, a), T = sampleT(m, b);
        double P = ANEOSPofRhoT(m, rho, T);
        double T2 = ANEOSTofRhoP(m, rho, P);
        if (isSentinel(T2)) { CHECK(0, "TofRhoP failed at rho=%g T=%g", rho, T); continue; }
        CHECK_REL(ANEOSPofRhoT(m, rho, T2), P, 1e-9, "P(rho,T(rho,P)) at rho=%g T=%g", rho, T);
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
            CHECK_REL(ANEOSPofRhoT(m, rhoP, T), P, 1e-9, "P(rho(P,T),T) at rho=%g T=%g", rho, T);
        if (isSentinel(rhoU)) CHECK(0, "RhoofUT failed at rho=%g T=%g", rho, T);
        else if (rhoU > rhoEdge * (1.0 + 1e-13))
            CHECK_REL(ANEOSUofRhoT(m, rhoU, T), u, 1e-9, "u(rho(u,T),T) at rho=%g T=%g", rho, T);
    }
    /* Exactly at a node the inversion must return that node. */
    {
        int i = m->nRho / 2, j = m->nT / 2;
        double T = m->TAxis[j];
        double P = m->pArray[i][j] * m->CGStoCodeUnitsforP;
        double rho = ANEOSRhoofPT(m, P, T);
        CHECK_REL(rho * m->CodeUnitstoCGSforRho, m->rhoAxis[i], 1e-12, "rho(P,T) at node (%d,%d)", i, j);
    }
    CHECK(m->nLookupFail == 0, "nLookupFail = %ld after density inversions", m->nLookupFail);
    END_TEST();
}

static void testFusedLookup(ANEOSMATERIAL *m)
{
    int failsBefore = nFails;
    TEST("ANEOSPCTofRhoU matches the separate calls");
    for (int a = 0; a < NSAMPLE_RHO; a++) for (int b = 0; b < NSAMPLE_T; b++) {
        double rho = sampleRho(m, a), T = sampleT(m, b);
        double u = ANEOSUofRhoT(m, rho, T);
        double c, Tf;
        double P = ANEOSPCTofRhoU(m, rho, u, &c, &Tf);
        CHECK_REL(P, ANEOSPofRhoU(m, rho, u), 1e-13, "fused P at rho=%g", rho);
        CHECK_REL(c, ANEOSCofRhoU(m, rho, u), 1e-13, "fused c at rho=%g", rho);
        CHECK_REL(Tf, ANEOSTofRhoU(m, rho, u), 1e-13, "fused T at rho=%g", rho);
    }
    END_TEST();
}

static void testDerivatives(ANEOSMATERIAL *m)
{
    int failsBefore = nFails;
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

            /* The finite difference below goes through T(rho,u), whose 1e-10
             * round-trip error is amplified where du/dT is small (cold cells),
             * so the tolerance is looser than for the fixed-T partials. Where
             * the cell has du/dT <= 0 (table defect) the documented sentinel
             * is expected instead. */
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
    END_TEST();
}

static void testIsentrope(ANEOSMATERIAL *m)
{
    int failsBefore = nFails;
    TEST("isentropic evolution conserves tabulated entropy");
    for (int a = 1; a < NSAMPLE_RHO - 1; a++) for (int b = 2; b < NSAMPLE_T - 2; b++) {
        double rho1 = sampleRho(m, a), T1 = sampleT(m, b);
        double u1 = ANEOSUofRhoT(m, rho1, T1);
        double S1 = ANEOSSofRhoU(m, rho1, u1);

        /* Identity: S is always reproduced; u only if S(T) is monotonic at
         * this rho (otherwise the lowest root of S is returned, a table defect). */
        double uSame = ANEOSisentropicU(m, rho1, u1, rho1);
        if (isSentinel(uSame)) { CHECK(0, "isentropicU failed rho=%g->%g", rho1, rho1); continue; }
        CHECK_REL(ANEOSSofRhoU(m, rho1, uSame), S1, 1e-10, "S after isentrope to same rho at rho=%g T=%g", rho1, T1);
        if (entropyMonotonicInColumn(m, findIndex(rho1 * m->CodeUnitstoCGSforRho, m->rhoAxis, m->nRho)))
            CHECK_REL(uSame, u1, 1e-10, "isentrope to same rho at rho=%g T=%g", rho1, T1);

        /* Small compression and expansion: S must be conserved where the state stays in the table. */
        for (int dir = -1; dir <= 1; dir += 2) {
            double rho2 = rho1 * (1.0 + dir * 0.05);
            double u2 = ANEOSisentropicU(m, rho1, u1, rho2);
            if (isSentinel(u2)) { CHECK(0, "isentropicU failed rho=%g->%g", rho1, rho2); continue; }
            double S2 = ANEOSSofRhoU(m, rho2, u2);
            /* The cold-curve clamp is allowed only if S1 is below the table at rho2. */
            double Smin2 = ANEOSSofRhoT(m, rho2, m->TAxis[0]);
            if (S1 >= Smin2) CHECK_REL(S2, S1, 1e-8, "S conserved rho=%g->%g", rho1, rho2);
        }
    }
    END_TEST();
}

static void testRhoofPU(ANEOSMATERIAL *m)
{
    int failsBefore = nFails;
    TEST("rho(P,u) and T(P,u)");
    for (int a = 1; a < NSAMPLE_RHO - 1; a += 3) for (int b = 1; b < NSAMPLE_T - 1; b += 3) {
        double rho = sampleRho(m, a), T = sampleT(m, b);
        double P = ANEOSPofRhoT(m, rho, T), u = ANEOSUofRhoT(m, rho, T);
        double rho2 = ANEOSRhoofPU(m, P, u);
        if (isSentinel(rho2)) { CHECK(0, "RhoofPU failed at rho=%g T=%g", rho, T); continue; }
        /* Possibly another root; the returned state must reproduce (P, u). */
        CHECK_REL(ANEOSPofRhoU(m, rho2, u), P, 1e-7, "P at rho(P,u), rho=%g T=%g", rho, T);
        double T2 = ANEOSTofPU(m, P, u);
        CHECK(!isSentinel(T2) && T2 > 0.0, "TofPU failed at rho=%g T=%g", rho, T);
    }
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

    CHECK(ANEOSUofRhoT(m, rhoLo, 300.0) == -1e50, "UofRhoT below rhoMin");
    CHECK(ANEOSUofRhoT(m, rhoHi, 300.0) == -1e50, "UofRhoT above rhoMax");
    CHECK(ANEOSUofRhoT(m, rhoOk, TLo) == -1e50, "UofRhoT below TMin");
    CHECK(ANEOSUofRhoT(m, rhoOk, THi) == -1e50, "UofRhoT above TMax");
    CHECK(ANEOSPofRhoT(m, rhoHi, 300.0) == -1e50, "PofRhoT above rhoMax");
    CHECK(ANEOSCofRhoT(m, rhoHi, 300.0) == -1e50, "CofRhoT above rhoMax");
    CHECK(ANEOSRhoofPT(m, 1.0, THi) == -1e50, "RhoofPT above TMax");
    CHECK(ANEOSRhoofUT(m, 1.0, TLo) == -1e50, "RhoofUT below TMin");
    CHECK(ANEOSTofRhoP(m, rhoHi, 1.0) == -1e50, "TofRhoP above rhoMax");

    /* TofRhoU falls back to TAxis[0], by contract. */
    CHECK(ANEOSTofRhoU(m, rhoHi, 1.0) == m->TAxis[0], "TofRhoU fallback above rhoMax");
    CHECK(ANEOSTofRhoU(m, rhoOk, uHuge) == m->TAxis[0], "TofRhoU fallback for huge u");

    /* Fused lookup out of rho range. */
    {
        double c, T;
        double P = ANEOSPCTofRhoU(m, rhoHi, 1.0, &c, &T);
        CHECK(P == -1e50 && c == -1e50 && T == m->TAxis[0], "PCTofRhoU out of range: P=%g c=%g T=%g", P, c, T);
    }
    /* Gamma with failed sound speed. */
    if (!m->fixedGamma) CHECK(ANEOSGammaofRhoT(m, rhoHi, 300.0) == -1e50, "Gamma with failed c");

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
        CHECK(ANEOSTmeltofRho(m, -1.0) == -1.0, "Tmelt for negative rho");
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

/* ------------------------------------------------------------------- main */

int main(int argc, char *argv[])
{
    int iMat = MAT_ANEOS_DUNITE;
    double dKpcUnit = 2.06701e-13;   /* Earth radius */
    double dMsolUnit = 4.80438e-08;  /* Earth mass */

    if (argc > 1) iMat = atoi(argv[1]);
    if (argc > 3) { dKpcUnit = atof(argv[2]); dMsolUnit = atof(argv[3]); }

    printf("ANEOSmaterial %s test suite: iMat = %d, dKpcUnit = %g, dMsolUnit = %g\n",
           ANEOS_VERSION_TEXT, iMat, dKpcUnit, dMsolUnit);

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
    testFailureContract(m);   /* last: deliberately increments nLookupFail */

    ANEOSfinalizeMaterial(m);

    printf("%d checks, %d failed\n", nChecks, nFails);
    return nFails;
}