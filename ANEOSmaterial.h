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
 * ANEOS material library
 * Header of ANEOSmaterial.c
 *
 */

#ifndef ANEOSMATERIAL_HINCLUDED
#define ANEOSMATERIAL_HINCLUDED
#include <stdio.h>
#include "interpBilinear/interpBilinear.h"   /* git subtree of the interpBilinear repository */

/*
 * Version history
 *   1.0.4  ANEOSRhoofPU: brackets roots next to the cold curve, rejects edge
 *          artefacts and keeps scanning past them; fused ANEOSPCTofRhoU;
 *          failure sentinel is exactly -1e50 in every unit system;
 *          interpolateValueNearest rejects NaN instead of aborting.
 *   1.0.5  interpolation through the header-only interpBilinear core (git
 *          subtree in interpBilinear/, replaces interpBilinear.c/.h); the
 *          table intervals are closed at the upper edge (rho == rhoAxis[nRho-1]
 *          and T == TAxis[nT-1] are inside) with the roundoff tolerance
 *          ANEOS_EDGE_TOL; the core's status codes are mapped to the ANEOS
 *          contract in ANEOSmaterial.c. Same search strategy and tolerances
 *          as before (INTERP_INVERT_TOLERANT), nothing else in the numerics changed.
 *   1.2.0  EOSlib failure contract spelled out: ANEOS_FAIL_VALUE and
 *          ANEOSFailed() replace the literal -1e50, ANEOS_SUCCESS/ANEOS_FAIL
 *          status codes; ANEOSTmeltofRho returns the sentinel for an
 *          invalid rho (a missing melt curve is fatal). Numerics unchanged.
 *          ANEOSReadStrengthParameters validates the parameters at load (fatal):
 *          finite values, model 0, 1 or 2, and for models 1/2 YM > Y0, mui > 0,
 *          xi > 0, a positive fixed Gamma, and mud >= 0 for model 2.
 *          The table is validated at load (fatal): rho0 > 0, both axes positive
 *          and strictly increasing, every entry of P, u, s and c finite, the
 *          melt curve finite and >= 0 (SCvH and REOS3 already checked theirs).
 *          An inconsistent unit pair (exactly one of dKpcUnit, dMsolUnit <= 0,
 *          NaN, inf) is fatal; it silently meant cgs. New ANEOSIsInTable,
 *          ANEOSTmin and ANEOSUCold, so that EOSlib no longer reads the axes.
 *          ANEOSUofPT, ANEOSUofRhoP and ANEOSPofUT return ANEOS_FAIL_VALUE when
 *          their intermediate lookup fails (counted once) instead of passing the
 *          sentinel on to a second lookup (counted twice). Requires interpBilinear
 *          1.0.3 with INTERP_ABI_VERSION 1 (ANEOS_INTERP_*, prefixed: the
 *          unprefixed INTERP_*_REQUIRED clashed with scvheos.h and reos3.h).
 *          ANEOSPrintMat prints nLookupFail. ANEOSTmeltofRho no longer has an
 *          unreachable branch returning -1.0 (outside the failure contract);
 *          its comment now says what it does (fatal without a melt curve,
 *          ANEOS_FAIL_VALUE for an invalid rho, clamped to the table).
 *          ANEOSRhoofPU scans every density node including the last one: a root in
 *          the last density cell was never bracketed and gave ANEOS_FAIL_VALUE.
 *          ANEOSisentropicU clamps an isentrope that leaves the table at the top to
 *          u(rho2, T_max), like scvheosIsentropicU; it was clamped to the cold curve
 *          like one leaving at the bottom, which dropped the thermal energy of a hot
 *          compressed state.
 */
#define ANEOS_VERSION_TEXT    "1.2.0"
#define ANEOS_VERSION_MAJOR   1
#define ANEOS_VERSION_MINOR   2
#define ANEOS_VERSION_PATCH   0

/*
 * interpBilinear: the minimum version and the layout of its types this library needs. The names
 * carry the library prefix because ANEOSmaterial, scvheos and reos3 all vendor the core and their
 * headers meet in one translation unit (EOSlib), where only the first copy of the core is
 * included (see INTERP_ABI_VERSION in interpBilinear.h). Checked at compile time.
 */
#define ANEOS_INTERP_VERSION_MAJOR_REQUIRED 1
#define ANEOS_INTERP_VERSION_MINOR_REQUIRED 0
#define ANEOS_INTERP_VERSION_PATCH_REQUIRED 3
#define ANEOS_INTERP_ABI_VERSION 1
#if (INTERP_VERSION_MAJOR * 10000 + INTERP_VERSION_MINOR * 100 + INTERP_VERSION_PATCH) < \
    (ANEOS_INTERP_VERSION_MAJOR_REQUIRED * 10000 + ANEOS_INTERP_VERSION_MINOR_REQUIRED * 100 + ANEOS_INTERP_VERSION_PATCH_REQUIRED)
#error "ANEOSmaterial: the interpBilinear in effect is older than the required version (in EOSlib: the copy of the first library header included)"
#endif
#if !defined(INTERP_ABI_VERSION) || (INTERP_ABI_VERSION != ANEOS_INTERP_ABI_VERSION)
#error "ANEOSmaterial: the interpBilinear in effect has another type layout (INTERP_ABI_VERSION) than this library; update the vendored copies together"
#endif

/*
 * Relative roundoff tolerance at the table edges: a (rho, T) within
 * ANEOS_EDGE_TOL*|bound| outside an axis is clamped onto it instead of
 * being rejected, so that a state on the edge survives the code <-> cgs
 * unit round trip. Far below any physical scale.
 */
#define ANEOS_EDGE_TOL 1e-13

#define MAT_ANEOS_ICE         52
#define MAT_ANEOS_DUNITE      54
#define MAT_ANEOS_IRON        55

#define MAT_MANEOS_QUARTZ     61
#define MAT_MANEOS_FORSTERITE 62
#define MAT_MANEOS_IRON       63
#define MAT_MANEOS_FE85SI15   64
#define MAT_MANEOS_SERPENTINE 65
#define MAT_MANEOS_WATER      66

#define ANEOS_N_MATERIAL_MAX 100

/*
 * Status codes (EOSlib convention: 0 is success, -99 a generic failure;
 * the ANEOSRead* functions return 1 for an optional file that is absent).
 */
#define ANEOS_SUCCESS 0
#define ANEOS_FAIL  -99

/*
 * Return codes of ANEOSIsInTable besides ANEOS_SUCCESS (the same values as
 * the other EOS libraries; EOSlib maps them one by one).
 */
#define ANEOS_OUTSIDE_RHOMIN -1
#define ANEOS_OUTSIDE_RHOMAX -2
#define ANEOS_OUTSIDE_TMIN   -3
#define ANEOS_OUTSIDE_TMAX   -4

/*
 * Failure contract (the EOSlib convention, shared by all EOS libraries).
 *
 * A failed lookup never exits. Every function returning a double returns
 * exactly ANEOS_FAIL_VALUE, in any unit system, when rho or T is outside the
 * closed table intervals or an inversion fails; ANEOSPCTofRhoU sets *c and
 * *T as well. The one exception is T(rho, u) at a valid rho: if u cannot be
 * inverted (below the cold curve, above the top row, not finite) T falls
 * back to TAxis[0], the lowest tabulated temperature, and all (rho, u)
 * functions then evaluate at (rho, T_min). The phase lookups return -1.
 *
 * Every failure and every fallback is counted in nLookupFail and reported
 * on stderr only when the library is compiled with -DEOSLIB_VERBOSE
 * (ANEOS_WARN). Test a returned value with ANEOSFailed(x) (true for the
 * sentinel and for NaN); EOSlib exposes the same test as EOSFailed().
 * Errors during initialisation (missing or defective table: short read,
 * non-positive rho0, axes not positive and strictly increasing, a non-finite
 * table entry, a negative or non-finite melt curve; unknown material) and a
 * missing melt curve in ANEOSTmeltofRho are fatal.
 */
#define ANEOS_FAIL_VALUE (-1e50)

static inline int ANEOSFailed(double x) {
    return !(x > -1e40);
}

/*
 * Diagnostic messages of the lookups, only with -DEOSLIB_VERBOSE. The
 * silent variant still type-checks the arguments but compiles to nothing.
 */
#ifdef EOSLIB_VERBOSE
#define ANEOS_WARN(...) fprintf(stderr, __VA_ARGS__)
#else
#define ANEOS_WARN(...) do { if (0) fprintf(stderr, __VA_ARGS__); } while (0)
#endif

typedef struct ANEOSmaterial
{
    int iMat; // Material number
    double rho0; // reference density
    int nRho; // Number of entries in the interpolation arrays in the rho dimension
    int nT; // Number of entries in the interpolation arrays in the T dimension
    char matName[256];
    char matstring[1024];

    // unit conversion
    double CodeUnitstoCGSforU;
    double CodeUnitstoCGSforP;
    double CodeUnitstoCGSforRho;
    double CodeUnitstoCGSforC;
    double CGStoCodeUnitsforU;
    double CGStoCodeUnitsforP;
    double CGStoCodeUnitsforRho;
    double CGStoCodeUnitsforC;

    /*
     * Interpolation tables. Each is one contiguous block, rho-major:
     *   xData[i*nT + j] = x(rhoAxis[i], TAxis[j])
     * and xArray[i] points at column i, so xArray[i][j] indexes the same
     * element. The inner index j (temperature) is contiguous, which is what
     * the backward T search and the bilinear cell evaluation walk.
     * The on-disk format is unchanged (T-major); the transpose happens at load.
     */
    double *rhoAxis;
    double *TAxis;
    INTERPGRID grid;   /* {nRho, nT, rhoAxis, TAxis}: the view the interpBilinear core works on */
    double *uData, *pData, *cData, *sData;
    double **uArray;
    double **pArray;
    double **cArray;
    double **sArray;

    // Optional phase table, same layout
    int *phaseData;
    int **PhaseArray;

    // Melt curve (optional)
    double *T_melt;

    // Yield parameters (optional)
    int yieldStrengthModel;
    double Y0; // Yield strength at 0 pressure
    double YM; // Yield strength at infinite pressure
    double mui; // coefficient of internal friction for intact material
    double mud; // coefficient of internal friction for damaged material
    double xi; // thermal softening parameter
    double Gamma; // Shear modulus
    int fixedGamma;

    /* Diagnostics: number of lookup failures since init (see ANEOS_WARN). */
    long nLookupFail;

} ANEOSMATERIAL;

// Initialization and finalization

ANEOSMATERIAL *ANEOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit);
ANEOSMATERIAL *ANEOSinitMaterialFromFile(int iMat, char *inputfile, double dKpcUnit, double dMsolUnit);
void ANEOSfinalizeMaterial(ANEOSMATERIAL *material);

int ANEOSReadExtendedTable(ANEOSMATERIAL *material, char *inputfile);
int ANEOSReadMeltCurve(ANEOSMATERIAL *material, char *inputfile);
int ANEOSReadStrengthParameters(ANEOSMATERIAL *material, char *inputfile);

// Access functions

double ANEOSgetRho0(ANEOSMATERIAL *material);
/* The lowest tabulated temperature (the T(rho, u) fallback), in K. */
double ANEOSTmin(ANEOSMATERIAL *material);
/* u(rho, T_min), the cold curve of the table (code units); ANEOS_FAIL_VALUE (counted) outside the table. */
double ANEOSUCold(ANEOSMATERIAL *material, double rho);
/*
 * Table status of (rho, u) (code units): ANEOS_SUCCESS inside the closed table
 * (edge tolerance ANEOS_EDGE_TOL, as the lookups), otherwise ANEOS_OUTSIDE_RHOMIN
 * (also for a NaN rho), _RHOMAX, _TMIN (u below u(rho, T_min), also NaN) or _TMAX
 * (u above u(rho, T_max)). A status query, not a failure: nothing is counted.
 */
int ANEOSIsInTable(ANEOSMATERIAL *material, double rho, double u);

double ANEOSTofPU(ANEOSMATERIAL *material, double p, double u);
double ANEOSTofRhoU(ANEOSMATERIAL *material, double rho, double u);
double ANEOSTofRhoP(ANEOSMATERIAL *material, double rho, double p);

double ANEOSUofPT(ANEOSMATERIAL *material, double p, double T);
double ANEOSUofRhoT(ANEOSMATERIAL *material, double rho, double T);
double ANEOSUofRhoP(ANEOSMATERIAL *material, double rho, double p);

double ANEOSPofUT(ANEOSMATERIAL *material, double u, double T);
double ANEOSPofRhoT(ANEOSMATERIAL *material, double rho, double T);
double ANEOSPofRhoU(ANEOSMATERIAL *material, double rho, double u);

double ANEOSRhoofUT(ANEOSMATERIAL *material, double u, double T);
double ANEOSRhoofPT(ANEOSMATERIAL *material, double p, double T);
double ANEOSRhoofPU(ANEOSMATERIAL *material, double p, double u);

double ANEOSCofRhoU(ANEOSMATERIAL *material, double rho, double u);
double ANEOSCofRhoT(ANEOSMATERIAL *material, double rho, double T);

/*
 * Fused lookup for the hydro code: P(rho,u), c(rho,u) and T(rho,u) from one
 * cell search. Returns P and fills *c and *T, all in code units. On failure
 * behaves like the separate calls: T falls back to TAxis[0]; P and c are the
 * table values at that T, or -1e50 if rho is outside the table.
 */
double ANEOSPCTofRhoU(ANEOSMATERIAL *material, double rho, double u, double *c, double *T);

// These functions are for extended EOS tables. Both return -1 if no extended
// table is loaded or (rho, T) is outside it.
int ANEOSPhaseofRhoU(ANEOSMATERIAL *material, double rho, double u);
int ANEOSPhaseofRhoT(ANEOSMATERIAL *material, double rho, double T);

// Melt temperature from the optional _tmelt table, clamped to the table's
// density range. Fatal if no melt curve is loaded; ANEOS_FAIL_VALUE (counted)
// for a rho that is not positive and finite.
double ANEOSTmeltofRho(ANEOSMATERIAL *material, double rho);

double ANEOSisentropicU(ANEOSMATERIAL *material, double rho1, double u1, double rho2);

/* Entropy: returned in cgs (erg/g/K) as tabulated, not converted to code units. */
double ANEOSSofRhoU(ANEOSMATERIAL *material, double rho, double u);
double ANEOSSofRhoT(ANEOSMATERIAL *material, double rho, double T);

double ANEOSdPdRhoofRhoU(ANEOSMATERIAL *material, double rho, double u);
double ANEOSdPdUofRhoU(ANEOSMATERIAL *material, double rho, double u);
double ANEOSdUdRhoofRhoU(ANEOSMATERIAL *material, double rho, double u);

double ANEOSdPdRhoofRhoT(ANEOSMATERIAL *material, double rho, double T);
double ANEOSdPdTofRhoT(ANEOSMATERIAL *material, double rho, double T);

int ANEOSYieldParameters(ANEOSMATERIAL *material, double *Y0, double *YM, double *mui, double *mud, double *xi);
double ANEOSGammaofRhoT(ANEOSMATERIAL *material, double rho, double T);

void ANEOSPrintMat(ANEOSMATERIAL *material, FILE *fp);
#endif