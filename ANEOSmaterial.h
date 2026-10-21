/*  This file is part of ANEOSmaterial.
 *  Copyright (c) 2020-2021 Thomas Meier & Christian Reinhardt
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

#define ANEOS_VERSION_TEXT    "1.0.2"
#define ANEOS_VERSION_MAJOR   1
#define ANEOS_VERSION_MINOR   0
#define ANEOS_VERSION_PATCH   2

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
 * Diagnostics. Lookup failures (rho or T outside the table, failed inversion)
 * are counted here and reported to stderr only if EOSLIB_VERBOSE is defined.
 * Return values on failure are unchanged: ANEOSTofRhoU falls back to TAxis[0];
 * all other functions return -1e50 (scaled to code units where applicable).
 */
#ifdef EOSLIB_VERBOSE
#define ANEOS_WARN(...) fprintf(stderr, __VA_ARGS__)
#else
#define ANEOS_WARN(...) ((void)0)
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
	
	// interpolation arrays, array of array of pointers
	double *rhoAxis;
	double *TAxis;
	double **uArray;
	double **pArray;
	double **cArray;
	double **sArray;
    
    // Optional array
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
 * table values at that T, or -1e50 (scaled) if rho is outside the table.
 */
double ANEOSPCTofRhoU(ANEOSMATERIAL *material, double rho, double u, double *c, double *T);

// These functions are for extended EOS tables. Both return -1 if no extended
// table is loaded or (rho, T) is outside it.
int ANEOSPhaseofRhoU(ANEOSMATERIAL *material, double rho, double u);
int ANEOSPhaseofRhoT(ANEOSMATERIAL *material, double rho, double T);

// Melt temperature from the optional _tmelt table, clamped to the table's
// density range. Returns -1.0 if no melt curve is loaded or rho is invalid.
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