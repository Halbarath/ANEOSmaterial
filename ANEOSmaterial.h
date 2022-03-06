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
 
#define ANEOS_VERSION_TEXT    "1.0.0"
#define ANEOS_VERSION_MAJOR   1
#define ANEOS_VERSION_MINOR   0
#define ANEOS_VERSION_PATCH   0

#define MAT_ANEOS_ICE         52
#define MAT_ANEOS_DUNITE      54
#define MAT_ANEOS_IRON        55

#define MAT_MANEOS_DUNITE     61
#define MAT_MANEOS_FOSTERITE  62
#define MAT_MANEOS_IRON       63

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
	
	// interpolation arrays, array of array of pointers
	double *rhoAxis;
	double *TAxis;
	double **uArray;
	double **pArray;
	double **cArray;
	double **sArray;
    
    // Optional array
	int **PhaseArray;

} ANEOSMATERIAL;

// Initialization and finalization

ANEOSMATERIAL *ANEOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit);
ANEOSMATERIAL *ANEOSinitMaterialFromFile(int iMat, char *inputfile, double dKpcUnit, double dMsolUnit);
void ANEOSfinalizeMaterial(ANEOSMATERIAL *material);

int ANEOSReadExtendedTable(ANEOSMATERIAL *material, char *inputfile);

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

// These functions are for extended EOS tables.
int ANEOSPhaseofRhoT(ANEOSMATERIAL *material, double rho, double T);

double ANEOSisentropicU(ANEOSMATERIAL *material, double rho1, double u1, double rho2);

double ANEOSdPdRhoofRhoU(ANEOSMATERIAL *material, double rho, double u);
double ANEOSdPdUofRhoU(ANEOSMATERIAL *material, double rho, double u);
double ANEOSdUdRhoofRhoU(ANEOSMATERIAL *material, double rho, double u);

double ANEOSdPdRhoofRhoT(ANEOSMATERIAL *material, double rho, double T);
double ANEOSdPdTofRhoT(ANEOSMATERIAL *material, double rho, double T);

void ANEOSPrintMat(ANEOSMATERIAL *material, FILE *fp);
