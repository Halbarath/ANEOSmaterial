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
 * ANEOS material model
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "ANEOSmaterial.h"
#include "interpBilinear.h"

/*
 * Initializes ANEOSMATERIAL given the material number iMat
 */
ANEOSMATERIAL *ANEOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit)
{
	char inputfile[256];
	
	switch(iMat)
	{
		case MAT_ANEOS_ICE:
			strcpy(inputfile, "ANEOStable_ice.in");
			break;
		case MAT_ANEOS_DUNITE:
			strcpy(inputfile, "ANEOStable_dunite.in");
			break;
		case MAT_ANEOS_IRON:
			strcpy(inputfile, "ANEOStable_iron.in");
			break;
		case MAT_MANEOS_QUARTZ:
			strcpy(inputfile, "MANEOStable_quartz.in");
			break;
		case MAT_MANEOS_FORSTERITE:
			strcpy(inputfile, "MANEOStable_forsterite.in");
			break;
		case MAT_MANEOS_IRON:
			strcpy(inputfile, "MANEOStable_iron.in");
			break;
		case MAT_MANEOS_FE85SI15:
			strcpy(inputfile, "MANEOStable_iron_Fe85Si15.in");
			break;
		case MAT_MANEOS_SERPENTINE:
            strcpy(inputfile, "MANEOStable_serpentine.in");
			break;
		case MAT_MANEOS_WATER:
            strcpy(inputfile, "MANEOStable_water.in");
			break;
		default:
			assert(0);
	}
    
    ANEOSMATERIAL *material;
	
	material = ANEOSinitMaterialFromFile(iMat, inputfile, dKpcUnit, dMsolUnit);

	return material;
}

/*
 * Read exactly n items or abort. A short read of a binary table is always a
 * corrupt or truncated file.
 */
static void ANEOSReadOrDie(void *ptr, size_t size, size_t n, FILE *file, const char *what, const char *inputfile)
{
	size_t nRead = fread(ptr, size, n, file);
	if (nRead != n) {
		fprintf(stderr, "ANEOSmaterial: failed to read %s from file %s (got %zu of %zu items)\n",
		        what, inputfile, nRead, n);
		assert(0);
		exit(1);   /* in case NDEBUG is set */
	}
}

ANEOSMATERIAL *ANEOSinitMaterialFromFile(int iMat, char *inputfile, double dKpcUnit, double dMsolUnit)
{
	double rho0;
	int nRho;
	int nT;
	
	FILE *file;

	if ((file= fopen(inputfile, "rb")) == NULL)
	{
		fprintf(stderr, "Could not open file %s\n", inputfile);
		assert(0);
		exit(1);
	}
	ANEOSReadOrDie(&rho0, sizeof(rho0), 1, file, "rho0", inputfile);
	ANEOSReadOrDie(&nRho, sizeof(nRho), 1, file, "nRho", inputfile);
	ANEOSReadOrDie(&nT,   sizeof(nT),   1, file, "nT",   inputfile);

	if (nRho < 2 || nT < 2) {
		fprintf(stderr, "ANEOSinitMaterialFromFile: invalid table size nRho = %d, nT = %d in %s\n", nRho, nT, inputfile);
		assert(0);
		exit(1);
	}

	ANEOSMATERIAL *material;
	material = (ANEOSMATERIAL *) calloc(1, sizeof(ANEOSMATERIAL));
	assert(material != NULL);
	material->iMat = iMat;
	material->rho0 = rho0;
	material->nRho = nRho;
	material->nT = nT;
	snprintf(material->matName, sizeof(material->matName), "%s", inputfile);
	
	material->CodeUnitstoCGSforU = 1.0;
	material->CodeUnitstoCGSforP = 1.0;
	material->CodeUnitstoCGSforRho = 1.0;
	material->CodeUnitstoCGSforC = 1.0;
	
	if (dKpcUnit > 0.0 && dMsolUnit > 0.0)
	{
    const double MSOLG = 1.99e33;        /* solar mass in grams */
    const double GCGS = 6.67e-8;         /* G in cgs */
    const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */
	
	// unit of u is erg/g = cm^2/s^2
	material->CodeUnitstoCGSforU = GCGS*dMsolUnit*MSOLG/(dKpcUnit*KPCCM);
	// unit of rho is g/cm^3
	material->CodeUnitstoCGSforRho = (dMsolUnit*MSOLG)/pow(dKpcUnit*KPCCM,3.0);
	// unit of p is erg/cm^3=1/(cm*s^2)
	// combine the upper two
	material->CodeUnitstoCGSforP = material->CodeUnitstoCGSforU*material->CodeUnitstoCGSforRho;
	// unit of c is cm/s
	material->CodeUnitstoCGSforC = dKpcUnit*KPCCM*sqrt((material->CodeUnitstoCGSforRho*GCGS));
	}

    material->CGStoCodeUnitsforU = 1.0 / material->CodeUnitstoCGSforU;
    material->CGStoCodeUnitsforP = 1.0 / material->CodeUnitstoCGSforP;
    material->CGStoCodeUnitsforRho = 1.0 / material->CodeUnitstoCGSforRho;
    material->CGStoCodeUnitsforC = 1.0 / material->CodeUnitstoCGSforC;

	double *rhoAxis = (double *)malloc(sizeof(double)*nRho);
	assert(rhoAxis != NULL);
	double *TAxis = (double *)malloc(sizeof(double)*nT);
	assert(TAxis != NULL);
	double **uArray = (double **)malloc(sizeof(double*)*nT);
	assert(uArray != NULL);
	double **pArray = (double **)malloc(sizeof(double*)*nT);
	assert(pArray != NULL);
	double **cArray = (double **)malloc(sizeof(double*)*nT);
	assert(cArray != NULL);
	double **sArray = (double **)malloc(sizeof(double*)*nT);
	assert(sArray != NULL);

    for (int i=0; i<nT; i++)
	{
		uArray[i] = (double *)malloc(nRho * sizeof(double));
		assert(uArray[i] != NULL);
		pArray[i] = (double *)malloc(nRho * sizeof(double));
		assert(pArray[i] != NULL);
		cArray[i] = (double *)malloc(nRho * sizeof(double));
		assert(cArray[i] != NULL);
		sArray[i] = (double *)malloc(nRho * sizeof(double));
		assert(sArray[i] != NULL);
	}
	
	material->rhoAxis = rhoAxis;
	material->TAxis = TAxis;
	material->uArray = uArray;
	material->pArray = pArray;
	material->cArray = cArray;
	material->sArray = sArray;

	ANEOSReadOrDie(rhoAxis, sizeof(rhoAxis[0]), nRho, file, "rhoAxis", inputfile);
	ANEOSReadOrDie(TAxis,   sizeof(TAxis[0]),   nT,   file, "TAxis",   inputfile);

	for (int i=0; i< nT; i++)
		ANEOSReadOrDie(pArray[i], sizeof(pArray[i][0]), nRho, file, "pArray", inputfile);
	for (int i=0; i< nT; i++)
		ANEOSReadOrDie(uArray[i], sizeof(uArray[i][0]), nRho, file, "uArray", inputfile);
	for (int i=0; i< nT; i++)
		ANEOSReadOrDie(sArray[i], sizeof(sArray[i][0]), nRho, file, "sArray", inputfile);
	for (int i=0; i< nT; i++)
		ANEOSReadOrDie(cArray[i], sizeof(cArray[i][0]), nRho, file, "cArray", inputfile);

	/* The material string is optional: older tables end here. */
	size_t nStr = fread(material->matstring, 1, sizeof(material->matstring) - 1, file);
	if (nStr == 0) {
		snprintf(material->matstring, sizeof(material->matstring), "No material string given");
	} else {
		material->matstring[nStr] = '\0';
	}
	
	fclose(file);

	/* Optional companion files: <inputfile>_ex, _tmelt, _strength. */
	char exinputfile[sizeof(material->matName) + 16];

	snprintf(exinputfile, sizeof(exinputfile), "%s_ex", inputfile);
	if (ANEOSReadExtendedTable(material, exinputfile)) {
		//fprintf(stderr, "No extended EOS tables found.\n");
	}

	snprintf(exinputfile, sizeof(exinputfile), "%s_tmelt", inputfile);
	if (ANEOSReadMeltCurve(material, exinputfile)) {
		//fprintf(stderr, "No melting curve table found.\n");
	}

	snprintf(exinputfile, sizeof(exinputfile), "%s_strength", inputfile);
	if (ANEOSReadStrengthParameters(material, exinputfile)) {
		//fprintf(stderr, "No strength parameters found.\n");
	}

	if (material->yieldStrengthModel >= 0 && material->T_melt == NULL) {
		material->yieldStrengthModel = -1;
	}

	return material;
}

/*
 * Releases the memory of ANEOSMATERIAL
 */
void ANEOSfinalizeMaterial(ANEOSMATERIAL *material)
{
	for (int i=0; i<material->nT; i++)
	{
		free(material->uArray[i]);
		free(material->pArray[i]); 
		free(material->cArray[i]); 
		free(material->sArray[i]); 
	}

    /* Extended tables. */
    if (material->PhaseArray != NULL) {
        for (int i=0; i<material->nT; i++)
        {
            if (material->PhaseArray[i] != NULL)
                free(material->PhaseArray[i]);
        }
	free(material->PhaseArray);
    }

	if (material->T_melt != NULL) {
		free(material->T_melt);
	}

	free(material->rhoAxis);
	free(material->TAxis);
	free(material->uArray);
	free(material->pArray);
	free(material->cArray);
	free(material->sArray);
	free(material);
}

/*
 * Read an extended EOS table.
 */
int ANEOSReadExtendedTable(ANEOSMATERIAL *material, char *inputfile)
{
	FILE *file;
	int nRho;
	int nT;

    if (material == NULL) {
        fprintf(stderr, "ANEOSReadExtendedTable: ANEOSmaterial not initialized.\n");
        assert(0);
    }

	if ((file= fopen(inputfile, "rb")) == NULL)
	{
		return 1;   /* optional file, absence is not an error */
	}

    /* If the file exists we assume it should be readable. */
	ANEOSReadOrDie(&nRho, sizeof(nRho), 1, file, "nRho", inputfile);
	ANEOSReadOrDie(&nT,   sizeof(nT),   1, file, "nT",   inputfile);

    if ((material->nRho != nRho) || (material->nT != nT)) {
		fprintf(stderr, "ANEOSReadExtendedTable: Inconsistent number of grid points in %s "
		        "(table: nRho = %d, nT = %d; main: nRho = %d, nT = %d).\n",
		        inputfile, nRho, nT, material->nRho, material->nT);
        assert(0);
        exit(1);
    }

	/*
	 * Read the axes and verify that they match the main table instead of
	 * skipping them. An extended table built from a different grid would
	 * otherwise be silently misaligned.
	 */
	double *axisBuf = (double *)malloc(sizeof(double) * (nRho > nT ? nRho : nT));
	assert(axisBuf != NULL);

	ANEOSReadOrDie(axisBuf, sizeof(double), nRho, file, "rhoAxis", inputfile);
	for (int i = 0; i < nRho; i++) {
		if (axisBuf[i] != material->rhoAxis[i]) {
			fprintf(stderr, "ANEOSReadExtendedTable: rhoAxis[%d] = %.15e in %s does not match main table (%.15e).\n",
			        i, axisBuf[i], inputfile, material->rhoAxis[i]);
			assert(0);
			exit(1);
		}
	}
	ANEOSReadOrDie(axisBuf, sizeof(double), nT, file, "TAxis", inputfile);
	for (int j = 0; j < nT; j++) {
		if (axisBuf[j] != material->TAxis[j]) {
			fprintf(stderr, "ANEOSReadExtendedTable: TAxis[%d] = %.15e in %s does not match main table (%.15e).\n",
			        j, axisBuf[j], inputfile, material->TAxis[j]);
			assert(0);
			exit(1);
		}
	}
	free(axisBuf);

	int **PhaseArray = (int **)malloc(sizeof(int*) * nT);
	assert(PhaseArray != NULL);

    for (int i=0; i<nT; i++)
	{
		PhaseArray[i] = (int *)malloc(nRho * sizeof(int));
		assert(PhaseArray[i] != NULL);
	}

	for (int i=0; i< nT; i++)
		ANEOSReadOrDie(PhaseArray[i], sizeof(PhaseArray[i][0]), nRho, file, "PhaseArray", inputfile);

	/* Only publish the array once it is complete. */
	material->PhaseArray = PhaseArray;
	fclose(file);

    return 0;
}

/*
 * Read the melt curve.
 */
int ANEOSReadMeltCurve(ANEOSMATERIAL *material, char *inputfile)
{
	FILE *file;
	int nRho;

    if (material == NULL) {
        fprintf(stderr, "ANEOSReadMeltCurve: ANEOSmaterial not initialized.\n");
        assert(0);
    }

	if ((file= fopen(inputfile, "rb")) == NULL)
	{
		return 1;   /* optional file */
	}

	ANEOSReadOrDie(&nRho, sizeof(nRho), 1, file, "nRho", inputfile);

    if (material->nRho != nRho) {
		fprintf(stderr, "ANEOSReadMeltCurve: Inconsistent number of grid points in %s "
		        "(table: %d, main: %d).\n", inputfile, nRho, material->nRho);
        assert(0);
        exit(1);
    }

	double *T_melt = (double *)malloc(sizeof(double)*nRho);
	assert(T_melt != NULL);

	ANEOSReadOrDie(T_melt, sizeof(T_melt[0]), nRho, file, "T_melt", inputfile);

	fclose(file);
	material->T_melt = T_melt;

    return 0;
}

#define ANEOS_STRENGTH_LINE_MAX 256

/*
 * Read the strength parameters from a text file of the form
 *
 *   # comment
 *   yieldStrengthModel = 2
 *   Y0  = 1e8
 *   YM  = 3.5e10
 *   mui = 1.5
 *   mud = 0.6
 *   xi  = 1.2
 *   Gamma = -1      # -1: compute shear modulus from the sound speed
 *
 * Whitespace around '=' is ignored; '#' starts a comment anywhere on a line.
 * Returns 1 if the file does not exist (no strength available), 0 on success.
 * Malformed files are fatal.
 */
int ANEOSReadStrengthParameters(ANEOSMATERIAL *material, char *inputfile)
{
	FILE *file;
	char buffer[ANEOS_STRENGTH_LINE_MAX];
	char key[ANEOS_STRENGTH_LINE_MAX];
	char value[ANEOS_STRENGTH_LINE_MAX];
	int iLine = 0;

	if (material == NULL) {
		fprintf(stderr, "ANEOSReadStrengthParameters: ANEOSmaterial not initialized.\n");
		assert(0);
	}

	if ((file = fopen(inputfile, "r")) == NULL)
	{
		material->yieldStrengthModel = -1; /* No yield strength available */
		return 1;
	}

	material->yieldStrengthModel = -1;
	material->Y0 = -1.0;
	material->YM = -1.0;
	material->mui = -1.0;
	material->mud = -1.0;
	material->xi = -1.0;
	material->Gamma = -2.0;
	material->fixedGamma = 1;

	while (fgets(buffer, sizeof(buffer), file)) {
		iLine++;

		/* Reject lines longer than the buffer rather than silently splitting them. */
		if (strchr(buffer, '\n') == NULL && !feof(file)) {
			fprintf(stderr, "ANEOSReadStrengthParameters: %s line %d exceeds %d characters.\n",
			        inputfile, iLine, ANEOS_STRENGTH_LINE_MAX - 1);
			fclose(file);
			assert(0);
			exit(1);
		}

		/* Strip comments. */
		char *hash = strchr(buffer, '#');
		if (hash) *hash = '\0';

		/* Skip blank lines (only whitespace left). */
		if (sscanf(buffer, " %255s", key) != 1) continue;

		/* key = value, whitespace around '=' allowed. */
		if (sscanf(buffer, " %255[^= \t\r\n] = %255s", key, value) != 2) {
			fprintf(stderr, "ANEOSReadStrengthParameters: cannot parse %s line %d: '%s'\n",
			        inputfile, iLine, buffer);
			fclose(file);
			assert(0);
			exit(1);
		}

		char *end = NULL;
		if (strcmp(key, "yieldStrengthModel") == 0) {
			material->yieldStrengthModel = (int) strtol(value, &end, 10);
		} else {
			double v = strtod(value, &end);
			if      (strcmp(key, "Y0")  == 0) material->Y0  = v;
			else if (strcmp(key, "YM")  == 0) material->YM  = v;
			else if (strcmp(key, "mui") == 0) material->mui = v;
			else if (strcmp(key, "mud") == 0) material->mud = v;
			else if (strcmp(key, "xi")  == 0) material->xi  = v;
			else if (strcmp(key, "Gamma") == 0) {
				if (v == -1.0) {
					/* Sentinel: derive the shear modulus from the sound speed. */
					material->fixedGamma = 0;
					material->Gamma = 0.0;
				} else {
					material->fixedGamma = 1;
					material->Gamma = v;
				}
			} else {
				fprintf(stderr, "ANEOSReadStrengthParameters: unknown parameter '%s' in %s line %d\n",
				        key, inputfile, iLine);
				fclose(file);
				assert(0);
				exit(1);
			}
		}

		if (end == value || *end != '\0') {
			fprintf(stderr, "ANEOSReadStrengthParameters: invalid value '%s' for '%s' in %s line %d\n",
			        value, key, inputfile, iLine);
			fclose(file);
			assert(0);
			exit(1);
		}
	}

	fclose(file);

	if ((material->yieldStrengthModel < 0) || (material->Y0 < 0.0) || (material->YM < 0.0) ||
	    (material->mui < 0.0) || (material->xi < 0.0) || (material->Gamma < 0.0)) {
		fprintf(stderr, "ANEOSReadStrengthParameters: missing or invalid parameter in %s "
		        "(yieldStrengthModel = %d Y0 = %g YM = %g mui = %g xi = %g Gamma = %g)\n",
		        inputfile, material->yieldStrengthModel, material->Y0, material->YM,
		        material->mui, material->xi, material->Gamma);
		assert(0);
		exit(1);
	}

	return 0;
}

/*
 * Converts the reference density in cgs to Code Units
 */
double ANEOSgetRho0(ANEOSMATERIAL *material)
{
	return material->rho0 * material->CGStoCodeUnitsforRho;
}

/*
 * Calculates the temperature T(p,u)
 */
double ANEOSTofPU(ANEOSMATERIAL *material, double p, double u)
{
	double rho = ANEOSRhoofPU(material, p, u);
	if (rho < -1e40) return -1e50;
	return ANEOSTofRhoU(material, rho, u);
}

/*
 * Calculates temperature T(rho,u)
 */
double ANEOSTofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	double T = backwardInterpolateTemperatureBilinear(rho*material->CodeUnitstoCGSforRho,u*material->CodeUnitstoCGSforU,material->nT,material->nRho,material->rhoAxis,material->TAxis,material->uArray); 
	if (T<-1e40)
	{
		/* Fallback (unchanged behaviour): return the lowest tabulated temperature. */
		material->nLookupFail++;
		ANEOS_WARN("ANEOSTofRhoU failed for rho = %.15e, u = %.15e, falling back to T = %.15e\n", rho, u, material->TAxis[0]);
		T = material->TAxis[0];
	}
	return T;
}

/*
 * Calculates temperature T(rho,p)
 */
double ANEOSTofRhoP(ANEOSMATERIAL *material, double rho, double p)
{
	double T = backwardInterpolateTemperatureBilinear(rho*material->CodeUnitstoCGSforRho,p*material->CodeUnitstoCGSforP,material->nT,material->nRho,material->rhoAxis,material->TAxis,material->pArray);
	if (T<-1e40) { material->nLookupFail++; ANEOS_WARN("ANEOSTofRhoP failed for rho = %.15e, p = %.15e\n", rho, p); }
	return T;
}

/*
 * Calculates internal energy u(p,T)
 */
double ANEOSUofPT(ANEOSMATERIAL *material, double p, double T)
{
	double rho = ANEOSRhoofPT(material,p,T);
	double u = ANEOSUofRhoT(material,rho,T);
	return u;
}

/*
 * Calculates internal energy u(rho,T)
 */
double ANEOSUofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
	double u = interpolateValueBilinear(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoAxis, material->TAxis, material->uArray);
	if (u<-1e40) { material->nLookupFail++; ANEOS_WARN("ANEOSUofRhoT failed for rho = %.15e, T = %.15e\n", rho, T); }
	u *= material->CGStoCodeUnitsforU;
	return u;
}

/*
 * Calculates internal energy u(rho,p)
 */
double ANEOSUofRhoP(ANEOSMATERIAL *material, double rho, double p)
{
	double T = ANEOSTofRhoP(material,rho,p);
	double u = ANEOSUofRhoT(material,rho,T);
	return u;
}

/*
 * Calculates pressure p(u,T)
 */
double ANEOSPofUT(ANEOSMATERIAL *material, double u, double T)
{
	double rho = ANEOSRhoofUT(material,u,T);
	double p = ANEOSPofRhoT(material,rho,T);
	return p;
}

/*
 * Calculates pressure p(rho,T)
 */
double ANEOSPofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
	double p = interpolateValueBilinear(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoAxis, material->TAxis, material->pArray);
	if (p<-1e40) { material->nLookupFail++; ANEOS_WARN("ANEOSPofRhoT failed for rho = %.15e, T = %.15e\n", rho, T); }
	p *= material->CGStoCodeUnitsforP;
	return p;
}

/*
 * Calculates pressure p(rho,u)
 */
double ANEOSPofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	double T = ANEOSTofRhoU(material,rho,u);
	double p = ANEOSPofRhoT(material,rho,T);
	return p;
}

/*
 * Calculates density rho(u,T)
 */
double ANEOSRhoofUT(ANEOSMATERIAL *material, double u, double T)
{
	double rho = backwardInterpolateDensityBilinear(T,u*material->CodeUnitstoCGSforU,material->nT,material->nRho,material->rhoAxis,material->TAxis,material->uArray);
	if (rho<-1e40) { material->nLookupFail++; ANEOS_WARN("ANEOSRhoofUT failed for u = %.15e, T = %.15e\n", u, T); }
	rho *= material->CGStoCodeUnitsforRho;
	return rho;
}

/*
 * Calculates density rho(p,T)
 */
double ANEOSRhoofPT(ANEOSMATERIAL *material, double p, double T)
{
	double rho = backwardInterpolateDensityBilinear(T,p*material->CodeUnitstoCGSforP,material->nT,material->nRho,material->rhoAxis,material->TAxis,material->pArray);
	if (rho<-1e40) { material->nLookupFail++; ANEOS_WARN("ANEOSRhoofPT failed for p = %.15e, T = %.15e\n", p, T); }
	rho *= material->CGStoCodeUnitsforRho;
	return rho;
}

/*
 * g(rho) = T(rho, P) - T(rho, u) in table (cgs) coordinates.
 * Returns -1e50 if either inversion fails at this rho.
 */
static double ANEOSRhoofPUResidual(ANEOSMATERIAL *material, double rho_cgs, double p_cgs, double u_cgs)
{
	double TP = backwardInterpolateTemperatureBilinear(rho_cgs, p_cgs, material->nT, material->nRho,
	            material->rhoAxis, material->TAxis, material->pArray);
	double Tu = backwardInterpolateTemperatureBilinear(rho_cgs, u_cgs, material->nT, material->nRho,
	            material->rhoAxis, material->TAxis, material->uArray);
	if (TP < -1e40 || Tu < -1e40) return -1e50;
	return TP - Tu;
}

/*
 * Calculates density rho(p,u): the density at which the isobar P and the
 * isoenergy u intersect, i.e. T(rho, P) = T(rho, u).
 *
 * The residual is evaluated at every tabulated density; the first interval
 * (scanning from low to high density) where it changes sign is then refined
 * by bisection to a relative tolerance. If several intersections exist the
 * lowest-density one is returned, matching the previous behaviour.
 *
 * Returns rho in code units, or -1e50 if no bracketed root exists.
 */
double ANEOSRhoofPU(ANEOSMATERIAL *material, double p, double u)
{
	const double tol = 1e-10;
	const int maxIter = 100;
	double p_cgs = p * material->CodeUnitstoCGSforP;
	double u_cgs = u * material->CodeUnitstoCGSforU;
	int nRho = material->nRho;

	/* Bracket: find consecutive nodes with valid residuals of opposite sign. */
	double a = -1.0, b = -1.0, ga = 0.0, gb = 0.0;
	double gprev = -1e50;
	for (int i = 0; i < nRho - 1; i++) {        /* rhoAxis[nRho-1] is exclusive */
		double g = ANEOSRhoofPUResidual(material, material->rhoAxis[i], p_cgs, u_cgs);
		if (g > -1e40 && gprev > -1e40) {
			if (g == 0.0) return material->rhoAxis[i] * material->CGStoCodeUnitsforRho;
			if ((g < 0.0) != (gprev < 0.0)) {
				a = material->rhoAxis[i-1]; ga = gprev;
				b = material->rhoAxis[i];   gb = g;
				break;
			}
		}
		gprev = g;
	}
	if (a < 0.0) {
		material->nLookupFail++;
		ANEOS_WARN("ANEOSRhoofPU: no bracket found for p = %.15e, u = %.15e\n", p, u);
		return -1e50;
	}

	/* Bisection inside [a, b]. */
	double c = 0.5 * (a + b);
	for (int n = 0; n < maxIter && 2.0 * (b - a) / (a + b) > tol; n++) {
		c = 0.5 * (a + b);
		double gc = ANEOSRhoofPUResidual(material, c, p_cgs, u_cgs);
		if (gc < -1e40) {
			material->nLookupFail++;
			ANEOS_WARN("ANEOSRhoofPU: residual undefined at rho = %.15e for p = %.15e, u = %.15e\n",
			           c * material->CGStoCodeUnitsforRho, p, u);
			return -1e50;
		}
		if (gc == 0.0) break;
		if ((gc < 0.0) == (ga < 0.0)) { a = c; ga = gc; }
		else                          { b = c; gb = gc; }
	}

	return c * material->CGStoCodeUnitsforRho;
}

/*
 * Calculates the sound speed c(rho,u)
 */
double ANEOSCofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	double T = ANEOSTofRhoU(material,rho,u);
	double c = ANEOSCofRhoT(material,rho,T);
	return c;
}

/*
 * Calculates the sound speed c(rho,T)
 */
double ANEOSCofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
	double c = interpolateValueBilinear(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoAxis, material->TAxis, material->cArray);
	if (c<-1e40) { material->nLookupFail++; ANEOS_WARN("ANEOSCofRhoT failed for rho = %.15e, T = %.15e\n", rho, T); }
	c *= material->CGStoCodeUnitsforC;
	return c;
}

/*
 * Calculates pressure phase(rho,u)
 */
int ANEOSPhaseofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	double T = ANEOSTofRhoU(material,rho,u);
	int phase = ANEOSPhaseofRhoT(material, rho, T);
	return phase;
}

/*
 * Calculates the phase(rho,T) from the extended EOS table (nearest neighbour).
 *
 * Returns -1 if no extended table was loaded for this material or if (rho, T)
 * lies outside the table.
 */
int ANEOSPhaseofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
	if (material->PhaseArray == NULL) {
		material->nLookupFail++;
		ANEOS_WARN("ANEOSPhaseofRhoT: no extended EOS table loaded for material %d.\n", material->iMat);
		return -1;
	}

	int phase = interpolateValueNearest(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoAxis, material->TAxis, material->PhaseArray);
	if (phase < 0) {
		material->nLookupFail++;
		ANEOS_WARN("ANEOSPhaseofRhoT failed for rho = %.15e, T = %.15e\n", rho, T);
		return -1;
	}
	return phase;
}

/*
 * Calculates the isentropic evolution u2(rho1,u1,rho2)
 */
double ANEOSisentropicU(ANEOSMATERIAL *material, double rho1, double u1, double rho2)
{
	double oldT = ANEOSTofRhoU(material,rho1,u1);
	double oldS = interpolateValueBilinear(rho1*material->CodeUnitstoCGSforRho, oldT, material->nT, material->nRho, material->rhoAxis, material->TAxis, material->sArray);
	if (oldS<-1e40) { material->nLookupFail++; ANEOS_WARN("ANEOSisentropicU oldS failed for rho1 = %.15e, u1 = %.15e, rho2 = %.15e\n", rho1, u1, rho2); }
	double newT = backwardInterpolateTemperatureBilinear(rho2*material->CodeUnitstoCGSforRho,oldS,material->nT,material->nRho,material->rhoAxis,material->TAxis,material->sArray);
	if (newT<-1e40) { material->nLookupFail++; ANEOS_WARN("ANEOSisentropicU newT failed for rho1 = %.15e, u1 = %.15e, rho2 = %.15e\n", rho1, u1, rho2); }
	double newu = ANEOSUofRhoT(material,rho2,newT);
	return newu;
}

/*
 * Calculates the entropy S(rho,u)
 */
double ANEOSSofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	double T = ANEOSTofRhoU(material, rho, u);
	double S = ANEOSSofRhoT(material, rho, T);
    return S;
}

/*
 * Calculates the entropy S(rho,T).
 * Returned in cgs (erg/g/K) directly from the table; no unit conversion by design.
 */
double ANEOSSofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
	double S = interpolateValueBilinear(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoAxis, material->TAxis, material->sArray);
	if (S<-1e40) { material->nLookupFail++; ANEOS_WARN("ANEOSSofRhoT failed for rho = %.15e, T = %.15e\n", rho, T); }
    return S;
}

/*
 * Calculates the melting temperature T_melt(rho).
 *
 * Outside the tabulated density range the melt curve is held constant at the
 * nearest edge value (conservative for the strength model: T_melt increases
 * with compression, so clamping above the table underestimates it).
 *
 * Returns -1.0 if the melt curve was not loaded or rho is not a valid density.
 */
double ANEOSTmeltofRho(ANEOSMATERIAL *material, double rho)
{
	if (material->T_melt == NULL) {
		fprintf(stderr, "ANEOSTmeltofRho: Melt curve was not loaded for material %d.\n", material->iMat);
		return -1.0;
	}

	double rho_cgs = rho * material->CodeUnitstoCGSforRho;
	if (!(rho_cgs > 0.0) || !isfinite(rho_cgs)) return -1.0;

	/* Clamp to the tabulated range (upper bound is exclusive in the interpolator). */
	double rhoMin = material->rhoAxis[0];
	double rhoMax = material->rhoAxis[material->nRho-1];
	if (rho_cgs <= rhoMin) return material->T_melt[0];
	if (rho_cgs >= rhoMax) return material->T_melt[material->nRho-1];

	double T = interpolateValueLinear(rho_cgs, material->nRho, material->rhoAxis, material->T_melt);
	if (T < -1e40) {
		/* Cannot happen after the clamp, but keep the contract explicit. */
		fprintf(stderr, "ANEOSTmeltofRho: interpolation failed for rho = %.15e\n", rho);
		return -1.0;
	}
	return T;
}

/*
 * Partial derivatives of P(rho,T) and u(rho,T) at fixed T and fixed rho,
 * in cgs. Shared by the public derivative functions below.
 * Returns 0 on success, -1 if (rho, T) is outside the table.
 */
static int ANEOSPartialsCGS(ANEOSMATERIAL *material, double rho, double T,
                            double *dPdrho_T, double *dPdT_rho,
                            double *dudrho_T, double *dudT_rho)
{
	double rho_cgs = rho * material->CodeUnitstoCGSforRho;
	int iret = 0;

	iret |= interpolateGradientBilinear(rho_cgs, T, material->nT, material->nRho,
	        material->rhoAxis, material->TAxis, material->pArray, dPdrho_T, dPdT_rho);
	if (dudrho_T && dudT_rho) {
		iret |= interpolateGradientBilinear(rho_cgs, T, material->nT, material->nRho,
		        material->rhoAxis, material->TAxis, material->uArray, dudrho_T, dudT_rho);
	}
	if (iret) {
		material->nLookupFail++;
		ANEOS_WARN("ANEOSPartialsCGS failed for rho = %.15e, T = %.15e\n", rho, T);
		return -1;
	}
	return 0;
}

/*
 * dP/drho at constant u (code units).
 *
 *   dP/drho|_u = dP/drho|_T - dP/dT|_rho * (du/drho|_T) / (du/dT|_rho)
 */
double ANEOSdPdRhoofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	double dPdrho_T, dPdT_rho, dudrho_T, dudT_rho;
	double T = ANEOSTofRhoU(material, rho, u);

	if (ANEOSPartialsCGS(material, rho, T, &dPdrho_T, &dPdT_rho, &dudrho_T, &dudT_rho)) {
		return -1e50;
	}
	if (!(dudT_rho > 0.0)) {
		/* Heat capacity must be positive; a non-positive du/dT is a table artefact. */
		material->nLookupFail++;
		ANEOS_WARN("ANEOSdPdRhoofRhoU: du/dT|rho = %.15e <= 0 at rho = %.15e, T = %.15e\n", dudT_rho, rho, T);
		return -1e50;
	}

	double dPdrho_u = dPdrho_T - dPdT_rho * dudrho_T / dudT_rho;

	/* cgs: erg cm^-3 / (g cm^-3) -> code units. */
	return dPdrho_u * material->CGStoCodeUnitsforP * material->CodeUnitstoCGSforRho;
}

/*
 * dP/du at constant rho (code units).
 *
 *   dP/du|_rho = dP/dT|_rho / du/dT|_rho
 */
double ANEOSdPdUofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	double dPdrho_T, dPdT_rho, dudrho_T, dudT_rho;
	double T = ANEOSTofRhoU(material, rho, u);

	if (ANEOSPartialsCGS(material, rho, T, &dPdrho_T, &dPdT_rho, &dudrho_T, &dudT_rho)) {
		return -1e50;
	}
	if (!(dudT_rho > 0.0)) {
		material->nLookupFail++;
		ANEOS_WARN("ANEOSdPdUofRhoU: du/dT|rho = %.15e <= 0 at rho = %.15e, T = %.15e\n", dudT_rho, rho, T);
		return -1e50;
	}

	/* cgs: erg cm^-3 / (erg g^-1) -> code units. */
	return (dPdT_rho / dudT_rho) * material->CGStoCodeUnitsforP * material->CodeUnitstoCGSforU;
}

/*
 * Calculates derivative dUdRho(rho,u)
 */
double ANEOSdUdRhoofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	double dUdRho = ANEOSPofRhoU(material, rho, u)/(rho*rho);
	return dUdRho;
}

/*
 * dP/drho at constant T (code units).
 */
double ANEOSdPdRhoofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
	double dPdrho_T, dPdT_rho;

	if (ANEOSPartialsCGS(material, rho, T, &dPdrho_T, &dPdT_rho, NULL, NULL)) {
		return -1e50;
	}
	return dPdrho_T * material->CGStoCodeUnitsforP * material->CodeUnitstoCGSforRho;
}

/*
 * dP/dT at constant rho (code units; T is in K in both unit systems).
 */
double ANEOSdPdTofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
	double dPdrho_T, dPdT_rho;

	if (ANEOSPartialsCGS(material, rho, T, &dPdrho_T, &dPdT_rho, NULL, NULL)) {
		return -1e50;
	}
	return dPdT_rho * material->CGStoCodeUnitsforP;
}

/*
 * Return yield strength parameters in code units
 */
int ANEOSYieldParameters(ANEOSMATERIAL *material, double *Y0, double *YM, double *mui, double *mud, double *xi) {
    if (Y0) *Y0 = material->Y0 * material->CGStoCodeUnitsforP;
    if (YM) *YM = material->YM * material->CGStoCodeUnitsforP;
    if (mui) *mui = material->mui;
    if (mud) *mud = material->mud;
    if (xi) *xi = material->xi;
    return material->yieldStrengthModel;
}

/*
 * Return shear modulus in code units
 */
double ANEOSGammaofRhoT(ANEOSMATERIAL *material, double rho, double T) {
    if (material->fixedGamma) {
        return material->Gamma * material->CGStoCodeUnitsforP;
    } else {
        // We guess nu = 0.25 as "normal" materials have nu \in [0.2,0.3]
        // and the difference in the resulting shear modulus is less then 25%
        double nu = 0.25;
        double c = ANEOSCofRhoT(material, rho, T);
        double bulkModulus = rho * c * c;
        return bulkModulus * 3.0 * (1.0 - 2.0 * nu) / (2.0 * (1.0 + nu));
    }
}

void ANEOSPrintMat(ANEOSMATERIAL *material, FILE *fp)
{
	fprintf(fp,"# Material: %i (%s)\n", material->iMat, material->matName);
	fprintf(fp,"# Reference density rho0: %g\n", material->rho0 * material->CGStoCodeUnitsforRho);
	fprintf(fp,"# Table size: nRho = %d, nT = %d\n", material->nRho, material->nT);
    fprintf(fp,"# Table boundaries (cgs units): minRho = %g, maxRho = %g, minT = %g, maxT = %g\n",material->rhoAxis[0],material->rhoAxis[material->nRho-1],material->TAxis[0],material->TAxis[material->nT-1]);
    fprintf(fp,"# %s\n",material->matstring);
}
