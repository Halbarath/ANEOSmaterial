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

/*
 * Read the strength parameters.
 */
int ANEOSReadStrengthParameters(ANEOSMATERIAL *material, char *inputfile)
{
	FILE *file;
    int bufferLength = 255;
    char buffer[bufferLength];
    char parameterString[bufferLength];
    char valueString[bufferLength];

    if (material == NULL) {
        fprintf(stderr, "ANEOSReadStrengthParameters: ANEOSmaterial not initialized.\n");
        assert(0);
    }

	if ((file= fopen(inputfile, "r")) == NULL)
	{
		// fprintf(stderr, "ANEOSReadStrengthParameters: Could not open file %s\n", inputfile);
        material->yieldStrengthModel = -1; // No Yield strength available
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

    while(fgets(buffer, bufferLength, file)) {
        if (buffer[0] == '#') continue;
        if (buffer[0] == '\n') continue;
        if (buffer[0] == '\r') continue;
        char *rest = buffer;
        char *token = strtok_r(rest, "=", &rest);
        strcpy(parameterString,token);
        token = strtok_r(rest, "=", &rest);
        strcpy(valueString,token);
        if (strcmp(parameterString,"yieldStrengthModel") == 0) {
            sscanf(valueString, "%d", &material->yieldStrengthModel);
        }
        else if (strcmp(parameterString,"Y0") == 0) {
            sscanf(valueString, "%lf", &material->Y0);
        }
        else if (strcmp(parameterString,"YM") == 0) {
            sscanf(valueString, "%lf", &material->YM);
        }
        else if (strcmp(parameterString,"mui") == 0) {
            sscanf(valueString, "%lf", &material->mui);
        }
        else if (strcmp(parameterString,"mud") == 0) {
            sscanf(valueString, "%lf", &material->mud);
        }
        else if (strcmp(parameterString,"xi") == 0) {
            sscanf(valueString, "%lf", &material->xi);
        }
        else if (strcmp(parameterString,"Gamma") == 0) {
            sscanf(valueString, "%lf", &material->Gamma);
            if (material->Gamma == -1.0) {
                material->fixedGamma = 0;
                material->Gamma = 0.0;
            }
        }
        else {
            fprintf(stderr,"ANEOSReadYieldParameters: Unknown parameter found: %s\n", parameterString);
            fclose(file);
            assert(0);
        }
    }

    fclose(file);

    if ((material->yieldStrengthModel < 0) || (material->Y0 < 0.0) || (material->YM < 0.0) || (material->mui < 0.0) || (material->xi < 0.0) || (material->Gamma < 0.0)) {
        fprintf(stderr, "ANEOSReadYieldParameters: Error while reading parameters\n");
        assert(0);
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
	double rho = ANEOSRhoofPU(material,p,u);
	double T = ANEOSTofRhoU(material,rho,u);
	return T;
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
 * Calculates density rho(p,u)
 */
double ANEOSRhoofPU(ANEOSMATERIAL *material, double p, double u)
{
	// This is not really tested and may not work properly
	double *Tdifference = (double *)malloc((material->nRho) * sizeof(double));
	assert(Tdifference != NULL);
	
	for (int i=0; i<material->nRho; i++) 
	{
		double T1 = backwardInterpolateTemperatureBilinear(material->rhoAxis[i],p*material->CodeUnitstoCGSforP,material->nT,material->nRho,material->rhoAxis,material->TAxis,material->pArray);
		double T2 = backwardInterpolateTemperatureBilinear(material->rhoAxis[i],u*material->CodeUnitstoCGSforU,material->nT,material->nRho,material->rhoAxis,material->TAxis,material->uArray); 
		if (T1>0 && T2>0)
		{
			Tdifference[i] = T1 - T2;
		} else {
			Tdifference[i] = -1e50;
		}
		
	}
	int k=0;
	for (int i=0; i< material->nRho-1; i++)
	{
		if (Tdifference[i]>-1e40)
		{
			if (Tdifference[i]/Tdifference[i+1]<0)
			{
				k=i;
				break;
			}
		}
	}
	
	double a = material->rhoAxis[k] * material->CGStoCodeUnitsforRho;
	double b = material->rhoAxis[k+1] * material->CGStoCodeUnitsforRho;

	double c;
	
	double rho = -1e50;
	int n = 0;
	while (n<1000)
	{
		c = (a+b)/2;
		if (ANEOSTofRhoP(material,c,p)-ANEOSTofRhoU(material,c,u) == 0 || (b-a)/2<1e-6)
		{
			rho = c;
			break;
		}
		n++;
		if (copysign(1,ANEOSTofRhoP(material,c,p)-ANEOSTofRhoU(material,c,u))== copysign(1,ANEOSTofRhoP(material,a,p)-ANEOSTofRhoU(material,a,u)))
		{
			a=c;
		} else {
			b=c;
		}
	}
	
	free(Tdifference);
	if (rho<-1e40) { material->nLookupFail++; ANEOS_WARN("ANEOSRhoofPU failed for p = %.15e, u = %.15e\n", p, u); }
	return rho;
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
 * Calculates derivative dPdRho(rho,u)
 */
double ANEOSdPdRhoofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	// Finite difference for derivative
    double h = 1e-6 * rho;
    double dPdRho = 0;

    double ucoldl = ANEOSUofRhoT(material, rho - h, material->TAxis[0]);
    double ucoldr = ANEOSUofRhoT(material, rho + h, material->TAxis[0]);

    if (ucoldl <= u && ucoldr <= u)
    {
        dPdRho = (-ANEOSPofRhoU(material, rho - h, u) + ANEOSPofRhoU(material, rho + h, u))/(2*h);
    } else if (ucoldl <= u)
    {
        dPdRho = (-ANEOSPofRhoU(material, rho - h, u) + ANEOSPofRhoU(material, rho, u))/h;
    } else if (ucoldr <= u)
    {
        dPdRho = (-ANEOSPofRhoU(material, rho, u) + ANEOSPofRhoU(material, rho + h, u))/h;
    } else
    {
        h = h*1e-4;
        dPdRho = (-ANEOSPofRhoU(material, rho - h, u) + ANEOSPofRhoU(material, rho + h, u))/(2*h);
    }

	return dPdRho;
}

/*
 * Calculates derivative dPdU(rho,u)
 */
double ANEOSdPdUofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	// Finite difference for derivative
	double h = fabs(1e-5 * u);
    double dPdU = (ANEOSPofRhoU(material, rho, u + h) - ANEOSPofRhoU(material, rho, u))/h;
	return dPdU;
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
 * Calculates derivative dPdRho(rho,T)
 */
double ANEOSdPdRhoofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
	// Finite difference for derivative
    double h = 1e-5 * rho;
    double dPdRho = 0;

    dPdRho = (ANEOSPofRhoT(material, rho + h, T) - ANEOSPofRhoT(material, rho - h, T))/(2*h);

	return dPdRho;
}

/*
 * Calculates derivative dPdU(rho,u)
 */
double ANEOSdPdTofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
	// Finite difference for derivative
	double h = 1e-5 * T;
    double dPdT = (ANEOSPofRhoT(material, rho, T + h) - ANEOSPofRhoT(material, rho, T - h))/(2*h);
	return dPdT;
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
