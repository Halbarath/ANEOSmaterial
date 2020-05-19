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
		case 51:
			strcpy(inputfile, "MANEOStable_water.in");
			break;
		case 52:
			strcpy(inputfile, "ANEOStable_ice.in");
			break;
		case 54:
			strcpy(inputfile, "ANEOStable_dunite.in");
			break;
		case 55:
			strcpy(inputfile, "ANEOStable_iron.in");
			break;
		default:
			assert(0);
	}

	double rho0;
	int nRho;
	int nT;
	
	FILE *file;

	if ((file= fopen(inputfile, "rb")) == NULL)
	{
		fprintf(stderr, "Could not open file %s\n", inputfile);
		assert(0);
	}
	if (fread(&rho0, sizeof(rho0), 1, file)==0)
	{
		fprintf(stderr, "Failed to read from file %s\n", inputfile);
		assert(0);
	}
	if (fread(&nRho, sizeof(nRho), 1, file)==0)
	{
		fprintf(stderr, "Failed to read from file %s\n", inputfile);
		assert(0);
	}
	if (fread(&nT, sizeof(nT), 1, file)==0)
	{
		fprintf(stderr, "Failed to read from file %s\n", inputfile);
		assert(0);
	}

	ANEOSMATERIAL *material;
	material = (ANEOSMATERIAL *) calloc(1, sizeof(ANEOSMATERIAL));
	material->iMat = iMat;
	material->rho0 = rho0;
	material->nRho = nRho;
	material->nT = nT;
	strcpy(material->matName,inputfile);
	
	material->CodeUnitstoCGSforU = 1;
	material->CodeUnitstoCGSforP = 1;
	material->CodeUnitstoCGSforRho = 1;
	material->CodeUnitstoCGSforC = 1;
	
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

	double *rhoAxis = (double *)malloc(sizeof(double)*nRho);
	double *TAxis = (double *)malloc(sizeof(double)*nT);
	double **uArray = (double **)malloc(sizeof(double*)*nT);
	double **pArray = (double **)malloc(sizeof(double*)*nT);
	double **cArray = (double **)malloc(sizeof(double*)*nT);
	double **sArray = (double **)malloc(sizeof(double*)*nT);

    for (int i=0; i<nT; i++)
	{
		uArray[i] = (double *)malloc(nRho * sizeof(double)); 
		pArray[i] = (double *)malloc(nRho * sizeof(double)); 
		cArray[i] = (double *)malloc(nRho * sizeof(double)); 
		sArray[i] = (double *)malloc(nRho * sizeof(double)); 
	}
	
	material->rhoAxis = rhoAxis;
	material->TAxis = TAxis;
	material->uArray = uArray;
	material->pArray = pArray;
	material->cArray = cArray;
	material->sArray = sArray;

	if (fread(rhoAxis, sizeof(rhoAxis[0]), nRho, file)==0)
	{
	fprintf(stderr, "Failed to read from file %s\n", inputfile);
	assert(0);
	}
	
	if (fread(TAxis, sizeof(TAxis[0]), nT, file)==0)
	{
	fprintf(stderr, "Failed to read from file %s\n", inputfile);
	assert(0);
	}
	
	for (int i=0; i< nT; i++)
	{
		if (fread(pArray[i], sizeof(pArray[i][0]), nRho, file)==0)
		{
		fprintf(stderr, "Failed to read from file %s\n", inputfile);
		assert(0);
		}
	}
	
	for (int i=0; i< nT; i++)
	{
		if (fread(uArray[i], sizeof(uArray[i][0]), nRho, file)==0)
		{
		fprintf(stderr, "Failed to read from file %s\n", inputfile);
		assert(0);
		}
	}
	
	for (int i=0; i< nT; i++)
	{
		if (fread(sArray[i], sizeof(sArray[i][0]), nRho, file)==0)
		{
		fprintf(stderr, "Failed to read from file %s\n", inputfile);
		assert(0);
		}
	}
	
	for (int i=0; i< nT; i++)
	{
		if (fread(cArray[i], sizeof(cArray[i][0]), nRho, file)==0)
		{
		fprintf(stderr, "Failed to read from file %s\n", inputfile);
		assert(0);
		}
	}
	
	fclose(file);

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
	free(material->rhoAxis);
	free(material->TAxis);
	free(material->uArray);
	free(material->pArray);
	free(material->cArray);
	free(material->sArray);
	free(material);
}

/*
 * Converts the reference density in cgs to Code Units
 */
double ANEOSgetRho0(ANEOSMATERIAL *material)
{
	return material->rho0/material->CodeUnitstoCGSforRho;
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
#ifdef EOSLIB_VERBOSE
		fprintf(stderr,"ANEOSTofRhoU failed for rho = %.15e, u = %.15e\n", rho, u);
#endif
		T = 1.0;
	}
	return T;
}

/*
 * Calculates temperature T(rho,p)
 */
double ANEOSTofRhoP(ANEOSMATERIAL *material, double rho, double p)
{
	double T = backwardInterpolateTemperatureBilinear(rho*material->CodeUnitstoCGSforRho,p*material->CodeUnitstoCGSforP,material->nT,material->nRho,material->rhoAxis,material->TAxis,material->pArray);
	if (T<-1e40){fprintf(stderr,"ANEOSTofRhoP failed for rho = %.15e, p = %.15e\n", rho, p);}
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
	if (u<-1e40){fprintf(stderr,"ANEOSUofRhoT failed for rho = %.15e, T = %.15e\n", rho, T);}
	u /= material->CodeUnitstoCGSforU;
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
	if (p<-1e40){fprintf(stderr,"ANEOSPofRhoT failed for rho = %.15e, T = %.15e\n", rho, T);}
	p /= material->CodeUnitstoCGSforP;
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
	if (rho<-1e40){fprintf(stderr,"ANEOSRhoofUT failed for u = %.15e, T = %.15e\n", u, T);}
	rho /= material->CodeUnitstoCGSforRho;
	return rho;
}

/*
 * Calculates density rho(p,T)
 */
double ANEOSRhoofPT(ANEOSMATERIAL *material, double p, double T)
{
	double rho = backwardInterpolateDensityBilinear(T,p*material->CodeUnitstoCGSforP,material->nT,material->nRho,material->rhoAxis,material->TAxis,material->pArray);
	if (rho<-1e40){fprintf(stderr,"ANEOSRhoofPT failed for p = %.15e, T = %.15e\n", p, T);}
	rho /= material->CodeUnitstoCGSforRho;
	return rho;
}

/*
 * Calculates density rho(p,u)
 */
double ANEOSRhoofPU(ANEOSMATERIAL *material, double p, double u)
{
	// This is not really tested and may not work properly
	double *Tdifference = (double *)malloc((material->nRho) * sizeof(double));
	
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
	
	double a = material->rhoAxis[k]/material->CodeUnitstoCGSforRho;
	double b = material->rhoAxis[k+1]/material->CodeUnitstoCGSforRho;

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
	if (rho<-1e40){fprintf(stderr,"ANEOSRhoofPU failed for p = %.15e, u = %.15e\n", p, u);}
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
	if (c<-1e40){fprintf(stderr,"ANEOSCofRhoT failed for rho = %.15e, T = %.15e\n", rho, T);}
	c /= material->CodeUnitstoCGSforC;
	return c;
}

/*
 * Calculates the isentropic evolution u2(rho1,u1,rho2)
 */
double ANEOSisentropicU(ANEOSMATERIAL *material, double rho1, double u1, double rho2)
{
	double oldT = ANEOSTofRhoU(material,rho1,u1);
	double oldS = interpolateValueBilinear(rho1*material->CodeUnitstoCGSforRho, oldT, material->nT, material->nRho, material->rhoAxis, material->TAxis, material->sArray);
	if (oldS<-1e40){fprintf(stderr,"ANEOSisentropicU oldS failed for rho1 = %.15e, u1 = %.15e, rho2 = %.15e\n", rho1, u1, rho2);}
	double newT = backwardInterpolateTemperatureBilinear(rho2*material->CodeUnitstoCGSforRho,oldS,material->nT,material->nRho,material->rhoAxis,material->TAxis,material->sArray);
	if (newT<-1e40){fprintf(stderr,"ANEOSisentropicU newT failed for rho1 = %.15e, u1 = %.15e, rho2 = %.15e\n", rho1, u1, rho2);}
	double newu = ANEOSUofRhoT(material,rho2,newT);
	return newu;
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

void ANEOSMatString(ANEOSMATERIAL *material, char *MatName)
{
	strcpy(MatName, material->matName);
}

void ANEOSPrintMat(ANEOSMATERIAL *material, FILE *fp)
{
	char MatName[256];
	ANEOSMatString(material, MatName);
	fprintf(fp,"# Material: %i (%s)\n", material->iMat, MatName);
	fprintf(fp,"# Reference density rho0: %g\n", material->rho0/material->CodeUnitstoCGSforRho);
	fprintf(fp,"# Table size: nRho = %d, nT = %d\n", material->nRho, material->nT);
    fprintf(fp,"# Table boundaries (cgs units): minRho = %g, maxRho = %g, minT = %g, maxT = %g\n",material->rhoAxis[0],material->rhoAxis[material->nRho-1],material->TAxis[0],material->TAxis[material->nT-1]);
}
