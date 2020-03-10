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



ANEOSMATERIAL *ANEOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit)
{
	char inputfile[256];
	
	switch(iMat)
	{
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

	fprintf(stderr, "Initializing material...\n");
	
	int nRho;
	int nT;
	
	FILE *file;

	if ((file= fopen(inputfile, "rb")) == NULL)
	{
		fprintf(stderr, "Could not open file %s\n", inputfile);
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
	material->nRho = nRho;
	material->nT = nT;
	
	material->CodeUnitstoCGSforU = 1;
	material->CodeUnitstoCGSforP = 1;
	material->CodeUnitstoCGSforRho = 1;
	material->CodeUnitstoCGSforC = 1;
	
    const double MSOLG = 1.99e33;        /* solar mass in grams */
    const double GCGS = 6.67e-8;         /* G in cgs */
    const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */
	
	// unit of u is erg/g = cm^2/s^2
	material->CodeUnitstoCGSforU = GCGS*dMsolUnit*MSOLG/(dKpcUnit*KPCCM);
	// unit of rho is g/cm^3
	material->CodeUnitstoCGSforRho = (dMsolUnit*MSOLG)/pow(dKpcUnit*KPCCM,3.0);
	// unit of p is erg/cm^3=1/(cm*s^2)
	// die oberen multiplizieren
	material->CodeUnitstoCGSforP = material->CodeUnitstoCGSforU*material->CodeUnitstoCGSforRho;
	// unit of c is cm/s
	material->CodeUnitstoCGSforC = dKpcUnit*KPCCM*sqrt((material->CodeUnitstoCGSforRho*GCGS));

	fprintf(stderr, "Initializing arrays\n");

	double *rhoAxis = (double *)malloc(sizeof(double)*nT);
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

	fprintf(stderr, "Initializing arrays finished\n");
	
	fprintf(stderr, "Filling arrays\n");

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
	fprintf(stderr, "Arrays filled\n");
	fprintf(stderr, "Material initialized\n");
	return material;
}


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

double ANEOSTofPU(ANEOSMATERIAL *material, double p, double u)
{
	double rho = ANEOSRhoofPU(material,p,u);
	double T = ANEOSTofRhoU(material,rho,u);
	return T;
}

double ANEOSTofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	double T = backwardInterpolateTemperatureBilinear(rho*material->CodeUnitstoCGSforRho,u*material->CodeUnitstoCGSforU,material->nT,material->nRho,material->rhoAxis,material->TAxis,material->uArray); 
	if (T<-1e40)
	{
		fprintf(stderr,"ANEOSTofRhoU failed for rho = %.15e, u = %.15e\n", rho, u);
		T = 1;
	}
	return T;
}

double ANEOSTofRhoP(ANEOSMATERIAL *material, double rho, double p)
{
	double T = backwardInterpolateTemperatureBilinear(rho*material->CodeUnitstoCGSforRho,p*material->CodeUnitstoCGSforP,material->nT,material->nRho,material->rhoAxis,material->TAxis,material->pArray);
	if (T<-1e40){fprintf(stderr,"ANEOSTofRhoP failed for rho = %.15e, p = %.15e\n", rho, p);}
	return T;
}

double ANEOSUofPT(ANEOSMATERIAL *material, double p, double T)
{
	double rho = ANEOSRhoofPT(material,p,T);
	double u = ANEOSUofRhoT(material,rho,T);
	return u;
}

double ANEOSUofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
	double u = interpolateValueBilinear(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoAxis, material->TAxis, material->uArray);
	if (u<-1e40){fprintf(stderr,"ANEOSUofRhoT failed for rho = %.15e, T = %.15e\n", rho, T);}
	u /= material->CodeUnitstoCGSforU;
	return u;
}

double ANEOSUofRhoP(ANEOSMATERIAL *material, double rho, double p)
{
	double T = ANEOSTofRhoP(material,rho,p);
	double u = ANEOSUofRhoT(material,rho,T);
	return u;
}

double ANEOSPofUT(ANEOSMATERIAL *material, double u, double T)
{
	double rho = ANEOSRhoofUT(material,u,T);
	double p = ANEOSPofRhoT(material,rho,T);
	return p;
}

double ANEOSPofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
	double p = interpolateValueBilinear(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoAxis, material->TAxis, material->pArray);
	if (p<-1e40){fprintf(stderr,"ANEOSPofRhoT failed for rho = %.15e, T = %.15e\n", rho, T);}
	p /= material->CodeUnitstoCGSforP;
	return p;
}

double ANEOSPofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	double T = ANEOSTofRhoU(material,rho,u);
	double p = ANEOSPofRhoT(material,rho,T);
	return p;
}

double ANEOSRhoofUT(ANEOSMATERIAL *material, double u, double T)
{
	double rho = backwardInterpolateDensityBilinear(T,u*material->CodeUnitstoCGSforU,material->nT,material->nRho,material->rhoAxis,material->TAxis,material->uArray);
	if (rho<-1e40){fprintf(stderr,"ANEOSRhoofUT failed for u = %.15e, T = %.15e\n", u, T);}
	rho /= material->CodeUnitstoCGSforRho;
	return rho;
}

double ANEOSRhoofPT(ANEOSMATERIAL *material, double p, double T)
{
	double rho = backwardInterpolateDensityBilinear(T,p*material->CodeUnitstoCGSforP,material->nT,material->nRho,material->rhoAxis,material->TAxis,material->pArray);
	if (rho<-1e40){fprintf(stderr,"ANEOSRhoofPT failed for p = %.15e, T = %.15e\n", p, T);}
	rho /= material->CodeUnitstoCGSforRho;
	return rho;
}

double ANEOSRhoofPU(ANEOSMATERIAL *material, double p, double u)
{
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

double ANEOSCofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	double T = ANEOSTofRhoU(material,rho,u);
	double c = interpolateValueBilinear(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoAxis, material->TAxis, material->cArray);
	if (c<-1e40){fprintf(stderr,"ANEOSCofRhoU failed for rho = %.15e, u = %.15e\n", rho, u);}
	c /= material->CodeUnitstoCGSforC;
	return c;
}

double ANEOSCofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
	double c = interpolateValueBilinear(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoAxis, material->TAxis, material->cArray);
	if (c<-1e40){fprintf(stderr,"ANEOSCofRhoU failed for rho = %.15e, T = %.15e\n", rho, T);}
	c /= material->CodeUnitstoCGSforC;
	return c;
}

double ANEOSisentropicU(ANEOSMATERIAL *material, double rho1, double u1, double rho2)
{
	double oldT = ANEOSTofRhoU(material,rho1,u1);
	double oldS = interpolateValueBilinear(rho1*material->CodeUnitstoCGSforRho, oldT, material->nT, material->nRho, material->rhoAxis, material->TAxis, material->sArray);
	double newT = backwardInterpolateTemperatureBilinear(rho2*material->CodeUnitstoCGSforRho,oldS,material->nT,material->nRho,material->rhoAxis,material->TAxis,material->sArray);
	if (oldS<-1e40){fprintf(stderr,"ANEOSisentropicU oldS failed for rho1 = %.15e, u1 = %.15e, rho2 = %.15e\n", rho1, u1, rho2);}
	if (newT<-1e40){fprintf(stderr,"ANEOSisentropicU newT failed for rho1 = %.15e, u1 = %.15e, rho2 = %.15e\n", rho1, u1, rho2);}
	double newu = ANEOSUofRhoT(material,rho2,newT);
	return newu;
}

double ANEOSdPdRhoofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	// Finite difference for derivative
	double h = 0.01 * rho;
	double dPdRho = (-ANEOSPofRhoU(material, rho - h, u) + ANEOSPofRhoU(material, rho + h, u))/(2*h);
	return dPdRho;
}

double ANEOSdPdUofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	// Finite difference for derivative
	double h = 0.01 * u;
	double dPdU=(-ANEOSPofRhoU(material, rho, u - h) + ANEOSPofRhoU(material, rho, u + h))/(2*h);
	return dPdU;
}

double ANEOSdUdRhoofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
	double dUdRho = ANEOSPofRhoU(material, rho, u)/(rho*rho);
	return dUdRho;
}

double backwardInterpolateTemperatureBilinear(double rho, double z, int nT, int nRho, double* rhoAxis,
	double* TAxis, double** zArray)
	{
		
		int i=0;
		for (int k=0; k<nRho; k++)
		{
			if (rhoAxis[k]>rho)
			{
				i = k-1;
				break;
			}
		}
		
		int *indices = (int *)malloc((nT-1) * sizeof(int));
		int qq =0;
		for (int k=0; k<nT-1; k++)
		{
			// calculate minimum of the four points
			double minimum = fmin(fmin(fmin(zArray[k][i],zArray[k+1][i]),zArray[k][i+1]),zArray[k+1][i+1]);
			// calculate maximum of the four points
			double maximum = fmax(fmax(fmax(zArray[k][i],zArray[k+1][i]),zArray[k][i+1]),zArray[k+1][i+1]);

			if (z <= maximum && z >= minimum)
			{
				qq=qq+1;
				indices[k] = 1;
			} else {
				indices[k] = 0;
			}
		}
		
		double T = -1e50;
		
		for (int j=0; j<(nT-1); j++)
		{
			if (indices[j]==0)
			{
				continue;
			}
			double x=(rho-rhoAxis[i])/(rhoAxis[i+1]-rhoAxis[i]);
		
			double f00=zArray[j][i];
			double f01=zArray[j+1][i];
			double f10=zArray[j][i+1];
			double f11=zArray[j+1][i+1];
			
			double y = -(z - f10*x + f00*(x - 1))/(f10*x - f11*x - f00*(x - 1) + f01*(x - 1));

			if (y>= 0 && y<=1)
			{
				T = (TAxis[j+1]-TAxis[j])*y+TAxis[j];
				break;
			}
			
		}
		
		free(indices);
		return T;
	}
	
double backwardInterpolateDensityBilinear(double T, double z, int nT, int nRho, double* rhoAxis,
	double* TAxis, double** zArray)
	{
		
		int i=0;
		for (int k=0; k<nT; k++)
		{
			if (TAxis[k]>T)
			{
				i = k-1;
				break;
			}
		}
		
		int *indices = (int *)malloc((nRho-1) * sizeof(int));
		int qq =0;
		for (int k=0; k<nRho-1; k++)
		{
			// calculate minimum of the four points
			double minimum = fmin(fmin(fmin(zArray[i][k],zArray[i+1][k]),zArray[i][k+1]),zArray[i+1][k+1]);
			// calculate maximum of the four points
			double maximum = fmax(fmax(fmax(zArray[i][k],zArray[i+1][k]),zArray[i][k+1]),zArray[i+1][k+1]);

			if (z <= maximum && z >= minimum)
			{
				qq=qq+1;;
				indices[k] = 1;
			} else {
				indices[k] = 0;
			}
		}

		double rho = -1e50;

		for (int j=0; j<(nRho-1); j++)
		{
			if (indices[j]==0)
			{
				continue;
			}
			double y=(T-TAxis[i])/(TAxis[i+1]-TAxis[i]);
		
			double f00=zArray[i][j];
			double f01=zArray[i+1][j];
			double f10=zArray[i][j+1];
			double f11=zArray[i+1][j+1];
			
			double x = -(z - f01*y + f00*(y - 1))/(y*(f01 - f11) - (f00 - f10)*(y - 1));
			
			if (x>= 0 && x<=1)
			{
				rho = (rhoAxis[j+1]-rhoAxis[j])*x+rhoAxis[j];
				break;
			}
		}
		
		free(indices);
		return rho;
	}
	
double interpolateValueBilinear(double rho, double T, int nT, int nRho, double* rhoAxis,
	double* TAxis, double** zArray)
	{
		int i=0;
		for (int k=0; k<nRho; k++)
		{
			if (rhoAxis[k]>rho)
			{
				i = k-1;
				break;
			}
		}
		int j=0;
		for (int k=0; k<nT; k++)
		{
			if (TAxis[k]>T)
			{
				j = k-1;
				break;
			}
		}
		
		if (i < 0 || j < 0)
		{
			return -1e50;
		}
		
		double x=(rho-rhoAxis[i])/(rhoAxis[i+1]-rhoAxis[i]);
		double y=(T-TAxis[j])/(TAxis[j+1]-TAxis[j]);
		
		double f00=zArray[j][i];
		double f01=zArray[j+1][i];
		double f10=zArray[j][i+1];
		double f11=zArray[j+1][i+1];
		
		double z = -1e50;
		
		z = y*(f11*x - f01*(x - 1)) - (f10*x - f00*(x - 1))*(y - 1);
		return z;
	};