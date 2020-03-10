/*
 * ANEOS material library
 *
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_poly.h>
#include "ANEOSmaterial.h"



ANEOSMATERIAL *ANEOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit)
{
	if (iMat == 101)
	{
		iMat = 4;
	}
	
	fprintf(stderr, "Initializing material...\n");
	
	FILE *fp;
	char str[1000];
	int nRho;
	int nT;

	fp = fopen("ANEOStable_dunite.in", "r");
	if (fp == NULL){
        fprintf(stderr,"Could not open file ANEOStable_dunite.in");
    }
	if (fgets(str, 1000, fp) != NULL)
	{
	nRho = (int) strtol(str, (char **)NULL, 10);
	}
	if (fgets(str, 1000, fp) != NULL)
	{
	nT = (int) strtol(str, (char **)NULL, 10);
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

	double **rhoArray = (double **)malloc(sizeof(double*)*nT);
	double **TArray = (double **)malloc(sizeof(double*)*nT);
	double **uArray = (double **)malloc(sizeof(double*)*nT);
	double **dudrhoArray = (double **)malloc(sizeof(double*)*nT);
	double **dudTArray = (double **)malloc(sizeof(double*)*nT);
	double **ddudTdrhoArray = (double **)malloc(sizeof(double*)*nT);
	double **pArray = (double **)malloc(sizeof(double*)*nT);
	double **dpdrhoArray = (double **)malloc(sizeof(double*)*nT);
	double **dpdTArray = (double **)malloc(sizeof(double*)*nT);
	double **ddpdTdrhoArray = (double **)malloc(sizeof(double*)*nT);
	double **cArray = (double **)malloc(sizeof(double*)*nT);
	double **dcdrhoArray = (double **)malloc(sizeof(double*)*nT);
	double **dcdTArray = (double **)malloc(sizeof(double*)*nT);
	double **ddcdTdrhoArray = (double **)malloc(sizeof(double*)*nT);
	double **sArray = (double **)malloc(sizeof(double*)*nT);
	double **dsdrhoArray = (double **)malloc(sizeof(double*)*nT);
	double **dsdTArray = (double **)malloc(sizeof(double*)*nT);
	double **ddsdTdrhoArray = (double **)malloc(sizeof(double*)*nT);
    for (int i=0; i<nT; i++)
	{
        rhoArray[i] = (double *)malloc(nRho * sizeof(double)); 
		TArray[i] = (double *)malloc(nRho * sizeof(double)); 
		uArray[i] = (double *)malloc(nRho * sizeof(double)); 
		dudrhoArray[i] = (double *)malloc(nRho * sizeof(double)); 
		dudTArray[i] = (double *)malloc(nRho * sizeof(double)); 
		ddudTdrhoArray[i] = (double *)malloc(nRho * sizeof(double)); 
		pArray[i] = (double *)malloc(nRho * sizeof(double)); 
		dpdrhoArray[i] = (double *)malloc(nRho * sizeof(double)); 
		dpdTArray[i] = (double *)malloc(nRho * sizeof(double)); 
		ddpdTdrhoArray[i] = (double *)malloc(nRho * sizeof(double)); 
		cArray[i] = (double *)malloc(nRho * sizeof(double)); 
		dcdrhoArray[i] = (double *)malloc(nRho * sizeof(double)); 
		dcdTArray[i] = (double *)malloc(nRho * sizeof(double)); 
		ddcdTdrhoArray[i] = (double *)malloc(nRho * sizeof(double)); 
		sArray[i] = (double *)malloc(nRho * sizeof(double)); 
		dsdrhoArray[i] = (double *)malloc(nRho * sizeof(double)); 
		dsdTArray[i] = (double *)malloc(nRho * sizeof(double)); 
		ddsdTdrhoArray[i] = (double *)malloc(nRho * sizeof(double)); 
	}
	
	material->rhoArray = rhoArray;
	material->TArray = TArray;
	material->uArray = uArray;
	material->dudrhoArray = dudrhoArray;
	material->dudTArray = dudTArray;
	material->ddudTdrhoArray = ddudTdrhoArray;
	material->pArray = pArray;
	material->dpdrhoArray = dpdrhoArray;
	material->dpdTArray = dpdTArray;
	material->ddpdTdrhoArray = ddpdTdrhoArray;
	material->cArray = cArray;
	material->dcdrhoArray = dcdrhoArray;
	material->dcdTArray = dcdTArray;
	material->ddcdTdrhoArray = ddcdTdrhoArray;
	material->sArray = sArray;
	material->dsdrhoArray = dsdrhoArray;
	material->dsdTArray = dsdTArray;
	material->ddsdTdrhoArray = ddsdTdrhoArray;	
	fprintf(stderr, "Initializing arrays finished\n");
	
	fprintf(stderr, "Filling arrays\n");

	for (int i = 0; i<nT; i++)
	{
		for (int j = 0; j<nRho; j++)
		{
			if (fgets(str, 1000, fp) != NULL)
			{
				TArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				rhoArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				uArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				dudrhoArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				dudTArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				ddudTdrhoArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				pArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				dpdrhoArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				dpdTArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				ddpdTdrhoArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				cArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				dcdrhoArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				dcdTArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				ddcdTdrhoArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				sArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				dsdrhoArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				dsdTArray[i][j] = (double) strtod(str, (char **)NULL);
			}
			if (fgets(str, 1000, fp) != NULL)
			{
				ddsdTdrhoArray[i][j] = (double) strtod(str, (char **)NULL);
			}
		}
	}
	
	fclose(fp);
	
	fprintf(stderr, "Arrays filled\n");
	fprintf(stderr, "Material initialized\n");
	return material;
}


void ANEOSfinalizeMaterial(ANEOSMATERIAL *material)
{
	for (int i=0; i<material->nT; i++)
	{
        free(material->rhoArray[i]); 
		free(material->TArray[i]); 
		free(material->uArray[i]);
		free(material->dudrhoArray[i]);
		free(material->dudTArray[i]); 
		free(material->ddudTdrhoArray[i]); 
		free(material->pArray[i]); 
		free(material->dpdrhoArray[i]); 
		free(material->dpdTArray[i]); 
		free(material->ddpdTdrhoArray[i]); 
		free(material->cArray[i]); 
		free(material->dcdrhoArray[i]); 
		free(material->dcdTArray[i]); 
		free(material->ddcdTdrhoArray[i]); 
		free(material->sArray[i]); 
		free(material->dsdrhoArray[i]); 
		free(material->dsdTArray[i]); 
		free(material->ddsdTdrhoArray[i]); 
	}
	free(material->rhoArray);
	free(material->TArray);
	free(material->uArray);
	free(material->dudrhoArray);
	free(material->dudTArray);
	free(material->ddudTdrhoArray);
	free(material->pArray);
	free(material->dpdrhoArray);
	free(material->dpdTArray);
	free(material->ddpdTdrhoArray);
	free(material->cArray);
	free(material->dcdrhoArray);
	free(material->dcdTArray);
	free(material->ddcdTdrhoArray);
	free(material->sArray);
	free(material->dsdrhoArray);
	free(material->dsdTArray);
	free(material->ddsdTdrhoArray);
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
	//double T = backwardInterpolateTemperature(rho*material->CodeUnitstoCGSforRho,u*material->CodeUnitstoCGSforU,material->nT,material->nRho,material->rhoArray,material->TArray,material->uArray,material->dudrhoArray,material->dudTArray,material->ddudTdrhoArray); 
	double T = backwardInterpolateTemperatureBilinear(rho*material->CodeUnitstoCGSforRho,u*material->CodeUnitstoCGSforU,material->nT,material->nRho,material->rhoArray,material->TArray,material->uArray); 
	if (T<-1e40)
	{
		fprintf(stderr,"ANEOSTofRhoU failed for rho = %.15e, u = %.15e\n", rho, u);
		T = 1;
	}
	return T;
}

double ANEOSTofRhoP(ANEOSMATERIAL *material, double rho, double p)
{
	double T = backwardInterpolateTemperature(rho*material->CodeUnitstoCGSforRho,p*material->CodeUnitstoCGSforP,material->nT,material->nRho,material->rhoArray,material->TArray,material->pArray,material->dpdrhoArray,material->dpdTArray,material->ddpdTdrhoArray);
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
	double u = interpolateValue(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoArray, material->TArray, material->uArray, material->dudrhoArray, material->dudTArray, material->ddudTdrhoArray);
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
	//double p = interpolateValue(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoArray, material->TArray, material->pArray, material->dpdrhoArray, material->dpdTArray, material->ddpdTdrhoArray);
	double p = interpolateValueBilinear(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoArray, material->TArray, material->pArray);
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
	double rho = backwardInterpolateDensity(T,u*material->CodeUnitstoCGSforU,material->nT,material->nRho,material->rhoArray,material->TArray,material->uArray,material->dudrhoArray,material->dudTArray,material->ddudTdrhoArray);
	if (rho<-1e40){fprintf(stderr,"ANEOSRhoofUT failed for u = %.15e, T = %.15e\n", u, T);}
	rho /= material->CodeUnitstoCGSforRho;
	return rho;
}

double ANEOSRhoofPT(ANEOSMATERIAL *material, double p, double T)
{
	double rho = backwardInterpolateDensity(T,p*material->CodeUnitstoCGSforP,material->nT,material->nRho,material->rhoArray,material->TArray,material->pArray,material->dpdrhoArray,material->dpdTArray,material->ddpdTdrhoArray);
	if (rho<-1e40){fprintf(stderr,"ANEOSRhoofPT failed for p = %.15e, T = %.15e\n", p, T);}
	rho /= material->CodeUnitstoCGSforRho;
	return rho;
}

double ANEOSRhoofPU(ANEOSMATERIAL *material, double p, double u)
{
	double *Tdifference = (double *)malloc((material->nRho) * sizeof(double));
	
	for (int i=0; i<material->nRho; i++) 
	{
		double T1 = backwardInterpolateTemperature(material->rhoArray[0][i],p*material->CodeUnitstoCGSforP,material->nT,material->nRho,material->rhoArray,material->TArray,material->pArray,material->dpdrhoArray,material->dpdTArray,material->ddpdTdrhoArray);
		double T2 = backwardInterpolateTemperature(material->rhoArray[0][i],u*material->CodeUnitstoCGSforU,material->nT,material->nRho,material->rhoArray,material->TArray,material->uArray,material->dudrhoArray,material->dudTArray,material->ddudTdrhoArray); 
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
	
	double a = material->rhoArray[0][k]/material->CodeUnitstoCGSforRho;
	double b = material->rhoArray[0][k+1]/material->CodeUnitstoCGSforRho;

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
	//double c = interpolateValue(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoArray, material->TArray, material->cArray, material->dcdrhoArray, material->dcdTArray, material->ddcdTdrhoArray);
	double c = interpolateValueBilinear(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoArray, material->TArray, material->cArray);
	if (c<-1e40){fprintf(stderr,"ANEOSCofRhoU failed for rho = %.15e, u = %.15e\n", rho, u);}
	c /= material->CodeUnitstoCGSforC;
	return c;
}

double ANEOSCofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
	//double c = interpolateValue(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoArray, material->TArray, material->cArray, material->dcdrhoArray, material->dcdTArray, material->ddcdTdrhoArray);
	double c = interpolateValueBilinear(rho*material->CodeUnitstoCGSforRho, T, material->nT, material->nRho, material->rhoArray, material->TArray, material->cArray);
	if (c<-1e40){fprintf(stderr,"ANEOSCofRhoU failed for rho = %.15e, T = %.15e\n", rho, T);}
	c /= material->CodeUnitstoCGSforC;
	return c;
}

double ANEOSisentropicU(ANEOSMATERIAL *material, double rho1, double u1, double rho2)
{
	double oldT = ANEOSTofRhoU(material,rho1,u1);
	double oldS = interpolateValue(rho1*material->CodeUnitstoCGSforRho, oldT, material->nT, material->nRho, material->rhoArray, material->TArray, material->sArray, material->dsdrhoArray, material->dsdTArray, material->ddsdTdrhoArray);
	double newT = backwardInterpolateTemperature(rho2*material->CodeUnitstoCGSforRho,oldS,material->nT,material->nRho,material->rhoArray,material->TArray,material->sArray,material->dsdrhoArray,material->dsdTArray,material->ddsdTdrhoArray);
	if (oldS<-1e40){fprintf(stderr,"ANEOSisentropicU oldS failed for rho1 = %.15e, u1 = %.15e, rho2 = %.15e\n", rho1, u1, rho2);}
	if (newT<-1e40){fprintf(stderr,"ANEOSisentropicU newT failed for rho1 = %.15e, u1 = %.15e, rho2 = %.15e\n", rho1, u1, rho2);}
	double newu = ANEOSUofRhoT(material,rho2,newT);
	return newu;
}



double backwardInterpolateTemperature(double rho, double z, int nT, int nRho, double** rhoArray,
	double** TArray, double** zArray, double** dzdrhoArray, double** dzdTArray, double** ddzdTdrhoArray)
	{
		
		int i=0;
		for (int k=0; k<nRho; k++)
		{
			if (rhoArray[0][k]>rho)
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

			if (z <= 1.5*maximum && z >= 0.75*minimum)
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
			double x=(rho-rhoArray[0][i])/(rhoArray[0][i+1]-rhoArray[0][i]);
			double drhodx = rhoArray[0][i+1]-rhoArray[0][i];
			double dTdy = TArray[j+1][0]-TArray[j][0];
		
			double f00=zArray[j][i];
			double f01=zArray[j+1][i];
			double f10=zArray[j][i+1];
			double f11=zArray[j+1][i+1];
			double fy00=dzdTArray[j][i]*dTdy;
			double fy01=dzdTArray[j+1][i]*dTdy;
			double fy10=dzdTArray[j][i+1]*dTdy;
			double fy11=dzdTArray[j+1][i+1]*dTdy;
			double fx00=dzdrhoArray[j][i]*drhodx;
			double fx01=dzdrhoArray[j+1][i]*drhodx;
			double fx10=dzdrhoArray[j][i+1]*drhodx;
			double fx11=dzdrhoArray[j+1][i+1]*drhodx;
			double fxy00=ddzdTdrhoArray[j][i]*drhodx*dTdy;
			double fxy01=ddzdTdrhoArray[j+1][i]*drhodx*dTdy;
			double fxy10=ddzdTdrhoArray[j][i+1]*drhodx*dTdy;
			double fxy11=ddzdTdrhoArray[j+1][i+1]*drhodx*dTdy;

			double d=(2*f00 - 2*f10 + fx00 + fx10)*pow(x,3) + (3*f10 - 3*f00 - 2*fx00 - fx10)*pow(x,2) + fx00*x + f00-z;
			double c=(fxy00 + fxy10 + 2*fy00 - 2*fy10)*pow(x,3) + (3*fy10 - fxy10 - 3*fy00 - 2*fxy00)*pow(x,2) + fxy00*x + fy00;
			double b=(6*f01 - 6*f00 + 6*f10 - 6*f11 - 3*fx00 + 3*fx01 - 3*fx10 + 3*fx11 - 2*fxy00 - fxy01 - 2*fxy10 - fxy11 - 4*fy00 - 2*fy01 + 4*fy10 + 2*fy11)*pow(x,3) + (9*f00 - 9*f01 - 9*f10 + 9*f11 + 6*fx00 - 6*fx01 + 3*fx10 - 3*fx11 + 4*fxy00 + 2*fxy01 + 2*fxy10 + fxy11 + 6*fy00 + 3*fy01 - 6*fy10 - 3*fy11)*pow(x,2) + (3*fx01 - 3*fx00 - 2*fxy00 - fxy01)*x - 3*f00 + 3*f01 - 2*fy00 - fy01;
			double a=(4*f00 - 4*f01 - 4*f10 + 4*f11 + 2*fx00 - 2*fx01 + 2*fx10 - 2*fx11 + fxy00 + fxy01 + fxy10 + fxy11 + 2*fy00 + 2*fy01 - 2*fy10 - 2*fy11)*pow(x,3) + (6*f01 - 6*f00 + 6*f10 - 6*f11 - 4*fx00 + 4*fx01 - 2*fx10 + 2*fx11 - 2*fxy00 - 2*fxy01 - fxy10 - fxy11 - 3*fy00 - 3*fy01 + 3*fy10 + 3*fy11)*pow(x,2) + (2*fx00 - 2*fx01 + fxy00 + fxy01)*x + 2*f00 - 2*f01 + fy00 + fy01;
			
			double x0 = -5;
			double x1 = -5;
			double x2 = -5;
			
			int roots = gsl_poly_solve_cubic(b/a, c/a, d/a, &x0, &x1, &x2);

			if (x0>= -0.1 && x0<=1.1)
			{
				T = (TArray[j+1][0]-TArray[j][0])*x0+TArray[j][0];
				break;
			}
			if (x1>= -0.1 && x1<=1.1)
			{
				T = (TArray[j+1][0]-TArray[j][0])*x1+TArray[j][0];
				break;
			}
			if (x2>= -0.1 && x2<=1.1)
			{
				T = (TArray[j+1][0]-TArray[j][0])*x2+TArray[j][0];
				break;
			}
		}
		
		free(indices);
		return T;
	}
	
double backwardInterpolateTemperatureBilinear(double rho, double z, int nT, int nRho, double** rhoArray,
	double** TArray, double** zArray)
	{
		
		int i=0;
		for (int k=0; k<nRho; k++)
		{
			if (rhoArray[0][k]>rho)
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
			double x=(rho-rhoArray[0][i])/(rhoArray[0][i+1]-rhoArray[0][i]);
			double drhodx = rhoArray[0][i+1]-rhoArray[0][i];
			double dTdy = TArray[j+1][0]-TArray[j][0];
		
			double f00=zArray[j][i];
			double f01=zArray[j+1][i];
			double f10=zArray[j][i+1];
			double f11=zArray[j+1][i+1];
			
			double y = -(z - f10*x + f00*(x - 1))/(f10*x - f11*x - f00*(x - 1) + f01*(x - 1));

			if (y>= 0 && y<=1)
			{
				T = (TArray[j+1][0]-TArray[j][0])*y+TArray[j][0];
				break;
			}
			
		}
		
		free(indices);
		return T;
	}
	
double backwardInterpolateDensity(double T, double z, int nT, int nRho, double** rhoArray,
	double** TArray, double** zArray, double** dzdrhoArray, double** dzdTArray, double** ddzdTdrhoArray)
	{
		
		int i=0;
		for (int k=0; k<nT; k++)
		{
			if (TArray[k][0]>T)
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

			if (z <= 1.5*maximum && z >= 0.75*minimum)
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
			double y=(T-TArray[i][0])/(TArray[i+1][0]-TArray[i][0]);
			double drhodx = rhoArray[0][j+1]-rhoArray[0][j];
			double dTdy = TArray[i+1][0]-TArray[i][0];
		
			double f00=zArray[i][j];
			double f01=zArray[i+1][j];
			double f10=zArray[i][j+1];
			double f11=zArray[i+1][j+1];
			double fy00=dzdTArray[i][j]*dTdy;
			double fy01=dzdTArray[i+1][j]*dTdy;
			double fy10=dzdTArray[i][j+1]*dTdy;
			double fy11=dzdTArray[i+1][j+1]*dTdy;
			double fx00=dzdrhoArray[i][j]*drhodx;
			double fx01=dzdrhoArray[i+1][j]*drhodx;
			double fx10=dzdrhoArray[i][j+1]*drhodx;
			double fx11=dzdrhoArray[i+1][j+1]*drhodx;
			double fxy00=ddzdTdrhoArray[i][j]*drhodx*dTdy;
			double fxy01=ddzdTdrhoArray[i+1][j]*drhodx*dTdy;
			double fxy10=ddzdTdrhoArray[i][j+1]*drhodx*dTdy;
			double fxy11=ddzdTdrhoArray[i+1][j+1]*drhodx*dTdy;

			double d=(2*f00 - 2*f01 + fy00 + fy01)*pow(y,3) + (3*f01 - 3*f00 - 2*fy00 - fy01)*pow(y,2) + fy00*y + f00-z;
			double c=(2*fx00 - 2*fx01 + fxy00 + fxy01)*pow(y,3) + (3*fx01 - 3*fx00 - 2*fxy00 - fxy01)*pow(y,2) + fxy00*y + fx00;
			double b=(6*f01 - 6*f00 + 6*f10 - 6*f11 - 4*fx00 + 4*fx01 - 2*fx10 + 2*fx11 - 2*fxy00 - 2*fxy01 - fxy10 - fxy11 - 3*fy00 - 3*fy01 + 3*fy10 + 3*fy11)*pow(y,3) + (9*f00 - 9*f01 - 9*f10 + 9*f11 + 6*fx00 - 6*fx01 + 3*fx10 - 3*fx11 + 4*fxy00 + 2*fxy01 + 2*fxy10 + fxy11 + 6*fy00 + 3*fy01 - 6*fy10 - 3*fy11)*pow(y,2) + (3*fy10 - fxy10 - 3*fy00 - 2*fxy00)*y - 3*f00 + 3*f10 - 2*fx00 - fx10;
			double a=(4*f00 - 4*f01 - 4*f10 + 4*f11 + 2*fx00 - 2*fx01 + 2*fx10 - 2*fx11 + fxy00 + fxy01 + fxy10 + fxy11 + 2*fy00 + 2*fy01 - 2*fy10 - 2*fy11)*pow(y,3) + (6*f01 - 6*f00 + 6*f10 - 6*f11 - 3*fx00 + 3*fx01 - 3*fx10 + 3*fx11 - 2*fxy00 - fxy01 - 2*fxy10 - fxy11 - 4*fy00 - 2*fy01 + 4*fy10 + 2*fy11)*pow(y,2) + (fxy00 + fxy10 + 2*fy00 - 2*fy10)*y + 2*f00 - 2*f10 + fx00 + fx10;
			
			double y0 = -5;
			double y1 = -5;
			double y2 = -5;

			int roots = gsl_poly_solve_cubic(b/a, c/a, d/a, &y0, &y1, &y2);

			if (y0>= 0 && y0<=1)
			{
				rho = (rhoArray[0][j+1]-rhoArray[0][j])*y0+rhoArray[0][j];
				break;
			}
			if (y1>= 0 && y1<=1)
			{
				rho = (rhoArray[0][j+1]-rhoArray[0][j])*y1+rhoArray[0][j];
				break;
			}
			if (y2>= 0 && y2<=1)
			{
				rho = (rhoArray[0][j+1]-rhoArray[0][j])*y2+rhoArray[0][j];
				break;
			}
		}
		
		free(indices);
		return rho;
	}

double interpolateValue(double rho, double T, int nT, int nRho, double** rhoArray,
	double** TArray, double** zArray, double** dzdrhoArray, double** dzdTArray, double** ddzdTdrhoArray)
	{
		int i=0;
		for (int k=0; k<nRho; k++)
		{
			if (rhoArray[0][k]>rho)
			{
				i = k-1;
				break;
			}
		}
		int j=0;
		for (int k=0; k<nT; k++)
		{
			if (TArray[k][0]>T)
			{
				j = k-1;
				break;
			}
		}
		
		if (i < 0 || j < 0)
		{
			return -1e50;
		}
		
		double x=(rho-rhoArray[0][i])/(rhoArray[0][i+1]-rhoArray[0][i]);
		double y=(T-TArray[j][0])/(TArray[j+1][0]-TArray[j][0]);
		double drhodx = rhoArray[0][i+1]-rhoArray[0][i];
		double dTdy = TArray[j+1][0]-TArray[j][0];
		
		double f00=zArray[j][i];
		double f01=zArray[j+1][i];
		double f10=zArray[j][i+1];
		double f11=zArray[j+1][i+1];
		double fy00=dzdTArray[j][i]*dTdy;
		double fy01=dzdTArray[j+1][i]*dTdy;
		double fy10=dzdTArray[j][i+1]*dTdy;
		double fy11=dzdTArray[j+1][i+1]*dTdy;
		double fx00=dzdrhoArray[j][i]*drhodx;
		double fx01=dzdrhoArray[j+1][i]*drhodx;
		double fx10=dzdrhoArray[j][i+1]*drhodx;
		double fx11=dzdrhoArray[j+1][i+1]*drhodx;
		double fxy00=ddzdTdrhoArray[j][i]*drhodx*dTdy;
		double fxy01=ddzdTdrhoArray[j+1][i]*drhodx*dTdy;
		double fxy10=ddzdTdrhoArray[j][i+1]*drhodx*dTdy;
		double fxy11=ddzdTdrhoArray[j+1][i+1]*drhodx*dTdy;
		
		double z = -1e50;
		
		z = f00 - pow(y,2)*((6*f00 - 6*f01 - 6*f10 + 6*f11 + 3*fx00 - 3*fx01 + 3*fx10 - 3*fx11 + 2*fxy00 + fxy01 + 2*fxy10 + fxy11 + 4*fy00 + 2*fy01 - 4*fy10 - 2*fy11)*pow(x,3) + (9*f01 - 9*f00 + 9*f10 - 9*f11 - 6*fx00 + 6*fx01 - 3*fx10 + 3*fx11 - 4*fxy00 - 2*fxy01 - 2*fxy10 - fxy11 - 6*fy00 - 3*fy01 + 6*fy10 + 3*fy11)*pow(x,2) + (3*fx00 - 3*fx01 + 2*fxy00 + fxy01)*x + 3*f00 - 3*f01 + 2*fy00 + fy01) + y*((fxy00 + fxy10 + 2*fy00 - 2*fy10)*pow(x,3) + (3*fy10 - fxy10 - 3*fy00 - 2*fxy00)*pow(x,2) + fxy00*x + fy00) + pow(y,3)*((4*f00 - 4*f01 - 4*f10 + 4*f11 + 2*fx00 - 2*fx01 + 2*fx10 - 2*fx11 + fxy00 + fxy01 + fxy10 + fxy11 + 2*fy00 + 2*fy01 - 2*fy10 - 2*fy11)*pow(x,3) + (6*f01 - 6*f00 + 6*f10 - 6*f11 - 4*fx00 + 4*fx01 - 2*fx10 + 2*fx11 - 2*fxy00 - 2*fxy01 - fxy10 - fxy11 - 3*fy00 - 3*fy01 + 3*fy10 + 3*fy11)*pow(x,2) + (2*fx00 - 2*fx01 + fxy00 + fxy01)*x + 2*f00 - 2*f01 + fy00 + fy01) + fx00*x - pow(x,2)*(3*f00 - 3*f10 + 2*fx00 + fx10) + pow(x,3)*(2*f00 - 2*f10 + fx00 + fx10);
		return z;
	};
	
	
double interpolateValueBilinear(double rho, double T, int nT, int nRho, double** rhoArray,
	double** TArray, double** zArray)
	{
		int i=0;
		for (int k=0; k<nRho; k++)
		{
			if (rhoArray[0][k]>rho)
			{
				i = k-1;
				break;
			}
		}
		int j=0;
		for (int k=0; k<nT; k++)
		{
			if (TArray[k][0]>T)
			{
				j = k-1;
				break;
			}
		}
		
		if (i < 0 || j < 0)
		{
			return -1e50;
		}
		
		double x=(rho-rhoArray[0][i])/(rhoArray[0][i+1]-rhoArray[0][i]);
		double y=(T-TArray[j][0])/(TArray[j+1][0]-TArray[j][0]);
		double drhodx = rhoArray[0][i+1]-rhoArray[0][i];
		double dTdy = TArray[j+1][0]-TArray[j][0];
		
		double f00=zArray[j][i];
		double f01=zArray[j+1][i];
		double f10=zArray[j][i+1];
		double f11=zArray[j+1][i+1];
		
		double z = -1e50;
		
		z = y*(f11*x - f01*(x - 1)) - (f10*x - f00*(x - 1))*(y - 1);
		return z;
	};