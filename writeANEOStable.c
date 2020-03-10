/*
 * ANEOS material library
 *
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "aneos.h"

int main(int argc, char *argv[])
{
	int iMat = 4;
	
	fprintf(stderr, "Initializing material...\n");
	char matFilename[256] = "aneos.input";
	char axesFilename[256] = "axes.in";
	
	FILE *fp;
	char str[1000];
	int nRho;
	int nT;

	fp = fopen(axesFilename, "r");
	if (fp == NULL){
        fprintf(stderr,"Could not open file %s",axesFilename);
    }
	if (fgets(str, 1000, fp) != NULL)
	{
	nRho = (int) strtol(str, (char **)NULL, 10);
	}
	if (fgets(str, 1000, fp) != NULL)
	{
	nT = (int) strtol(str, (char **)NULL, 10);
	}
	
	double *rhoAxis = (double *)malloc(nRho * sizeof(double));
	double *TAxis = (double *)malloc(nT * sizeof(double));
	
	for (int i=0; i<nRho; i++)
	{
		if (fgets(str, 1000, fp) != NULL)
		{
			rhoAxis[i] = (double) strtod(str, (char **)NULL);
		}
	}

	for (int i=0; i<nT; i++)
	{
		if (fgets(str, 1000, fp) != NULL)
		{
			TAxis[i] = (double) strtod(str, (char **)NULL);
		}
	}
	
	fclose(fp);
	
	initaneos(matFilename);
	fprintf(stderr, "Material initializing finished...\n");
	
	fprintf(stderr, "Initializing arrays\n");

	double **rhoArray = (double **)malloc(sizeof(double*)*nT);
	double **TArray = (double **)malloc(sizeof(double*)*nT);
	double **uArray = (double **)malloc(sizeof(double*)*nT);
	double **pArray = (double **)malloc(sizeof(double*)*nT);
	double **cArray = (double **)malloc(sizeof(double*)*nT);
	double **sArray = (double **)malloc(sizeof(double*)*nT);
    for (int i=0; i<nT; i++)
	{
        rhoArray[i] = (double *)malloc(nRho * sizeof(double)); 
		TArray[i] = (double *)malloc(nRho * sizeof(double)); 
		uArray[i] = (double *)malloc(nRho * sizeof(double)); 
		pArray[i] = (double *)malloc(nRho * sizeof(double)); 
		cArray[i] = (double *)malloc(nRho * sizeof(double)); 
		sArray[i] = (double *)malloc(nRho * sizeof(double)); 
	}
	
	
	double T;
    double rho;
    double cv;
	double dPdT;
	double dPdrho;
    double fkros;
    int iPhase;
    double rhoL;
    double rhoH;
    double ion;
	
	fprintf(stderr, "Filling arrays\n");
	
	for (int i = 0; i<nT; i++)
	{
		if ( (i%50)==0){
			fprintf(stderr, "Iterating over T for %d of %d\n",i,nT);
		}
		for (int j = 0; j<nRho; j++)
		{
			T = TAxis[i];
			rho = rhoAxis[j];
			
			double uStencil[25];
			double pStencil[25];
			double cStencil[25];
			double sStencil[25];
						
			callaneos_cgs(T, rho, iMat, &pArray[i][j], &uArray[i][j], &sArray[i][j], &cv, &dPdT, &dPdrho, &fkros, &cArray[i][j], &iPhase, &rhoL, &rhoH, &ion);
			
			TArray[i][j] = T;
			rhoArray[i][j] = rho;
		}
	}
	fprintf(stderr, "Arrays filled\n");
	fprintf(stderr, "Material initialized\n");
	fprintf(stderr, "Write file\n");
	
	FILE *file = fopen("ANEOStable_dunite.in", "wb");
	
	fwrite(&nRho, sizeof(nRho), 1, file);
	fwrite(&nT, sizeof(nT), 1, file);
	
	fwrite(rhoAxis, sizeof(rhoAxis[0]), nRho, file);
	fwrite(TAxis, sizeof(TAxis[0]), nT, file);
	
	for (int i=0; i< nT; i++)
	{
		fwrite(pArray[i], sizeof(pArray[i][0]), nRho, file);
	}
	
	for (int i=0; i< nT; i++)
	{
		fwrite(uArray[i], sizeof(uArray[i][0]), nRho, file);
	}
	
	for (int i=0; i< nT; i++)
	{
		fwrite(sArray[i], sizeof(sArray[i][0]), nRho, file);
	}
	
	for (int i=0; i< nT; i++)
	{
		fwrite(cArray[i], sizeof(cArray[i][0]), nRho, file);
	}
	
	fclose(file);
	
    return 0;
}
