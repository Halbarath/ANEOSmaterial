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
	
	double dudrho;
	double dudT;
	double ddudTdrho;
	double dpdrho;
	double dpdT;
	double ddpdTdrho;
	double dcdrho;
	double dcdT;
	double ddcdTdrho;
	double dsdrho;
	double dsdT;
	double ddsdTdrho;
	
	double delta = 2;
	
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
			
			double h = 0;
			double k = 0;
			
			if (i == 0)
			{
				h = fmin((TAxis[i+1]-TAxis[i])/delta,(TAxis[i]-1e-8)/2);
			} else if (i == (nT-1)){
				h = (TAxis[i]-TAxis[i-1])/delta;
			} else {
				h = fmin(TAxis[i+1]-TAxis[i],TAxis[i]-TAxis[i-1])/delta;
			}
			
			if (j == 0)
			{
				k = fmin((rhoAxis[j+1]-rhoAxis[j])/delta,(rhoAxis[j]-1e-8)/2);
			} else if (j == (nRho-1)){
				k = (rhoAxis[j]-rhoAxis[j-1])/delta;
			} else {
				k = fmin(rhoAxis[j+1]-rhoAxis[j],rhoAxis[j]-rhoAxis[j-1])/delta;
			}
			
			callaneos_cgs(T, rho, iMat, &pStencil[0], &uStencil[0], &sStencil[0], &cv, &dPdT, &dPdrho, &fkros, &cStencil[0], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T, rho+k, iMat, &pStencil[1], &uStencil[1], &sStencil[1], &cv, &dPdT, &dPdrho, &fkros, &cStencil[1], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T+h, rho+k, iMat, &pStencil[2], &uStencil[2], &sStencil[2], &cv, &dPdT, &dPdrho, &fkros, &cStencil[2], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T+h, rho, iMat, &pStencil[3], &uStencil[3], &sStencil[3], &cv, &dPdT, &dPdrho, &fkros, &cStencil[3], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T+h, rho-k, iMat, &pStencil[4], &uStencil[4], &sStencil[4], &cv, &dPdT, &dPdrho, &fkros, &cStencil[4], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T, rho-k, iMat, &pStencil[5], &uStencil[5], &sStencil[5], &cv, &dPdT, &dPdrho, &fkros, &cStencil[5], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T-h, rho-k, iMat, &pStencil[6], &uStencil[6], &sStencil[6], &cv, &dPdT, &dPdrho, &fkros, &cStencil[6], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T-h, rho, iMat, &pStencil[7], &uStencil[7], &sStencil[7], &cv, &dPdT, &dPdrho, &fkros, &cStencil[7], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T-h, rho+k, iMat, &pStencil[8], &uStencil[8], &sStencil[8], &cv, &dPdT, &dPdrho, &fkros, &cStencil[8], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T, rho+2*k, iMat, &pStencil[9], &uStencil[9], &sStencil[9], &cv, &dPdT, &dPdrho, &fkros, &cStencil[9], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T+h, rho+2*k, iMat, &pStencil[10], &uStencil[10], &sStencil[10], &cv, &dPdT, &dPdrho, &fkros, &cStencil[10], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T+2*h, rho+2*k, iMat, &pStencil[11], &uStencil[11], &sStencil[11], &cv, &dPdT, &dPdrho, &fkros, &cStencil[11], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T+2*h, rho+k, iMat, &pStencil[12], &uStencil[12], &sStencil[12], &cv, &dPdT, &dPdrho, &fkros, &cStencil[12], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T+2*h, rho, iMat, &pStencil[13], &uStencil[13], &sStencil[13], &cv, &dPdT, &dPdrho, &fkros, &cStencil[13], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T+2*h, rho-k, iMat, &pStencil[14], &uStencil[14], &sStencil[14], &cv, &dPdT, &dPdrho, &fkros, &cStencil[14], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T+2*h, rho-2*k, iMat, &pStencil[15], &uStencil[15], &sStencil[15], &cv, &dPdT, &dPdrho, &fkros, &cStencil[15], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T+h, rho-2*k, iMat, &pStencil[16], &uStencil[16], &sStencil[16], &cv, &dPdT, &dPdrho, &fkros, &cStencil[16], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T, rho-2*k, iMat, &pStencil[17], &uStencil[17], &sStencil[17], &cv, &dPdT, &dPdrho, &fkros, &cStencil[17], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T-h, rho-2*k, iMat, &pStencil[18], &uStencil[18], &sStencil[18], &cv, &dPdT, &dPdrho, &fkros, &cStencil[18], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T-2*h, rho-2*k, iMat, &pStencil[19], &uStencil[19], &sStencil[19], &cv, &dPdT, &dPdrho, &fkros, &cStencil[19], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T-2*h, rho-k, iMat, &pStencil[20], &uStencil[20], &sStencil[20], &cv, &dPdT, &dPdrho, &fkros, &cStencil[20], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T-2*h, rho, iMat, &pStencil[21], &uStencil[21], &sStencil[21], &cv, &dPdT, &dPdrho, &fkros, &cStencil[21], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T-2*h, rho+k, iMat, &pStencil[22], &uStencil[22], &sStencil[22], &cv, &dPdT, &dPdrho, &fkros, &cStencil[22], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T-2*h, rho+2*k, iMat, &pStencil[23], &uStencil[23], &sStencil[23], &cv, &dPdT, &dPdrho, &fkros, &cStencil[23], &iPhase, &rhoL, &rhoH, &ion);
			callaneos_cgs(T-h, rho+2*k, iMat, &pStencil[24], &uStencil[24], &sStencil[24], &cv, &dPdT, &dPdrho, &fkros, &cStencil[24], &iPhase, &rhoL, &rhoH, &ion);			
			
			
			/*// second order derivatives
			dudT = (-uStencil[7]+uStencil[3])/(2*h);
			dudrho = (-uStencil[5]+uStencil[1])/(2*k);
			ddudTdrho = (uStencil[2]-uStencil[4]-uStencil[8]+uStencil[6])/(4*h*k);
			
			dpdT = (-pStencil[7]+pStencil[3])/(2*h);
			dpdrho = (-pStencil[5]+pStencil[1])/(2*k);
			ddpdTdrho = (pStencil[2]-pStencil[4]-pStencil[8]+pStencil[6])/(4*h*k);
			
			dcdT = (-cStencil[7]+cStencil[3])/(2*h);
			dcdrho = (-cStencil[5]+cStencil[1])/(2*k);
			ddcdTdrho = (cStencil[2]-cStencil[4]-cStencil[8]+cStencil[6])/(4*h*k);
			
			dsdT = (-sStencil[7]+sStencil[3])/(2*h);
			dsdrho = (-sStencil[5]+sStencil[1])/(2*k);
			ddsdTdrho = (sStencil[2]-sStencil[4]-sStencil[8]+sStencil[6])/(4*h*k);*/
					
			// third order derivatives
			dudT = (uStencil[21]-8*uStencil[7]+8*uStencil[3]-uStencil[13])/(12*h);
			dudrho = (uStencil[17]-8*uStencil[5]+8*uStencil[1]-uStencil[9])/(12*k);
			ddudTdrho = (64*uStencil[2] - 64*uStencil[4] + 64*uStencil[6] - 64*uStencil[8] - 8*uStencil[10] + uStencil[11] - 8*uStencil[12] + 8*uStencil[14] - uStencil[15] + 8*uStencil[16] - 8*uStencil[18] + uStencil[19] - 8*uStencil[20] + 8*uStencil[22] - uStencil[23] + 8*uStencil[24])/(144*h*k);
			
			dpdT = (pStencil[21]-8*pStencil[7]+8*pStencil[3]-pStencil[13])/(12*h);
			dpdrho = (pStencil[17]-8*pStencil[5]+8*pStencil[1]-pStencil[9])/(12*k);
			ddpdTdrho = (64*pStencil[2] - 64*pStencil[4] + 64*pStencil[6] - 64*pStencil[8] - 8*pStencil[10] + pStencil[11] - 8*pStencil[12] + 8*pStencil[14] - pStencil[15] + 8*pStencil[16] - 8*pStencil[18] + pStencil[19] - 8*pStencil[20] + 8*pStencil[22] - pStencil[23] + 8*pStencil[24])/(144*h*k);
			
			dcdT = (cStencil[21]-8*cStencil[7]+8*cStencil[3]-cStencil[13])/(12*h);
			dcdrho = (cStencil[17]-8*cStencil[5]+8*cStencil[1]-cStencil[9])/(12*k);
			ddcdTdrho = (64*cStencil[2] - 64*cStencil[4] + 64*cStencil[6] - 64*cStencil[8] - 8*cStencil[10] + cStencil[11] - 8*cStencil[12] + 8*cStencil[14] - cStencil[15] + 8*cStencil[16] - 8*cStencil[18] + cStencil[19] - 8*cStencil[20] + 8*cStencil[22] - cStencil[23] + 8*cStencil[24])/(144*h*k);
			
			dsdT = (sStencil[21]-8*sStencil[7]+8*sStencil[3]-sStencil[13])/(12*h);
			dsdrho = (sStencil[17]-8*sStencil[5]+8*sStencil[1]-sStencil[9])/(12*k);
			ddsdTdrho = (64*sStencil[2] - 64*sStencil[4] + 64*sStencil[6] - 64*sStencil[8] - 8*sStencil[10] + sStencil[11] - 8*sStencil[12] + 8*sStencil[14] - sStencil[15] + 8*sStencil[16] - 8*sStencil[18] + sStencil[19] - 8*sStencil[20] + 8*sStencil[22] - sStencil[23] + 8*sStencil[24])/(144*h*k);
			
			TArray[i][j] = T;
			rhoArray[i][j] = rho;
			uArray[i][j] = uStencil[0];
			dudrhoArray[i][j] = dudrho;
			dudTArray[i][j] = dudT;
			ddudTdrhoArray[i][j] = ddudTdrho;
			pArray[i][j] = pStencil[0];
			dpdrhoArray[i][j] = dpdrho;
			dpdTArray[i][j] = dpdT;
			ddpdTdrhoArray[i][j] = ddpdTdrho;
			cArray[i][j] = cStencil[0];
			dcdrhoArray[i][j] = dcdrho;
			dcdTArray[i][j] = dcdT;
			ddcdTdrhoArray[i][j] = ddcdTdrho;
			sArray[i][j] = sStencil[0];
			dsdrhoArray[i][j] = dsdrho;
			dsdTArray[i][j] = dsdT;
			ddsdTdrhoArray[i][j] = ddsdTdrho;
		}
	}
	fprintf(stderr, "Arrays filled\n");
	fprintf(stderr, "Material initialized\n");
	fprintf(stderr, "Write file\n");
	
	fp = fopen("ANEOStable_dunite.in", "w");
	
	fprintf(fp," %d\n",nRho);
	fprintf(fp," %d\n",nT);
	for (int i = 0; i<nT; i++)
	{
		for (int j = 0; j<nRho; j++)
		{
			fprintf(fp," %.15e\n", TArray[i][j]);
			fprintf(fp," %.15e\n", rhoArray[i][j]);
			fprintf(fp," %.15e\n", uArray[i][j]);
			fprintf(fp," %.15e\n", dudrhoArray[i][j]);
			fprintf(fp," %.15e\n", dudTArray[i][j]);
			fprintf(fp," %.15e\n", ddudTdrhoArray[i][j]);
			fprintf(fp," %.15e\n", pArray[i][j]);
			fprintf(fp," %.15e\n", dpdrhoArray[i][j]);
			fprintf(fp," %.15e\n", dpdTArray[i][j]);
			fprintf(fp," %.15e\n", ddpdTdrhoArray[i][j]);
			fprintf(fp," %.15e\n", cArray[i][j]);
			fprintf(fp," %.15e\n", dcdrhoArray[i][j]);
			fprintf(fp," %.15e\n", dcdTArray[i][j]);
			fprintf(fp," %.15e\n", ddcdTdrhoArray[i][j]);
			fprintf(fp," %.15e\n", sArray[i][j]);
			fprintf(fp," %.15e\n", dsdrhoArray[i][j]);
			fprintf(fp," %.15e\n", dsdTArray[i][j]);
			fprintf(fp," %.15e\n", ddsdTdrhoArray[i][j]);
		}
	}
	fclose(fp);	
	
	
    return 0;
}
