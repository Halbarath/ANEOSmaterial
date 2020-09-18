/*
 * ANEOS material library
 * Writes the phase data on a grid for the ANEOS material model
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "aneos.h"

int main(int argc, char *argv[])
{
	if (argc != 2) {
        fprintf(stderr,"Usage: writePressureTable <iMat>\n");
        exit(1);
    }
	
	int iMat = atoi(argv[1]);
	
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

	char matFilename[256] = "aneos.input";
	initaneos(matFilename);
	
	double T;
    double rho;
	double p;
	double u;
	double s;
    double cv;
	double dPdT;
	double dPdrho;
    double fkros;
	double c;
    double rhoL;
    double rhoH;
    double ion;
    int phase;
	
    fp = fopen("pressure.txt", "w");
    
	for (int i = 0; i<nT; i++)
	{
		if ( (i%50)==0){
			fprintf(stderr, "Iterating over T for %d of %d\n",i,nT);
		}
		for (int j = 0; j<nRho; j++)
		{
			T = TAxis[i];
			rho = rhoAxis[j];
					
			callaneos_cgs(T, rho, iMat, &p, &u, &s, &cv, &dPdT,
				&dPdrho, &fkros, &c, &phase, &rhoL, &rhoH, &ion);
            fprintf(fp,"%g %g %g\n",rho, T, p);
		}
	}
    fclose(fp);

    return 0;
}
