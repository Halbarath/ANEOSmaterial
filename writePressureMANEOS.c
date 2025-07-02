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
        fprintf(stderr,"Usage: writePressureMANEOS <aneos.input>\n");
        exit(1);
    }
	
	/* In M-ANEOS each input file contains one material so the iMat number is irrelevant */
	int iMat = 1;
	initaneos(argv[1]);
	
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

    double **pArray = (double **)malloc(sizeof(double*)*nT);
    for (int i=0; i<nT; i++)
	{
		pArray[i] = (double *)malloc(nRho * sizeof(double)); 
	}
	
	double T;
    double rho;
	double u;
	double s;
    double cv;
	double dPdT;
	double dPdrho;
    double fkros;
	double c;
	int phase;
    double rhoL;
    double rhoH;
    double ion;
	
	for (int i = 0; i<nT; i++)
	{
		if ( (i%50)==0){
			fprintf(stderr, "Iterating over T for %d of %d\n",i,nT);
		}
		for (int j = 0; j<nRho; j++)
		{
			T = TAxis[i];
			rho = rhoAxis[j];
					
			callaneos_cgs(T, rho, iMat, &pArray[i][j], &u, &s, &cv, &dPdT,
				&dPdrho, &fkros, &c, &phase, &rhoL, &rhoH, &ion);
		}
	}
	
	fp = fopen("pressure_maneos.txt", "w");
	fprintf(fp, "# Material: %s\n", argv[1]);

	for (int i = 0; i<nT; i++)
	{
		for (int j = 0; j<nRho; j++)
		{
			fprintf(fp," %.7E", pArray[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

    return 0;
}
