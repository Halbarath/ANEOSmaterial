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
 * Writes the binary input file for the M-ANEOS material model.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "aneos.h"

int main(int argc, char *argv[])
{
    char inputfile[256] = "";
    char outputfile[256] = "";
    char exoutputfile[256] = "";
    char matstring[1024] = "";
    double rho0;

    int iMat = 1;
    if (argc != 4 && argc != 5) {
        fprintf(stderr,"Usage: writeMANEOStable <inputfile> <rho0> <outputfile> (\"material string\")\n");
        exit(1);
    }

    strcpy(inputfile, argv[1]);
    rho0 = atof(argv[2]);
    strcpy(outputfile, argv[3]);

    if (argc == 5) {
        strcpy(matstring, argv[4]);
    }

    fprintf(stderr, "iMat: %d rho0: %15.7E\n", iMat, rho0);
    fprintf(stderr, "inputfile: %s\n", inputfile);
    fprintf(stderr, "outputfile: %s\n", outputfile);
    fprintf(stderr, "matstring: %s\n", matstring);

    fprintf(stderr, "Read axes file\n");

    char axesFilename[256] = "axes.in";

    FILE *fp = NULL;
    char str[1000];
    int nRho;
    int nT;

    fp = fopen(axesFilename, "r");
    if (fp == NULL){
        fprintf(stderr,"Could not open file %s", axesFilename);
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

    fprintf(stderr, "Initializing material...\n");
    initaneos(inputfile);
    fprintf(stderr, "Material initializing finished...\n");

    fprintf(stderr, "Initializing arrays\n");

    /* Basic arrays that are required for generating IC and impact simulations. */
    double **rhoArray = (double **)malloc(sizeof(double*)*nT);
    double **TArray = (double **)malloc(sizeof(double*)*nT);
    double **uArray = (double **)malloc(sizeof(double*)*nT);
    double **pArray = (double **)malloc(sizeof(double*)*nT);
    double **cArray = (double **)malloc(sizeof(double*)*nT);
    double **sArray = (double **)malloc(sizeof(double*)*nT);

    /* Extended arrays that can be useful for interpreting simulation results. */
    int **PhaseArray = (int **)malloc(sizeof(int*)*nT);

    for (int i=0; i<nT; i++)
    {
        rhoArray[i] = (double *)malloc(nRho * sizeof(double)); 
        TArray[i] = (double *)malloc(nRho * sizeof(double)); 
        uArray[i] = (double *)malloc(nRho * sizeof(double)); 
        pArray[i] = (double *)malloc(nRho * sizeof(double)); 
        cArray[i] = (double *)malloc(nRho * sizeof(double)); 
        sArray[i] = (double *)malloc(nRho * sizeof(double)); 

        /* Extended arrays. */ 
        PhaseArray[i] = (int *)malloc(nRho * sizeof(int)); 
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

            callaneos_cgs(T, rho, iMat, &pArray[i][j], &uArray[i][j], &sArray[i][j], &cv, &dPdT,
                    &dPdrho, &fkros, &cArray[i][j], &PhaseArray[i][j], &rhoL, &rhoH, &ion);

            TArray[i][j] = T;
            rhoArray[i][j] = rho;
        }
    }

    fprintf(stderr, "Arrays filled\n");

    fprintf(stderr, "Correcting Pressure array\n");
    double *workP = (double *)malloc(nT * sizeof(double));
    for (int i = 0; i<nT; i++) {
        // assign line to work array
        for (int j = 0; j<nRho; j++) {
            workP[j] = pArray[i][j];
        }
        // set first value
        workP[0] = workP[0] > 1e-20 ? workP[0] : 1e-20;
        // force cummax
        double maximum = workP[0];
        for (int k = 1; k<nRho; k++) {
            if (workP[k] <= maximum) {
                workP[k] = maximum;
            }
            maximum = workP[k] > maximum ? workP[k] : maximum;
        }
        // Correct pressure
        int foundRegion = 0;
        int index1 = -1;
        int index2 = -1;
        for (int k = 1; k<nRho; k++) {
            if (!foundRegion && (fabs(workP[k] - workP[k-1]) / workP[k] <= 1e-8)) {
                index1 = k-1;
                foundRegion = 1;
            }
            if (foundRegion && (fabs(workP[k] - workP[k-1]) / workP[k] > 1e-8)) {
                index2 = k-1;
                double x1 = log10(rhoAxis[index1]);
                double x2 = log10(rhoAxis[index2+1]);
                double y1 = log10(workP[index1]);
                double y2 = log10((workP[index2] * 99.0 + workP[index2+1])/100.0);
                for (int m = index1; m<=index2; m++) {
                    workP[m] = pow(10.0,y1 + (log10(rhoAxis[m]) - x1) * (y2 - y1) / (x2 - x1));
                }
                foundRegion = 0;
            }
        }
        // assign corrected values back
        for (int j = 0; j<nRho; j++) {
            pArray[i][j] = workP[j];
        }
    }

    pArray[0][0] = pArray[0][0] > 1e-20 ? pArray[0][0] : 1e-20;
    for (int i = 0; i<nT; i++)
    {
        for (int j = 0; j<nRho; j++)
        {
            if (i > 0) {
                if (pArray[i][j] < 1.0001 * pArray[i-1][j]) {
                    pArray[i][j] = 1.0001 * pArray[i-1][j];
                }
            }
            if (j > 0) {
                if (pArray[i][j] < 1.0001 * pArray[i][j-1]) {
                    pArray[i][j] = 1.0001 * pArray[i][j-1];
                }
            }
        }
    }
    fprintf(stderr, "Pressure array corrected\n");

    fprintf(stderr, "Write file\n");

    FILE *file = fopen(outputfile, "wb");

    /* Write header. */
    fwrite(&rho0, sizeof(rho0), 1, file);
    fwrite(&nRho, sizeof(nRho), 1, file);
    fwrite(&nT, sizeof(nT), 1, file);

    fwrite(rhoAxis, sizeof(rhoAxis[0]), nRho, file);
    fwrite(TAxis, sizeof(TAxis[0]), nT, file);

    /* Arrays. */
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
    
    if (argc == 5) {
        fwrite(matstring, sizeof(matstring), 1, file);
    }

    fclose(file);
    
    /* Writing extended arrays. */
    strcpy(exoutputfile, outputfile);
    strcat(exoutputfile, "_ex");

    file = fopen(exoutputfile, "wb");

    fwrite(&nRho, sizeof(nRho), 1, file);
    fwrite(&nT, sizeof(nT), 1, file);

    fwrite(rhoAxis, sizeof(rhoAxis[0]), nRho, file);
    fwrite(TAxis, sizeof(TAxis[0]), nT, file);

    for (int i=0; i< nT; i++)
    {
        fwrite(PhaseArray[i], sizeof(PhaseArray[i][0]), nRho, file);
    }
 
    if (argc == 5) {
        fwrite(matstring, sizeof(matstring), 1, file);
    }

    fclose(file);


    fprintf(stderr, "Finished, exiting\n");
    return 0;
}
