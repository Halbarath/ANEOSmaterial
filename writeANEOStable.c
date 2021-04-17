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
 * Writes the binary input file for the ANEOS material model
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "aneos.h"
#include <string.h>

int main(int argc, char *argv[])
{
    if (argc != 3 && argc != 4) {
        fprintf(stderr,"Usage: writeANEOStable <iMat> <outputfile> (\"material string\")\n");
        exit(1);
    }

    int iMat = atoi(argv[1]);
    char outputfile[256] = "";
    strcpy(outputfile, argv[2]);
    char matstring[1024] = "";
    if (argc == 4) {
        strcpy(matstring, argv[3]);
    }

    double rho0;

    switch(iMat)
    {
        case 1:
            rho0 = 1.1;
            break;
        case 2:
            rho0 = 1.11;
            break;
        case 4:
            rho0 = 3.32;
            break;
        case 5:
            rho0 = 7.85;
            break;
        default:
            fprintf(stderr, "Unknown material id: %i\n", iMat);
            exit(1);
    }


    fprintf(stderr, "iMat: %d\n",iMat);
    fprintf(stderr, "outputfile: %s\n",outputfile);

    fprintf(stderr, "Read axes file\n");

    char axesFilename[256] = "axes.in";

    FILE *fp = NULL;
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

    fprintf(stderr, "Initializing material...\n");
    char matFilename[256] = "aneos.input";
    // Uncomment below to use M-ANEOS
    //char matFilename[256] = "maneos.in";
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

            callaneos_cgs(T, rho, iMat, &pArray[i][j], &uArray[i][j], &sArray[i][j], &cv, &dPdT,
                    &dPdrho, &fkros, &cArray[i][j], &iPhase, &rhoL, &rhoH, &ion);

            TArray[i][j] = T;
            rhoArray[i][j] = rho;
        }
    }
    fprintf(stderr, "Arrays filled\n");

    fprintf(stderr, "Write file\n");

    FILE *file = fopen(outputfile, "wb");

    fwrite(&rho0, sizeof(rho0), 1, file);
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
    
    if (argc == 4) {
        fwrite(matstring, sizeof(matstring), 1, file);
    }

    fclose(file);

    fprintf(stderr, "Finished, exiting\n");
    return 0;
}
