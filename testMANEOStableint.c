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
 *
 * Read an EOS table generated with writeMANEOStable and verify that the data
 * ist consistent with M-ANEOS direct calls.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "ANEOSmaterial.h"
#include "aneos.h"

int PrintArrayDouble(double **Array, int nRho, int nT, FILE *fp) {
    if (Array == NULL)
        return 1;

    if (fp == NULL)
        return 1;
    
    if ((nRho < 1) || (nT < 1))
        return 1;

    for (int i=0; i<nT; i++) {
        for (int j=1; j<nRho; j++) {
            fprintf(fp, "%15.7E", Array[i][j]);
        }
        fprintf(fp, "\n");
    }

    return 0;
}

int PrintArrayInt(int **Array, int nRho, int nT, FILE *fp) {
    if (Array == NULL)
        return 1;

    if (fp == NULL)
        return 1;
    
    if ((nRho < 1) || (nT < 1))
        return 1;
    
    for (int i=0; i<nT; i++) {
        for (int j=1; j<nRho; j++) {
            fprintf(fp, "%4i", Array[i][j]);
        }
        fprintf(fp, "\n");
    }

    return 0;
}

int main(int argc, char *argv[])
{
    ANEOSMATERIAL *Mat;
    char inputfile[256] = "";
    int iMat;
    double rho;
    double T;
    /* Values returned by ANEOS. */
    double P;
    double u;
    double s;
    double cv;
    double dPdT;
    double dPdrho;
    double fkros;
    double cs;
    int Phase;
    double rhoL;
    double rhoH;
    double ion;
    /* Store differences. */
    double **pDiff;
    double **uDiff;
    double **sDiff;
    double **cDiff;
    /* Extended EOS tables. */
    int **PhaseDiff;
    /* Interpolate between two grid points. */
    int nSteps = 5;
    double eps = 1e-6;
    FILE *fp;

    if (argc != 3) {
        fprintf(stderr,"Usage: testMANEOStable <inputfile> <iMat>\n");
        exit(1);
    }

    strcpy(inputfile, argv[1]);
    iMat = atoi(argv[2]); 
    assert(iMat >= 0);

    fprintf(stderr, "Initialize M-ANEOS with input file: %s\n", inputfile);
    initaneos(inputfile);
    
    fprintf(stderr, "Initialize ANEOSmaterial with iMat: %i\n", iMat);
    /* Use cgs units. */
	Mat = ANEOSinitMaterial(iMat, 0.0, 0.0);
    ANEOSPrintMat(Mat, stderr);
    fprintf(stderr, "\n");

    fprintf(stderr, "CodeUnitstoCGSforRho = %15.7E\n", Mat->CodeUnitstoCGSforRho);
    fprintf(stderr, "CodeUnitstoCGSforU   = %15.7E\n", Mat->CodeUnitstoCGSforU);
    fprintf(stderr, "CodeUnitstoCGSforP   = %15.7E\n", Mat->CodeUnitstoCGSforP);
    fprintf(stderr, "CodeUnitstoCGSforC   = %15.7E\n", Mat->CodeUnitstoCGSforC);
    fprintf(stderr, "\n");

    /* Allocate memory. */
    pDiff = (double **)malloc(sizeof(double*)*Mat->nT);
    uDiff = (double **)malloc(sizeof(double*)*Mat->nT);
    sDiff = (double **)malloc(sizeof(double*)*Mat->nT);
    cDiff = (double **)malloc(sizeof(double*)*Mat->nT);

    PhaseDiff = (int **)malloc(sizeof(int*)*Mat->nT);

    for (int i=0; i<Mat->nT; i++)
    {
        pDiff[i] = (double *)malloc(Mat->nRho * sizeof(double));
        uDiff[i] = (double *)malloc(Mat->nRho * sizeof(double));
        sDiff[i] = (double *)malloc(Mat->nRho * sizeof(double));
        cDiff[i] = (double *)malloc(Mat->nRho * sizeof(double));

        PhaseDiff[i] = (int *)malloc(Mat->nRho * sizeof(int));
    }

#if 0
    /* Check if the data agree. */
    int nTArray = 0;
    for (int i=0; i<Mat->nT-1; i++) {
        for (int k=0; k<nSteps; k++) {
            T = Mat->TAxis[i] + k*(Mat->TAxis[i+1]-Mat->TAxis[i])/(nSteps-1);

            fprintf(stderr, "T_i= %15.7E T_i+1= %15.7E T=%15.7E\n", Mat->TAxis[i], Mat->TAxis[i+1], T);
            nTArray++;
        }
    }

    fprintf(stderr, "nTArray = %i nT= %i nSteps= %i\n", nTArray, Mat->nT, nSteps);
    fprintf(stderr, "T_min = %15.7E T_max = %15.7E\n", Mat->TAxis[0], Mat->TAxis[Mat->nT-1]);

    int nRhoArray = 0;
    for (int j=0; j<Mat->nRho-1; j++) {
        for (int l=0; l<nSteps; l++) {
            rho = Mat->rhoAxis[j] + l*(Mat->rhoAxis[j+1]-Mat->rhoAxis[j])/(nSteps-1);

            fprintf(stderr, "rho_j= %15.7E rho_j+1= %15.7E rho=%15.7E\n", Mat->rhoAxis[j], Mat->rhoAxis[j+1], rho);
            nRhoArray++;
        }
    }

    fprintf(stderr, "nRhoArray = %i nRho= %i nSteps= %i\n", nRhoArray, Mat->nRho, nSteps);
    fprintf(stderr, "rho_min = %15.7E rho_max = %15.7E\n", Mat->rhoAxis[0], Mat->rhoAxis[Mat->nRho-1]);
#endif

    for (int i=0; i<Mat->nT-1; i++) {
        for (int k=0; k<nSteps; k++) {
            T = Mat->TAxis[i] + k*(Mat->TAxis[i+1]-Mat->TAxis[i])/(nSteps-1);

            for (int j=0; j<Mat->nRho-1; j++) {
                for (int l=0; l<nSteps; l++) {
                    rho = Mat->rhoAxis[j] + l*(Mat->rhoAxis[j+1]-Mat->rhoAxis[j])/(nSteps-1);

                    /* Assume that only one material is stored in an M-ANEOS input file. */
                    callaneos_cgs(T, rho, 1, &P, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &cs,
                                  &Phase, &rhoL, &rhoH, &ion);

#if 0
                pDiff[i][j] = fabs((P-Mat->pArray[i][j])/P);
                uDiff[i][j] = fabs((u-Mat->uArray[i][j])/u);
                sDiff[i][j] = fabs((s-Mat->sArray[i][j])/s);
                cDiff[i][j] = fabs((cs-Mat->cArray[i][j])/cs);

                /* Check extended EOS tables. */
                if (Mat->PhaseArray != NULL)
                    PhaseDiff[i][j] = fabs((Phase-Mat->PhaseArray[i][j])/Phase);
#endif
                    //fprintf(stderr, "rho = %15.7E T = %15.7E\n", rho, T);
                }
            }
        }
    }

    exit(1);
    /* Print differences to a file. */
    fp = fopen("diff_press.txt", "w");
    PrintArrayDouble(pDiff, Mat->nRho, Mat->nT, fp);
    fclose(fp);

    fp = fopen("diff_u.txt", "w");
    PrintArrayDouble(uDiff, Mat->nRho, Mat->nT, fp);
    fclose(fp);

    /* Print extended EOS tables. */
    if (Mat->PhaseArray != NULL) {
        fp = fopen("diff_phase.txt", "w");
        PrintArrayInt(PhaseDiff, Mat->nRho, Mat->nT, fp);
        fclose(fp);
    }

    /* Free memory. */
	ANEOSfinalizeMaterial(Mat);

    fprintf(stderr, "Finished, exiting\n");
    return 0;
}
