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


    /* Check if the data agree. */
    for (int i=0; i<Mat->nT; i++) {
        for (int j=1; j<Mat->nRho; j++) {
            T = Mat->TAxis[i];
            rho = Mat->rhoAxis[j];

            /* Assume that only one material is stored in an M-ANEOS input file. */
            callaneos_cgs(T, rho, 1, &P, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &cs,
                          &Phase, &rhoL, &rhoH, &ion);

            pDiff[i][j] = fabs((P-Mat->pArray[i][j])/P);
            uDiff[i][j] = fabs((u-Mat->uArray[i][j])/u);
            sDiff[i][j] = fabs((s-Mat->sArray[i][j])/s);
            cDiff[i][j] = fabs((cs-Mat->cArray[i][j])/cs);

            /* Check extended EOS tables. */
            if (Mat->PhaseArray != NULL)
                PhaseDiff[i][j] = fabs((Phase-Mat->PhaseArray[i][j])/Phase);
 
#if 0
            if (dDiffP >= eps)
                fprintf(stderr, "rho= %g T= %g P= %g %g dDiffP= %g\n", rho, T, P, Mat->pArray[i][j], dDiffP);

            if (dDiffu >= eps)
                fprintf(stderr, "rho= %g T= %g u= %g %g dDiffu= %g\n", rho, T, u, Mat->uArray[i][j], dDiffu);

            assert(dDiffP < eps);
            assert(dDiffu < eps);
            assert(dDiffs < eps);
            assert(dDiffcs < eps);
#endif
            //PhaseArray[i][j]
        }
    }

    /* Print differences to a file. */
    fp = fopen("diff_press.txt", "w");

    for (int i=0; i<Mat->nT; i++) {
        for (int j=1; j<Mat->nRho; j++) {
            fprintf(fp, "%15.7E", pDiff[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);

    fp = fopen("diff_u.txt", "w");

    for (int i=0; i<Mat->nT; i++) {
        for (int j=1; j<Mat->nRho; j++) {
            fprintf(fp, "%15.7E", uDiff[i][j]);
        }
        fprintf(fp, "\n");
    }

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
