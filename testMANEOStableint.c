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
    /* Axis for interpolation. */
    double *rhoAxisInt;
    double *TAxisInt;
    /* Store differences. */
    double **pDiff;
    double **uDiff;
    double **sDiff;
    double **cDiff;
    /* Extended EOS tables. */
    int **PhaseDiff;
    /* Number of grid points on the interpolated grid. */
    int nGridT;
    int nGridRho;
    int nGridInt = 2;
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

    /* Between two grid point of the EOS table we have nGridInt interpolation points.*/
    nGridRho = (nGridInt+1)*(Mat->nRho-1)+1;
    nGridT = (nGridInt+1)*(Mat->nT-1)+1;

    fprintf(stderr, "nGridRho = %i nGridT = %i nGridInt = %i\n", nGridRho, nGridT, nGridInt);
    fprintf(stderr, "nRho = %i nT = %i\n", Mat->nRho, Mat->nT);
    
    /* Allocate memory. */
    rhoAxisInt = (double *)malloc(nGridRho * sizeof(double));
    TAxisInt = (double *)malloc(nGridT * sizeof(double));

    pDiff = (double **)malloc(sizeof(double*)*nGridT);
    uDiff = (double **)malloc(sizeof(double*)*nGridT);
    sDiff = (double **)malloc(sizeof(double*)*nGridT);
    cDiff = (double **)malloc(sizeof(double*)*nGridT);

    PhaseDiff = (int **)malloc(sizeof(int*)*nGridT);

    for (int i=0; i<nGridT; i++)
    {
        pDiff[i] = (double *)malloc(nGridRho * sizeof(double));
        uDiff[i] = (double *)malloc(nGridRho * sizeof(double));
        sDiff[i] = (double *)malloc(nGridRho * sizeof(double));
        cDiff[i] = (double *)malloc(nGridRho * sizeof(double));

        PhaseDiff[i] = (int *)malloc(nGridRho * sizeof(int));
    }

    /* Generate rho and T axis for interpolation. */
    for (int i=0; i<Mat->nT-1; i++) {
        for (int k=0; k<=nGridInt; k++) {
            T = Mat->TAxis[i] + k*(Mat->TAxis[i+1]-Mat->TAxis[i])/(nGridInt+2);
            TAxisInt[i*(nGridInt+1)+k] = T;
        }
    }
    
    for (int j=0; j<Mat->nRho-1; j++) {
        for (int l=0; l<=nGridInt; l++) {
            rho = Mat->rhoAxis[j] + l*(Mat->rhoAxis[j+1]-Mat->rhoAxis[j])/(nGridRho+2);
            rhoAxisInt[j*(nGridInt+1)+l] = rho;
        }
    }

    /* Include the last grid point. */
    TAxisInt[nGridT-1] = Mat->TAxis[Mat->nT-1];
    rhoAxisInt[nGridRho-1] = Mat->rhoAxis[Mat->nRho-1]; 

    fprintf(stderr, "rho_min = %15.7E rho_max = %15.7E\n", rhoAxisInt[0], rhoAxisInt[nGridRho-1]);
    fprintf(stderr, "T_min   = %15.7E T_max   = %15.7E\n", TAxisInt[0], TAxisInt[nGridT-1]);

    /* Write axis to a file. */
    fp = fopen("T_axis_int.txt", "w");

    for (int i=0; i<nGridT; i++) {
        fprintf(fp, "%15.7E\n", TAxisInt[i]);
    }

    fclose(fp);

    fp = fopen("rho_axis_int.txt", "w");

    for (int j=0; j<nGridRho; j++) {
        fprintf(fp, "%15.7E\n", rhoAxisInt[j]);
    }

    fclose(fp);

    fprintf(stderr, "Check EOS table.\n");

    /* Do not include the end points of the EOS table because interpolation fails there. */
    for (int i=0; i<nGridT-1; i++) {
        for (int j=0; j<nGridRho-1; j++) {
            T = TAxisInt[i];
            rho = rhoAxisInt[j];

            /* Assume that only one material is stored in an M-ANEOS input file. */
            callaneos_cgs(T, rho, 1, &P, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &cs,
                          &Phase, &rhoL, &rhoH, &ion);
            
            //fprintf(stderr, "i=% i j= %i\n", i, j);
            
            pDiff[i][j] = fabs((ANEOSPofRhoT(Mat, rho, T)-P)/P);
            uDiff[i][j] = fabs((ANEOSUofRhoT(Mat, rho, T)-u)/u);
            //sDiff[i][j] = fabs((ANEOSSofRhoT(Mat, rho, T)-s)/s);
            cDiff[i][j] = fabs((ANEOSCofRhoT(Mat, rho, T)-cs)/cs);

            /* Check extended EOS tables. */
            if (Mat->PhaseArray != NULL) {
                if (Phase == ANEOSPhaseofRhoT(Mat, rho, T)) {
                    PhaseDiff[i][j] = 0;
                } else {
                    PhaseDiff[i][j] = 1;
                }
            }
        }
    }

    fprintf(stderr, "\n");

    /* Print differences to a file. */
    fp = fopen("diff_press_int.txt", "w");
    PrintArrayDouble(pDiff, nGridRho-1, nGridT-1, fp);
    fclose(fp);

    fp = fopen("diff_u_int.txt", "w");
    PrintArrayDouble(uDiff, nGridRho-1, nGridT-1, fp);
    fclose(fp);

    /* Print extended EOS tables. */
    if (Mat->PhaseArray != NULL) {
        fp = fopen("diff_phase_int.txt", "w");
        PrintArrayInt(PhaseDiff, nGridRho-1, nGridT-1, fp);
        fclose(fp);
    }

    /* Free memory. */
	ANEOSfinalizeMaterial(Mat);

    fprintf(stderr, "Finished, exiting\n");
    return 0;
}
