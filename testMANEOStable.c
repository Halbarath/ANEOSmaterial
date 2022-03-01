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
#include "ANEOSmaterial.h"
#include "aneos.h"

int main(int argc, char *argv[])
{
    ANEOSmaterial *Mat;
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    char inputfile[256] = "";
    /* Assume that each M-ANEOS input file only contains one material.*/
    int iMatMANEOS = 1;
    int iMat;

    if (argc != 3) {
        fprintf(stderr,"Usage: testMANEOStable <inputfile> <iMat>\n");
        exit(1);
    }

    strcpy(inputfile, argv[1]);
    iMat = atoi(argv[2]); 
    assert(iMat >= 0)

    fprintf(stderr, "Initialize M-ANEOS with input file: %s\n", inputfile);
    initaneos(inputfile);
    
    fprintf(stderr, "Initialize ANEOSmaterial with iMat: %i\n", iMat);
	Mat = ANEOSinitMaterial(iMat, dKpcUnit, dMsolUnit);	
    ANEOSPrintMat(Mat, stderr);
    fprintf(stderr, "\n");

#if 0
    /* Convert all data to cgs. */
    for (int i=0; i<Mat->nT; i++) {
        for (int j=0; j<nRho; j++) {
            rho = Mat->rhoAxis[j]*;

            callaneos_cgs(T, rho, iMat, &pArray[i][j], &uArray[i][j], &sArray[i][j], &cv, &dPdT,
                    &dPdrho, &fkros, &cArray[i][j], &PhaseArray[i][j], &rhoL, &rhoH, &ion);

            TArray[i][j] = T;
            rhoArray[i][j] = rho;
        }
    }

    /* Check if the data agree. */
    for (int i=0; i<Mat->nT; i++) {
        for (int j=0; j<nRho; j++) {
            T = TAxis[i];
            rho = rhoAxis[j];

            callaneos_cgs(T, rho, iMat, &pArray[i][j], &uArray[i][j], &sArray[i][j], &cv, &dPdT,
                    &dPdrho, &fkros, &cArray[i][j], &PhaseArray[i][j], &rhoL, &rhoH, &ion);

            TArray[i][j] = T;
            rhoArray[i][j] = rho;
        }
    }
#endif

	ANEOSfinalizeMaterial(Mat);

    fprintf(stderr, "Finished, exiting\n");
    return 0;
}
