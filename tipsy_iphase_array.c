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
 * This program reads a Tipsy binary file and generates a tipsy array file that contains the phase
 * of each particle calculated from ANEOS.
 *
 * Author:   Christian Reinhardt
 * Created:  05.12.2020
 * Modified: 06.03.2022 
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "tipsy.h"
#include "ANEOSmaterial.h"

int main(int argc, char **argv) {
    /* Tipsy library. */
    TCTX in;
    int N;
    /* ANEOS material. */
    ANEOSMATERIAL **Mat;
    int iMat;
    double dKpcUnit = 2.06701e-13;
    double dMsolUnit = 4.80438e-08;
    double rho, u;
    int i;

    if (argc != 2) {
        fprintf(stderr,"Usage: tipsy_iphase_array tipsy.std\n");
        exit(1);
    }
    
    char inputfile[256] = "";
    strcpy(inputfile, argv[1]);
    
    /* Allocate memory. */
    Mat = (ANEOSMATERIAL **) calloc(ANEOS_N_MATERIAL_MAX, sizeof(ANEOSMATERIAL *));
    assert(Mat != NULL);

    /* Initialize tipsy library */
    TipsyInitialize(&in, 0, inputfile);

    N = iTipsyNumParticles(in);

    /* Read all particles. */
	TipsyReadAll(in);
    
    /* Print header for the array file. */
    printf("%i\n",N);

    for (i = 0; i < N; i++) {
        iMat = (int) in->gp[i].metals;
    
        if (Mat[iMat] == NULL) {
            fprintf(stderr, "Initalizing iMat %i.\n", iMat);
            Mat[iMat] = ANEOSinitMaterial(iMat, dKpcUnit, dMsolUnit);
            assert(Mat[iMat] != NULL);
            assert(Mat[iMat]->PhaseArray != NULL);
        }

        rho = in->gp[i].rho;
        u = in->gp[i].temp;

        
        printf("%i\n", ANEOSPhaseofRhoU(Mat[iMat], rho, u));
    }

    for (int i=0; i<ANEOS_N_MATERIAL_MAX; i++) {
        if (Mat[i] != NULL) free(Mat[i]);
    }
    free(Mat);

    TipsyFinish(in);

    return 0;
}
