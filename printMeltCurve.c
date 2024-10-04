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
 * Print the melt curve T_melt(rho) for a material.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "ANEOSmaterial.h"

int main(int argc, char *argv[]) {
	double dKpcUnit = 0.0;
	double dMsolUnit = 0.0;
	
	if (argc != 2) {
        fprintf(stderr,"Usage: printMeltCurve <iMat>\n");
        exit(1);
    }
	
	int iMat = atoi(argv[1]);

	// Version check
	if (ANEOS_VERSION_MAJOR != 1) {
        fprintf(stderr, "ANEOS library has the wrong version (%s)\n", ANEOS_VERSION_TEXT);
        exit(1);
    }
	
	// Initialization
	ANEOSMATERIAL *material;
	material = ANEOSinitMaterial(iMat, dKpcUnit, dMsolUnit);

    printf("# iMat= %i nRho= %i\n", material->iMat, material->nRho);
    printf("#%14s%15s%15s\n", "rho", "T", "P");

    for (int i=0; i<material->nRho; i++) {
		double rho = material->rhoAxis[i];
		double T_melt = ANEOSTmeltofRho(material, rho);
        double P_melt = ANEOSPofRhoT(material, rho, T_melt);

        // Set P = 0 if T_melt = 0
        if (P_melt < 0.0) P_melt = 0.0;

        printf("%15.7E%15.7E%15.7E\n", rho, T_melt, P_melt);
    }

    return 0;
}
