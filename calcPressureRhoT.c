/*  This file is part of EOSlib.
 *  Copyright (c) 2020-2021 Thomas Meier & Christian Reinhardt
 *
 *  EOSlib is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EOSlib is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EOSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Calculate the presssure P(rho, T).
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "ANEOSmaterial.h"

int main(int argc, char *argv[])
{
	//double dKpcUnit = 2.06701e-13;
	//double dMsolUnit = 4.80438e-08;
	double dKpcUnit = 0.0;
	double dMsolUnit = 0.0;
	
	if (argc != 4) {
        fprintf(stderr,"Usage: calcPressureRhoT <iMat> <rho> <T>\n");
        exit(1);
    }
	
	int iMat = atoi(argv[1]);
	double rho = atof(argv[2]);
	double T = atof(argv[3]);

    assert(rho > 0.0);
    assert(T >= 1.0);
    
	// Version check
	if (ANEOS_VERSION_MAJOR != 1) {
        fprintf(stderr, "ANEOS library has the wrong version (%s)\n", ANEOS_VERSION_TEXT);
        exit(1);
    }
	
	// Initialization
	ANEOSMATERIAL *material;
	material = ANEOSinitMaterial(iMat, dKpcUnit, dMsolUnit);

    printf("Input:\n");
    printf("rho= %15.7E g/cm^3\n", rho);
    printf("T= %15.7E K\n", T);
    printf("\n");

    // Convert input parameters to code units
    rho /= material->CodeUnitstoCGSforRho;

	double P = ANEOSPofRhoT(material, rho, T);

    printf("Output:\n");
    printf("P= %15.7E code units = %15.7E erg/cm^3\n", P, P*material->CodeUnitstoCGSforP);
    printf("\n");
		
	// Finalize
    ANEOSfinalizeMaterial(material);
	return 0;
}
