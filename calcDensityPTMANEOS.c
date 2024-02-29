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
 * Calculate the density rho(P, T) directly from M-ANEOS.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "aneos.h"
#include "aneosdirect.h"

int main(int argc, char *argv[]) {
	/* In M-ANEOS each input file contains one material so the iMat number is irrelevant */
    int iMat = 1;

	if (argc != 4) {
        fprintf(stderr,"Usage: calcDensityPTMANEOS <P> <T> <aneos.input>\n");
        fprintf(stderr,"Note: Units are cgs and K\n");
        exit(1);
    }
	
	double P = atof(argv[1]);
	double T = atof(argv[2]);
 
	initaneos(argv[3]);

    printf("Input:\n");
    printf("P= %15.7E erg/cm^3\n", P);
    printf("T= %15.7E K\n", T);
    printf("\n");

	double rho = ANEOSRhoofPTDirect(iMat, P, T, 1e-25, 1e3);
	
    printf("Output:\n");
    printf("rho= %15.7E g/cm^3\n", rho);
    printf("\n");
    
    // Check if the results agree
#ifdef DEBUG
    double P_inv = ANEOSPofRhoTDirect(iMat, rho, T);
    double P_err = fabs((P - P_inv)/P);
    printf("rho=%15.7E T=%15.7E P=%15.7E err=%15.7E\n", rho, T, P_inv, P_err);
#endif

	return 0;
}
