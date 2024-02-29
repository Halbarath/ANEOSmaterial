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
 * Test the function in aneosdirect.c that directly call M-ANEOS.
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
    int nRho;
    int nT;
    char inputFile[256] = "aneos-iron-stewart-2020.input";   
	char axesFile[256] = "axes.in";
	char str[1000];
    FILE *fp;

	initaneos(inputFile);

	fp = fopen(axesFile, "r");
	
	if (fp == NULL) {
        fprintf(stderr,"Could not open file %s",axesFile);
	}

	if (fgets(str, 1000, fp) != NULL) {
		nRho = (int) strtol(str, (char **)NULL, 10);
	}

	if (fgets(str, 1000, fp) != NULL) {
		nT = (int) strtol(str, (char **)NULL, 10);
	}

	double *rho = (double *) calloc(nRho, sizeof(double));
	double *T = (double *) calloc(nT, sizeof(double));

    for (int i=0; i<nRho; i++)
    {
        if (fgets(str, 1000, fp) != NULL)
        {
            rho[i] = (double) strtod(str, (char **)NULL);
        }
    }

    for (int i=0; i<nT; i++)
    {
        if (fgets(str, 1000, fp) != NULL)
        {
            T[i] = (double) strtod(str, (char **)NULL);
        }
    }

    fclose(fp);

	fprintf(stderr, "Read file %s: nRho= %i nT= %i\n", axesFile, nRho, nT);
    fprintf(stderr, "rho_min=%15.7E rho_max=%15.7E\n", rho[0], rho[nRho-1]);
    fprintf(stderr, "T_min=%15.7E T_max=%15.7E\n", T[0], T[nT-1]);

    for (int i=0; i<nT; i++) {
        for (int j=1; j<nRho; j++) {
            /* For low densities root finding can fail at phase boundaries. */
            if (rho[i] < 1e0) continue;

            double p = ANEOSPofRhoTDirect(iMat, rho[i], T[i]);
            double rho_inv = ANEOSRhoofPTDirect(iMat, p, T[i], rho[1], rho[nRho-1]);
            double p_inv = ANEOSPofRhoTDirect(iMat, rho_inv, T[i]);
            double err = (rho[i] - rho_inv)/rho[i];
            double err_P = (p - p_inv)/p;

            assert (rho_inv > 0.0);

            if (err_P > 1e-6) {
                fprintf(stderr, "Inversion failed:\n");
                fprintf(stderr, "rho=%15.7E T=%15.7E p=%15.7E\n", rho[i], T[i], p);
                fprintf(stderr, "rho_inv=%15.7E p_inv=%15.7E err_P=%15.7E\n", rho_inv, p_inv, err_P);
                assert(0);
            }

            if (err > 1e-6) {
                fprintf(stderr, "Inversion failed:\n");
                fprintf(stderr, "rho=%15.7E T=%15.7E p=%15.7E\n", rho[i], T[i], p);
                fprintf(stderr, "rho_inv=%15.7E p_inv=%15.7E err=%15.7E\n", rho_inv, p_inv, err);
                assert(0);
            }
        }
    }

	return 0;
}
