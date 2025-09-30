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
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	
	if (argc != 4) {
        fprintf(stderr,"Usage: calcTemperatureRhoP <matfilename> <rhofilename> <P>\n");
        exit(1);
    }
	
	double P = atof(argv[3]);

	// Version check
	if (ANEOS_VERSION_MAJOR != 1) {
        fprintf(stderr, "ANEOS library has the wrong version (%s)\n", ANEOS_VERSION_TEXT);
        exit(1);
    }
	
	// Initialization
	ANEOSMATERIAL *material1;
	ANEOSMATERIAL *material2;
	material1 = ANEOSinitMaterial(62, dKpcUnit, dMsolUnit);
	material2 = ANEOSinitMaterial(63, dKpcUnit, dMsolUnit);

	FILE *matfile, *rhofile, *outfile;

	int n_lines, i;
	double rho;
	int mat;

	matfile = fopen(argv[1], "r");
	rhofile = fopen(argv[2], "r");
	outfile = fopen("temperature.txt", "w");
	fscanf(rhofile, "%d", &n_lines);
	fscanf(matfile, "%d", &n_lines);
	fprintf(outfile, "%d\n",n_lines);
	for (i = 0; i < n_lines; i++) {
	    fscanf(rhofile, "%lf", &rho);
	    fscanf(matfile, "%d", &mat);
	    double T = 0.0;
	    if (mat == 62) {
		T = ANEOSTofRhoP(material1, rho, P);
	    } else if (mat == 63) {
                T = ANEOSTofRhoP(material2, rho, P);
	    }

	    fprintf(outfile, "%.15e\n", T);
	}

	fclose(matfile);
	fclose(rhofile);
	fclose(outfile);

	// Finalize
    	ANEOSfinalizeMaterial(material1);
	ANEOSfinalizeMaterial(material2);
	return 0;
}
