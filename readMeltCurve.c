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
 * ANEOS material model
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

int main(int argc, char *argv[]) {
    double *T;
	int nRho;
    char inputfile[256] = "";
	FILE *file;

    if (argc != 2) {
        fprintf(stderr,"Usage: readMeltCurve <melt_curve.in>\n");
        exit(1);
    }

    strcpy(inputfile, argv[1]);

	if ((file= fopen(inputfile, "rb")) == NULL) {
		fprintf(stderr, "Could not open file %s\n", inputfile);
		assert(0);
	}

	if (fread(&nRho, sizeof(nRho), 1, file) == 0) {
		fprintf(stderr, "Failed to read from file %s\n", inputfile);
		assert(0);
	}

    T = (double *) calloc(nRho, sizeof(double));

	if (fread(T, sizeof(T[0]), nRho, file) == 0) {
	    fprintf(stderr, "Failed to read from file %s\n", inputfile);
	    assert(0);
	}
	
	fclose(file);

    printf("# nRho= %i\n", nRho);

    for (int i=0; i<nRho; i++) {
        printf("%15.7E\n", T[i]);
    }

    return 0;
}
