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
 * Calculate P(rho, T) from M-ANEOS directly.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "aneos.h"

int main(int argc, char *argv[]) {
	double T;
    double rho;
	double p;
	double u;
	double s;
    double cv;
	double dPdT;
	double dPdrho;
    double fkros;
	double c;
	int phase;
    double rhoL;
    double rhoH;
    double ion;
	/* In M-ANEOS each input file contains one material so the iMat number is irrelevant */
	int iMat = 1;

	if (argc != 4) {
        fprintf(stderr,"Usage: calcPressureRhoTMANEOS <rho> <T> <aneos.input>\n");
        exit(1);
    }

	rho = atof(argv[1]);
	T = atof(argv[2]);

	assert(rho > 0.0);
	assert(T > 0.0);

	initaneos(argv[3]);
					
	callaneos_cgs(T, rho, iMat, &p, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &c, &phase, &rhoL, &rhoH, &ion);

	printf("rho=%15.7E T=%15.7E P=%15.7E (cgs)\n", rho, T, p);

    return 0;
}
