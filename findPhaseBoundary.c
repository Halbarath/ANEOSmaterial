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
 * Find phase boundary from ANEOS
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "aneos.h"

int main(int argc, char *argv[])
{
	/* Variables needed for ANEOS */
	//double T;
    //double rho;
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
	/* Phase id for a solid */
	const int phase_id_solid = 4;
	/* rho is read from axes.in */
	char axesFilename[256] = "axes.in";
	double *rhoAxis;
	double T_min = 1.0;
	int nRho;
	char str[1000];
	FILE *fp;
	/* Store T_melt(rho) */
	double *rho_melt;
	double *T_melt;
	int nMelt = 0;

	if (argc != 2) {
        fprintf(stderr,"Usage: findPhaseBoundary <aneos.input>\n");
        exit(1);
    }
	
	/* In M-ANEOS each input file contains one material so the iMat number is irrelevant */
	int iMat = 1;
	initaneos(argv[1]);

	fp = fopen(axesFilename, "r");
	
	if (fp == NULL) {
        fprintf(stderr,"Could not open file %s",axesFilename);
	}

	if (fgets(str, 1000, fp) != NULL) {
		nRho = (int) strtol(str, (char **)NULL, 10);
	}

	/* nT read from the input file is not used here! */
	if (fgets(str, 1000, fp) != NULL) {
		int tmp = (int) strtol(str, (char **)NULL, 10);
	}

	rhoAxis = (double *)malloc(nRho * sizeof(double));
	T_melt = (double *)malloc(nRho * sizeof(double));
	rho_melt = (double *)malloc(nRho * sizeof(double));
	
	for (int i=0; i<nRho; i++) {
		if (fgets(str, 1000, fp) != NULL) {
			rhoAxis[i] = (double) strtod(str, (char **)NULL);
		}
	}

	fclose(fp);

	fprintf(stderr, "Read file %s: nRho= %i\n", axesFilename, nRho);

	for (int i=0; i<nRho; i++) {
		double T = T_min;
		double Ta = T_min;
		double Tb = T_min;
		double rho = rhoAxis[i];

		callaneos_cgs(T, rho, iMat, &p, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &c, &phase, &rhoL, &rhoH, &ion);

		fprintf(stderr, "rho= %15.7E T= %15.7E phase= %i\n", rho, T, phase);

		/* Check if the is no solid phase even at T_min*/
		if (phase != phase_id_solid) continue;

		/* Find upper limit */
		while (phase == phase_id_solid) {
			Ta = Tb;
			Tb *= 2.0;

			callaneos_cgs(Tb, rho, iMat, &p, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &c, &phase, &rhoL, &rhoH, &ion);

			fprintf(stderr, "rho= %15.7E Tb= %15.7E phase= %i\n", rho, Tb, phase);
		}

		fprintf(stderr, "Root bracketed: rho=%15.7E Ta=%15.7E Tb=%15.7E\n", rho, Ta, Tb);

		while ((Tb - Ta)/T > 1e-8) {
			T = 0.5*(Ta + Tb);

			callaneos_cgs(T, rho, iMat, &p, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &c, &phase, &rhoL, &rhoH, &ion);

			if (phase == phase_id_solid) {
				Ta = T;
			} else {
				Tb = T;
			}
		
			fprintf(stderr, "Root bracketed: rho=%15.7E T=%15.7E Ta=%15.7E Tb=%15.7E phase= %i\n", rho, T, Ta, Tb, phase);
		}

		T_melt[nMelt] = T;
		rho_melt[nMelt] = rho;
		nMelt++;
	}

	fp = fopen("phase_solid.txt", "w");

	fprintf(fp, "# n = %i\n", nMelt);
	fprintf(fp, "#%14s%15s\n", "rho_m [g/cm^3]", "T_m [K]");

	for (int i=0; i<nMelt; i++) {
		fprintf(fp, "%15.7E%15.7E\n", rho_melt[i], T_melt[i]);
	}

	fclose(fp);

    return 0;
}
