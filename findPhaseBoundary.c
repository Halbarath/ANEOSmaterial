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
	/* Variables needes for ANEOS */
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
	/* Define limits from axes.in */
	//double rho_min = 1e-25;
	//double rho_max = 1e3;
	double rho_min = 2.0;
	double rho_max = 1e3;
	double T_min = 1.0;
	int nRho = 10000;
	//int nRho = 1402;
	//int nT = 1601;

	if (argc != 2) {
        fprintf(stderr,"Usage: findPhaseBoundary <aneos.input>\n");
        exit(1);
    }
	
	/* In M-ANEOS each input file contains one material so the iMat number is irrelevant */
	int iMat = 1;
	initaneos(argv[1]);

	//rho_min = atof(argv[1]);
	//rho_max = atof(argv[2]);

	double *rhoAxis = (double *)malloc(nRho * sizeof(double));
	double *Tm = (double *)malloc(nRho * sizeof(double));
	//double *TAxis = (double *)malloc(nT * sizeof(double));
	
	
	for (int i=0; i<nRho; i++) {
		rhoAxis[i] = rho_min + i*(rho_max - rho_min)/(nRho-1);
	}

	for (int i=0; i<nRho; i++) {
		fprintf(stderr, "%15.7E\n", rhoAxis[i]);
	}

	for (int i=0; i<nRho; i++) {
		double T = T_min;
		double Ta = T_min;
		double Tb = T_min;
		double rho = rhoAxis[i];

		callaneos_cgs(T, rho, iMat, &p, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &c, &phase, &rhoL, &rhoH, &ion);

		fprintf(stderr, "rho= %15.7E T= %15.7E phase= %i\n", rho, T, phase);

		/* Check if the is no solid phase even at T_min*/
		if (phase != phase_id_solid) {
			Tm[i] = 0;
			continue;
		}

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

		Tm[i] = T;
	}

	FILE *fp;

	fp = fopen("phase_solid.txt", "w");

	fprintf(fp, "#%14s%15s\n", "rho_m", "T_m");
	
	for (int i=0; i<nRho; i++) {
		fprintf(fp, "%15.7E%15.7E\n", rhoAxis[i], Tm[i]);
	}

	fclose(fp);

    return 0;
}
