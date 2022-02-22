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
 * Calculate all thermodynamic variables from M-ANEOS.
 *
 * Author:   Christian Reinhardt
 * Created:  12.06.2021
 * Modified: 22.02.2021 
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "aneos.h"

int main(int argc, char **argv) {
    double rho;
    double T;
    /* In M-ANEOS files there is only one material per file. */
    int iMat = 1;
    char MatFile[256];
    double p;
    double u;
    double s;
    double cv;
    double dPdT;
    double dPdrho;
    double fkros;
    double cs;
    int iPhase;
    double rhoL;
    double rhoH;
    double ion;

    if (argc != 4) {
        fprintf(stderr, "Usage: maneoscall <rho> <T> <input file>\n");
        exit(1);
    }

    rho = atof(argv[1]);
    T = atof(argv[2]);
    strcpy(MatFile, argv[3]);

    assert(rho > 0.0);
    assert(T > 0.0);
    assert(MatFile != NULL);

    assert(strlen(MatFile) > 0);

    fprintf(stderr, "ANEOS: Initializing material...\n");
    initaneos(MatFile);

    callaneos_cgs(T, rho, iMat, &p, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &cs, &iPhase, &rhoL, &rhoH,
                  &ion);

    printf("Input:\n");
    printf("rho = %15.7E\n", rho);
    printf("T   = %15.7E\n", T);
    printf("iMat= %i\n", iMat);
    printf("\n");

    printf("Output:\n");
    printf("p      = %15.7E\n", p);
    printf("u      = %15.7E\n", u);
    printf("s      = %15.7E\n", s);
    printf("cv     = %15.7E\n", cv);
    printf("dPdT   = %15.7E\n", dPdT);
    printf("dPdrho = %15.7E\n", dPdrho);
    printf("fkros  = %15.7E\n", fkros);
    printf("cs     = %15.7E\n", cs);
    printf("iPhase = %15i\n", iPhase);
    printf("rhoL   = %15.7E\n", rhoL);
    printf("rhoH   = %15.7E\n", rhoH);
    printf("ion    = %15.7E\n", ion);
    printf("\n");

    return 0;
}
