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
 * Writes the phase data on a grid for the ANEOS material model
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "ANEOSmaterial.h"

int main(int argc, char *argv[])
{
    ANEOSMATERIAL *Mat;
    // Use cgs units.
    double dKpcUnit = 0.0;
    double dMsolUnit = 0.0;
    double rho;
    double T;
    double p;
    FILE *fp;

	if (argc != 2) {
        fprintf(stderr,"Usage: writePhaseTable <iMat>\n");
        exit(1);
    }
	
	int iMat = atoi(argv[1]);
	
    Mat =  ANEOSinitMaterial(iMat, dKpcUnit, dMsolUnit);
    assert(Mat != NULL);

    // Check if the extended EOS table was loaded
    if (Mat->PhaseArray == NULL) {
        fprintf(stderr, "Extended EOS table not loaded.\n");
        exit(1);
    }

    fprintf(stderr, "iMat = %i:  %s (%s)\n", Mat->iMat, Mat->matName, Mat->matstring);
    fprintf(stderr, "nRho = %i nT = %i\n", Mat->nRho, Mat->nT);
    
    fp = fopen("phase.txt", "w");
    assert(fp != NULL);

    fprintf(fp, "# %s\n", Mat->matstring);
    fprintf(fp, "# iMat = %i nRho = %i nT = %i\n", iMat, Mat->nRho, Mat->nT);

	for (int i=0; i<Mat->nT; i++)
	{
		if ( (i%50)==0){
			fprintf(stderr, "Iterating over T for %d of %d\n", i, Mat->nT);
		}
		for (int j=0; j<Mat->nRho; j++)
		{
            fprintf(fp,"%15.7E %15.7E %i\n", Mat->rhoAxis[j], Mat->TAxis[i], Mat->PhaseArray[i][j]);
		}
	}
    fclose(fp);

    return 0;
}
