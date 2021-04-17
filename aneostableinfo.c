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
 * Display information, e.g., number of data points or density and temperature values, of an ANEOS
 * table generated with writeANEOStable.c.
 *
 * Author:   Christian Reinhardt
 * Created:  12.03.2021
 * Modified:  
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ANEOSmaterial.h"

int main(int argc, char **argv) { 
	char inputfile[256];
    double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;

    if (argc != 2) {
        fprintf(stderr, "Usage: aneostableinfo <ANEOStable.in>\n");
        exit(1);
    }

	strcpy(inputfile, argv[1]);

	ANEOSMATERIAL *material;
	
	material = ANEOSinitMaterialFromFile(1, inputfile, dKpcUnit, dMsolUnit);
    ANEOSPrintMat(material, stderr);

    return 0;
}