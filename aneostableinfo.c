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