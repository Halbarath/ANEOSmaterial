#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ANEOSmaterial.h"

int main(int argc, char *argv[])
{
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	
	if (ANEOS_VERSION_MAJOR != 1) {
		fprintf(stderr, "main: ANEOS library has the wrong version (%s)\n", ANEOS_VERSION_TEXT);
		exit(1);
	}

	ANEOSMATERIAL *material;
	
	material = ANEOSinitMaterial(54, dKpcUnit, dMsolUnit);
	
	double P = ANEOSPofRhoU(material, 8/material->CodeUnitstoCGSforRho, 1e12/material->CodeUnitstoCGSforU);
	printf("Pressure %.15e\n", P*material->CodeUnitstoCGSforP);
	
	printf("CodeUnitstoCGSforRho %.15e\n", material->CodeUnitstoCGSforRho);
	printf("CodeUnitstoCGSforU %.15e\n", material->CodeUnitstoCGSforU);
	
	double U = ANEOSUofRhoP(material, 8/material->CodeUnitstoCGSforRho,P);
	printf("Energy %.15e\n", U*material->CodeUnitstoCGSforU);
	
	P = ANEOSPofRhoT(material, 8/material->CodeUnitstoCGSforRho, 273);
	printf("Pressure %.15e\n", P*material->CodeUnitstoCGSforP);
	
	double C = ANEOSCofRhoU(material, 8/material->CodeUnitstoCGSforRho, 1e12/material->CodeUnitstoCGSforU);
	printf("Sound speed %.15e\n", C*material->CodeUnitstoCGSforC);
	
	double newu = ANEOSisentropicU(material, 8/material->CodeUnitstoCGSforRho, 1e12/material->CodeUnitstoCGSforU, 10/material->CodeUnitstoCGSforRho);
	printf("new U %.15e\n", newu*material->CodeUnitstoCGSforU);
	
	
	
	double testU = ANEOSUofRhoT(material,8/material->CodeUnitstoCGSforRho,273); //7.8623669
	double testT = ANEOSTofRhoU(material,8/material->CodeUnitstoCGSforRho,testU);
	double testU1 = ANEOSUofRhoT(material,8/material->CodeUnitstoCGSforRho,testT);
	printf("testU %.15e\n", testU*material->CodeUnitstoCGSforU);
	printf("testT %.15e\n", testT);
	printf("testU1 %.15e should be %.15e\n", testU1*material->CodeUnitstoCGSforU, testU*material->CodeUnitstoCGSforU);
	
	double testRho = ANEOSRhoofUT(material,testU,273);
	printf("testRho %.15e, should be 8\n", testRho*material->CodeUnitstoCGSforRho);
	
	double testU2 = ANEOSUofRhoT(material,testRho,273);
	printf("testU2 %.15e should be %.15e\n", testU2*material->CodeUnitstoCGSforU, testU*material->CodeUnitstoCGSforU);
	
	
	double peeee = ANEOSPofRhoT(material,8/material->CodeUnitstoCGSforRho,4000);
	double uuuuu = ANEOSUofRhoT(material,8/material->CodeUnitstoCGSforRho,4000);
	double rhoooo = ANEOSRhoofPU(material,peeee,uuuuu);
	double Teeee = ANEOSTofPU(material,peeee,uuuuu);
	
	printf("rhoooo %.15e, should be 8\n", rhoooo*material->CodeUnitstoCGSforRho);
	printf("Teeee %.15e, should be 4000\n", Teeee);
	
	double dPdRho = ANEOSdPdRhoofRhoU(material, 8/material->CodeUnitstoCGSforRho, 1e12/material->CodeUnitstoCGSforU);
	printf("dPdRho %.15e\n", dPdRho);
	
	double dPdU = ANEOSdPdUofRhoU(material, 8/material->CodeUnitstoCGSforRho, 1e12/material->CodeUnitstoCGSforU);
	printf("dPdU %.15e\n", dPdU);
	
	double dUdRho = ANEOSdUdRhoofRhoU(material, 8/material->CodeUnitstoCGSforRho, 1e12/material->CodeUnitstoCGSforU);
	printf("dUdRho %.15e\n", dUdRho);
	
	printf("\n");
	double testfail = ANEOSTofRhoU(material, 0/material->CodeUnitstoCGSforRho, 1e12/material->CodeUnitstoCGSforU);
	testfail = ANEOSTofRhoU(material, 1000/material->CodeUnitstoCGSforRho, 1e12/material->CodeUnitstoCGSforU);
	testfail = ANEOSTofRhoU(material, 8/material->CodeUnitstoCGSforRho, 0/material->CodeUnitstoCGSforU);
	testfail = ANEOSTofRhoU(material, 8/material->CodeUnitstoCGSforRho, 1e25/material->CodeUnitstoCGSforU);
	printf("\n");

	testfail = ANEOSRhoofUT(material, 1e12/material->CodeUnitstoCGSforU, 0);
	testfail = ANEOSRhoofUT(material, 1e12/material->CodeUnitstoCGSforU, 1e10);
	testfail = ANEOSRhoofUT(material, 0/material->CodeUnitstoCGSforU, 100);
	testfail = ANEOSRhoofUT(material, 1e25/material->CodeUnitstoCGSforU, 100);
	printf("\n");

	testfail = ANEOSUofRhoT(material, 0/material->CodeUnitstoCGSforRho, 100);
	testfail = ANEOSUofRhoT(material, 1000/material->CodeUnitstoCGSforRho, 100);
	testfail = ANEOSUofRhoT(material, 8/material->CodeUnitstoCGSforRho, 0);
	testfail = ANEOSUofRhoT(material, 8/material->CodeUnitstoCGSforRho, 1e10);

	char matName[256];
	ANEOSMatString(material, matName);
	printf("Mat name %s\n",matName);

	ANEOSprintMat(material, stderr);

	ANEOSfinalizeMaterial(material);
	return 0;
}

