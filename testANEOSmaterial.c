#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ANEOSmaterial.h"

int main(int argc, char *argv[])
{
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	
	ANEOSMATERIAL *material;
	
	material = ANEOSinitMaterial(4, dKpcUnit, dMsolUnit);
	
	double P = ANEOSPofRhoU(material, 8/material->CodeUnitstoCGSforRho, 1e12/material->CodeUnitstoCGSforU);
	printf("Pressure %.15e\n", P*material->CodeUnitstoCGSforP);
	
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
	printf("testU1 %.15e\n", testU1*material->CodeUnitstoCGSforU);
	
	double testRho = ANEOSRhoofUT(material,testU,273);
	printf("testRho %.15e\n", testRho*material->CodeUnitstoCGSforRho);
	
	double testU2 = ANEOSUofRhoT(material,testRho,273);
	printf("testU2 %.15e\n", testU2*material->CodeUnitstoCGSforU);
	
	
	double peeee = ANEOSPofRhoT(material,8/material->CodeUnitstoCGSforRho,4000);
	double uuuuu = ANEOSUofRhoT(material,8/material->CodeUnitstoCGSforRho,4000);
	double rhoooo = ANEOSRhoofPU(material,peeee,uuuuu);
	double Teeee = ANEOSTofPU(material,peeee,uuuuu);
	
	printf("rhoooo %.15e\n", rhoooo*material->CodeUnitstoCGSforRho);
	printf("Teeee %.15e\n", Teeee);
	
	
	ANEOSfinalizeMaterial(material);
	return 0;
}

