/*
 * The header file for the ANEOS material library
 */
 
typedef struct ANEOSmaterial
{
	int iMat; // Material number
	int nRho; // Number of entries in the interpolation arrays in the rho dimension
	int nT; // Number of entries in the interpolation arrays in the T dimension
	
	// unit conversion
	double CodeUnitstoCGSforU;
	double CodeUnitstoCGSforP;
	double CodeUnitstoCGSforRho;
	double CodeUnitstoCGSforC;
	
	// interpolation arrays, array of array of pointers
	double **rhoArray;
	double **TArray;
	double **uArray;
	double **dudrhoArray;
	double **dudTArray;
	double **ddudTdrhoArray;
	double **pArray;
	double **dpdrhoArray;
	double **dpdTArray;
	double **ddpdTdrhoArray;
	double **cArray;
	double **dcdrhoArray;
	double **dcdTArray;
	double **ddcdTdrhoArray;
	double **sArray;
	double **dsdrhoArray;
	double **dsdTArray;
	double **ddsdTdrhoArray;
} ANEOSMATERIAL;

// Initialization and finalization

ANEOSMATERIAL *ANEOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit);
void ANEOSfinalizeMaterial(ANEOSMATERIAL *material);

// Access functions

double ANEOSTofPU(ANEOSMATERIAL *material, double p, double u);
double ANEOSTofRhoU(ANEOSMATERIAL *material, double rho, double u);
double ANEOSTofRhoP(ANEOSMATERIAL *material, double rho, double p);

double ANEOSUofPT(ANEOSMATERIAL *material, double p, double T);
double ANEOSUofRhoT(ANEOSMATERIAL *material, double rho, double T);
double ANEOSUofRhoP(ANEOSMATERIAL *material, double rho, double p);

double ANEOSPofUT(ANEOSMATERIAL *material, double u, double T);
double ANEOSPofRhoT(ANEOSMATERIAL *material, double rho, double T);
double ANEOSPofRhoU(ANEOSMATERIAL *material, double rho, double u);

double ANEOSRhoofUT(ANEOSMATERIAL *material, double u, double T);
double ANEOSRhoofPT(ANEOSMATERIAL *material, double p, double T);
double ANEOSRhoofPU(ANEOSMATERIAL *material, double p, double u);

double ANEOSCofRhoU(ANEOSMATERIAL *material, double rho, double u);
double ANEOSCofRhoT(ANEOSMATERIAL *material, double rho, double T);
double ANEOSisentropicU(ANEOSMATERIAL *material, double rho1, double u1, double rho2);

// Internal functions

double backwardInterpolateTemperature(double rho, double z, int nT, int nRho, double** rhoArray,
	double** TArray, double** zArray, double** dudrhoArray, double** dudTArray, double** ddudTdrhoArray);
double backwardInterpolateTemperatureBilinear(double rho, double z, int nT, int nRho, double** rhoArray,
	double** TArray, double** zArray);
double backwardInterpolateDensity(double T, double z, int nT, int nRho, double** rhoArray,
	double** TArray, double** uArray, double** dudrhoArray, double** dudTArray, double** ddudTdrhoArray);
double interpolateValue(double rho, double T, int nT, int nRho, double** rhoArray,
	double** TArray, double** zArray, double** dzdrhoArray, double** dzdTArray, double** ddzdTdrhoArray);
double interpolateValueBilinear(double rho, double T, int nT, int nRho, double** rhoArray,
	double** TArray, double** zArray);