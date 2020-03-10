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
	double *rhoAxis;
	double *TAxis;
	double **uArray;
	double **pArray;
	double **cArray;
	double **sArray;
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

double ANEOSdPdRhoofRhoU(ANEOSMATERIAL *material, double rho, double u);
double ANEOSdPdUofRhoU(ANEOSMATERIAL *material, double rho, double u);
double ANEOSdUdRhoofRhoU(ANEOSMATERIAL *material, double rho, double u);

// Internal functions

double backwardInterpolateTemperatureBilinear(double rho, double z, int nT, int nRho, double* rhoAxis,
	double* TAxis, double** zArray);
double backwardInterpolateDensityBilinear(double T, double z, int nT, int nRho, double* rhoAxis,
	double* TAxis, double** uArray);
double interpolateValueBilinear(double rho, double T, int nT, int nRho, double* rhoAxis,
	double* TAxis, double** zArray);