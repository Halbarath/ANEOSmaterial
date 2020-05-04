/*
 * ANEOS material library
 * Header of interpBilinear.c
 *
 */

#ifndef INTERPBILINEAR_HINCLUDED
#define INTERPBILINEAR_HINCLUDED

double backwardInterpolateTemperatureBilinear(double rho, double z, int nT, int nRho, double* rhoAxis,
	double* TAxis, double** zArray);
double backwardInterpolateDensityBilinear(double T, double z, int nT, int nRho, double* rhoAxis,
	double* TAxis, double** uArray);
double interpolateValueBilinear(double rho, double T, int nT, int nRho, double* rhoAxis,
	double* TAxis, double** zArray);
int findIndex(double x, double* xAxis, int nX);
#endif