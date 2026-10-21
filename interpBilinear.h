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
 * Header of interpBilinear.c
 *
 */

#ifndef INTERPBILINEAR_HINCLUDED
#define INTERPBILINEAR_HINCLUDED

double backwardInterpolateTemperatureBilinear(double rho, double z, int nT, int nRho, double* rhoAxis,
	double* TAxis, double** zArray);
double backwardInterpolateTemperatureBilinear_Scan(double rho, double z, int nT, int nRho, double* rhoAxis,
	double* TAxis, double** zArray);
double backwardInterpolateDensityBilinear(double T, double z, int nT, int nRho, double* rhoAxis,
	double* TAxis, double** uArray);
double interpolateValueBilinear(double rho, double T, int nT, int nRho, double* rhoAxis,
	double* TAxis, double** zArray);
int interpolateValueNearest(double rho, double T, int nT, int nRho, double* rhoAxis,
        double* TAxis, int** zArray);
int findIndex(double x, double* xAxis, int nX);
double interpolateValueLinear(double x, int nX, double* xAxis, double* yAxis);
int interpolateGradientBilinear(double rho, double T, int nT, int nRho, double* rhoAxis,
        double* TAxis, double** zArray, double *dzdrho, double *dzdT);
/*
 * Evaluate the bilinear interpolant in cell (i, j) at local coordinates
 * (x, y) in [0,1]^2. No bounds checking: the caller owns the cell.
 */
double interpolateCellBilinear(double** zArray, int i, int j, double x, double y);

/*
 * Backward solve for T at fixed rho, returning the cell and local coordinates
 * so the caller can evaluate other tables in the same cell without searching
 * again. Returns 0 on success (i, j, x, y, T filled; y may lie slightly
 * outside [0,1] within the inversion tolerance), -1 on failure.
 */
int backwardInterpolateTemperatureCell(double rho, double z, int nT, int nRho, double* rhoAxis,
        double* TAxis, double** zArray, int *pi, int *pj, double *px, double *py, double *pT);
#endif
