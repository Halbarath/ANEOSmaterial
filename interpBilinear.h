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

/*
 * All 2D tables are rho-major: zArray[i][j] = z(rhoAxis[i], TAxis[j]),
 * with the T index j contiguous in memory.
 */

/*
 * Evaluate the bilinear interpolant in cell (i, j) at local coordinates
 * (x, y) in [0,1]^2, where i is the rho index and j the T index.
 * No bounds checking: the caller owns the cell.
 */
static inline double interpolateCellBilinear(double* const* zArray, int i, int j, double x, double y)
{
    double f00=zArray[i][j];
    double f01=zArray[i][j+1];
    double f10=zArray[i+1][j];
    double f11=zArray[i+1][j+1];
    return y*(f11*x - f01*(x - 1.0)) - (f10*x - f00*(x - 1.0))*(y - 1.0);
}
/*
 * Find i such that xAxis[i] <= x < xAxis[i+1] for a strictly increasing axis.
 *
 * The result is clamped to [0, nX-2], i.e. x < xAxis[0] gives 0 and
 * x >= xAxis[nX-1] gives nX-2, so the returned cell is always valid for
 * interpolation (which then extrapolates linearly). Callers that must reject
 * out-of-range input check the bounds themselves before calling.
 *
 * O(log nX), no assumption on nX.
 */
static inline int findIndex(double x, const double* xAxis, int nX)
{
    int lo = 0;
    int hi = nX - 1;

    if (x < xAxis[0]) return 0;
    if (x >= xAxis[nX - 1]) return nX - 2;

    /* Invariant: xAxis[lo] <= x < xAxis[hi]. */
    while (hi - lo > 1) {
        int mid = lo + (hi - lo) / 2;
        if (xAxis[mid] <= x) lo = mid;
        else                 hi = mid;
    }
    return lo;
}

double backwardInterpolateTemperatureBilinear(double rho, double z, int nT, int nRho, const double* rhoAxis, const double* TAxis, double* const* zArray);
double backwardInterpolateDensityBilinear(double T, double z, int nT, int nRho, const double* rhoAxis, const double* TAxis, double* const* zArray);
double interpolateValueBilinear(double rho, double T, int nT, int nRho, const double* rhoAxis, const double* TAxis, double* const* zArray);
int interpolateValueNearest(double rho, double T, int nT, int nRho, const double* rhoAxis, const double* TAxis, int* const* zArray);
double interpolateValueLinear(double x, int nX, const double* xAxis, const double* yAxis);
int interpolateGradientBilinear(double rho, double T, int nT, int nRho, const double* rhoAxis, const double* TAxis, double* const* zArray, double *dzdrho, double *dzdT);

/*
 * Backward solve for T at fixed rho, returning the cell and local coordinates
 * so the caller can evaluate other tables in the same cell without searching
 * again. Returns 0 on success (i, j, x, y, T filled; y in [0,1]), -1 on failure.
 */
int backwardInterpolateTemperatureCell(double rho, double z, int nT, int nRho, const double* rhoAxis, const double* TAxis, double* const* zArray, int *pi, int *pj, double *px, double *py, double *pT);
#endif
