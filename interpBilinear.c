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
 * ANEOS material bilinear interpolation functions
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "interpBilinear.h"

/*
 * Backward interpolate the temperature using bilinear interpolation
 */
double backwardInterpolateTemperatureBilinear(double rho, double z, int nT, int nRho, double* rhoAxis,
        double* TAxis, double** zArray)
{
    // check if (rho,T) is out of bounds
    if (rho < rhoAxis[0])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS backwardInterpolateTemperatureBilinear failed, rho = %.15e is smaller than minRho = %.15e\n", rho, rhoAxis[0]);
#endif
        return -1e50;
    }
    if (rho >= rhoAxis[nRho-1])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS backwardInterpolateTemperatureBilinear failed, rho = %.15e is larger than maxRho = %.15e\n", rho, rhoAxis[nRho-1]);
#endif
        return -1e50;
    }

    int i=findIndex(rho, rhoAxis, nRho);
    double x=(rho-rhoAxis[i])/(rhoAxis[i+1]-rhoAxis[i]);
    double T = -1e50;
    int a = 0;
    int b = nT - 2;
    int c = (a + b) / 2;
    int dontexit = 1;
    while (dontexit) {
        if ((b - a) < 2) dontexit = 0;
        double f00=zArray[c][i];
        double f01=zArray[c+1][i];
        double f10=zArray[c][i+1];
        double f11=zArray[c+1][i+1];
        double ylow = -0.0001;
        double yhigh = 1.0001;
        double zlow = ylow*(f11*x - f01*(x - 1.0)) - (f10*x - f00*(x - 1.0))*(ylow - 1.0);
        double zhigh = yhigh*(f11*x - f01*(x - 1.0)) - (f10*x - f00*(x - 1.0))*(yhigh - 1.0);
        if (z > zhigh) {
            a = c;
            c = (a + b) / 2;
        } else if (z < zlow) {
            b = c;
            c = (a + b) / 2;
        } else if ((zlow < zhigh && (z <= zhigh && z >= zlow)) || (zhigh <= zlow && (z >= zhigh && z <= zlow))) {
            double y = -(z - f10*x + f00*(x - 1.0))/(f10*x - f11*x - f00*(x - 1.0) + f01*(x - 1.0));
            T = (TAxis[c+1]-TAxis[c])*y+TAxis[c];
            break;
        } else {
            break;
        }
    }
    if (T < -1e40) {
        T = backwardInterpolateTemperatureBilinearOld(rho, z, nT, nRho, rhoAxis, TAxis, zArray);
    }
    return T;
}

/*
 * Backward interpolate the temperature using bilinear interpolation
 */
double backwardInterpolateTemperatureBilinearOld(double rho, double z, int nT, int nRho, double* rhoAxis,
        double* TAxis, double** zArray)
{
    // check if (rho,T) is out of bounds
    if (rho < rhoAxis[0])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS backwardInterpolateTemperatureBilinear failed, rho = %.15e is smaller than minRho = %.15e\n", rho, rhoAxis[0]);
#endif
        return -1e50;
    }
    if (rho >= rhoAxis[nRho-1])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS backwardInterpolateTemperatureBilinear failed, rho = %.15e is larger than maxRho = %.15e\n", rho, rhoAxis[nRho-1]);
#endif
        return -1e50;
    }

    // searching the rho interval containing the rho value
    int i=findIndex(rho, rhoAxis, nRho);

    // calculating the inverted bilinear interpolation for each of the selected rectangles
    // until the value is found
    double T = -1e50;

    for (int j=0; j<(nT-1); j++)
    {
        double x=(rho-rhoAxis[i])/(rhoAxis[i+1]-rhoAxis[i]);

        double f00=zArray[j][i];
        double f01=zArray[j+1][i];
        double f10=zArray[j][i+1];
        double f11=zArray[j+1][i+1];

        double y = -(z - f10*x + f00*(x - 1.0))/(f10*x - f11*x - f00*(x - 1.0) + f01*(x - 1.0));

        if (y>= -0.0001 && y<=1.0001)
        {
            T = (TAxis[j+1]-TAxis[j])*y+TAxis[j];
            break;
        }

    }
#ifdef EOSLIB_VERBOSE
    if (T < -1e40)
    {
        // nothing found, find out why
        // this may not always give the correct answer
        // here we assume that z is somewhat proportional to T for a given rho
        // and that the extrema lie on the grid boundaries
        int j = 0;
        double x=(rho-rhoAxis[i])/(rhoAxis[i+1]-rhoAxis[i]);
        double y = 0.0;
        double f00=zArray[j][i];
        double f01=zArray[j+1][i];
        double f10=zArray[j][i+1];
        double f11=zArray[j+1][i+1];
        double zbottom = y*(f11*x - f01*(x - 1.0)) - (f10*x - f00*(x - 1.0))*(y - 1.0);

        if (z < zbottom)
        {
            // below the grid
            fprintf(stderr,"ANEOS backwardInterpolateTemperatureBilinear failed, z = %.15e is smaller than minz = %.15e\n", z, zbottom);
        }
        j = nT-2;
        y = 1.0;
        f00=zArray[j][i];
        f01=zArray[j+1][i];
        f10=zArray[j][i+1];
        f11=zArray[j+1][i+1];
        double ztop = y*(f11*x - f01*(x - 1.0)) - (f10*x - f00*(x - 1.0))*(y - 1.0);
        if (z > ztop)
        {
            // above the grid
            fprintf(stderr,"ANEOS backwardInterpolateTemperatureBilinear failed, z = %.15e is bigger than maxz = %.15e\n", z, ztop);
        }
    }
#endif
    return T;
}

/*
 * Backward interpolate the density using bilinear interpolation
 */
double backwardInterpolateDensityBilinear(double T, double z, int nT, int nRho, double* rhoAxis,
        double* TAxis, double** zArray)
{
    // check if T is out of bounds
    if (T < TAxis[0])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS backwardInterpolateDensityBilinear failed, T = %.15e is smaller than minT = %.15e\n", T, TAxis[0]);
#endif
        return -1e50;
    }
    if (T >= TAxis[nT-1])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS backwardInterpolateDensityBilinear failed, T = %.15e is larger than maxT = %.15e\n", T, TAxis[nT-1]);
#endif
        return -1e50;
    }

    double a = rhoAxis[0];
    double b = rhoAxis[nRho-1]*0.999;
    double c = -1e50;
    double zc = 0;

    while (2*(b-a)/(b+a) > 1e-10) {
        c = 0.5*(a + b);
        zc = interpolateValueBilinear(c, T, nT, nRho, rhoAxis, TAxis, zArray);
        if (zc > z) {
            b = c;
        } else {
            a = c;
        }
    }
    return c;
}

/*
 * Interpolate a value using nearest neighbor interpolation
 */
int interpolateValueNearest(double rho, double T, int nT, int nRho, double* rhoAxis,
        double* TAxis, int** zArray)
{
    // check if (rho,T) is out of bounds
    if (rho < rhoAxis[0])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS interpolateValueBilinear failed, rho = %.15e is smaller than minRho = %.15e\n", rho, rhoAxis[0]);
#endif
        return -1;
    }
    if (rho >= rhoAxis[nRho-1])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS interpolateValueBilinear failed, rho = %.15e is larger than maxRho = %.15e\n", rho, rhoAxis[nRho-1]);
#endif
        return -1;
    }
    if (T < TAxis[0])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS interpolateValueBilinear failed, T = %.15e is smaller than minT = %.15e\n", T, TAxis[0]);
#endif
        return -1;
    }
    if (T >= TAxis[nT-1])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS interpolateValueBilinear failed, T = %.15e is larger than maxT = %.15e\n", T, TAxis[nT-1]);
#endif
        return -1;
    }

    // searching for grid rectangle containing the point
    // could be calculated if grid is guarantied to be logarithmic
    // as we do not assume that, we search for the grid rectangle
    int i=findIndex(rho, rhoAxis, nRho);

    int j=findIndex(T, TAxis, nT);

    // scaling grid rectangle to unit square
    double x=(rho-rhoAxis[i])/(rhoAxis[i+1]-rhoAxis[i]);
    double y=(T-TAxis[j])/(TAxis[j+1]-TAxis[j]);

    // calculating bilinear interpolation
    int f00=zArray[j][i];
    int f01=zArray[j+1][i];
    int f10=zArray[j][i+1];
    int f11=zArray[j+1][i+1];
    
    if ((x <= 0.5) && (y <= 0.5)) {
        return f00;
    } else if ((x <= 0.5) && (y > 0.5)) {
        return f01;
    } else if ((x > 0.5) && (y <= 0.5)) {
        return f10;
    } else if ((x > 0.5) && (y > 0.5)) {
        return f11;
    } else {
        //return f11;
        fprintf(stderr, "x= %g y= %g\n", x, y);
        exit(1);
    }
}


/*
 * Interpolate a value using bilinear interpolation
 */
double interpolateValueBilinear(double rho, double T, int nT, int nRho, double* rhoAxis,
        double* TAxis, double** zArray)
{
    // check if (rho,T) is out of bounds
    if (rho < rhoAxis[0])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS interpolateValueBilinear failed, rho = %.15e is smaller than minRho = %.15e\n", rho, rhoAxis[0]);
#endif
        return -1e50;
    }
    if (rho >= rhoAxis[nRho-1])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS interpolateValueBilinear failed, rho = %.15e is larger than maxRho = %.15e\n", rho, rhoAxis[nRho-1]);
#endif
        return -1e50;
    }
    if (T < TAxis[0])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS interpolateValueBilinear failed, T = %.15e is smaller than minT = %.15e\n", T, TAxis[0]);
#endif
        return -1e50;
    }
    if (T >= TAxis[nT-1])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS interpolateValueBilinear failed, T = %.15e is larger than maxT = %.15e\n", T, TAxis[nT-1]);
#endif
        return -1e50;
    }

    // searching for grid rectangle containing the point
    // could be calculated if grid is guarantied to be logarithmic
    // as we do not assume that, we search for the grid rectangle
    int i=findIndex(rho, rhoAxis, nRho);

    int j=findIndex(T, TAxis, nT);

    // scaling grid rectangle to unit square
    double x=(rho-rhoAxis[i])/(rhoAxis[i+1]-rhoAxis[i]);
    double y=(T-TAxis[j])/(TAxis[j+1]-TAxis[j]);

    // calculating bilinear interpolation
    double f00=zArray[j][i];
    double f01=zArray[j+1][i];
    double f10=zArray[j][i+1];
    double f11=zArray[j+1][i+1];

    double z = -1e50;

    z = y*(f11*x - f01*(x - 1.0)) - (f10*x - f00*(x - 1.0))*(y - 1.0);
    return z;
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
int findIndex(double x, double* xAxis, int nX)
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

/*
 * Interpolate a value using linear interpolation
 */
double interpolateValueLinear(double x, int nX, double* xAxis, double* yAxis)
{
    // check if x is out of bounds
    if (x < xAxis[0])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS interpolateValueLinear failed, x = %.15e is smaller than minx = %.15e\n", x, xAxis[0]);
#endif
        return -1e50;
    }
    if (x >= xAxis[nX-1])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS interpolateValueLinear failed, x = %.15e is larger than maxx = %.15e\n", x, xAxis[nX-1]);
#endif
        return -1e50;
    }

    // searching for grid rectangle containing the point
    // could be calculated if grid is guarantied to be logarithmic
    // as we do not assume that, we search for the grid rectangle
    int i=findIndex(x, xAxis, nX);

    // scaling grid rectangle to unit square
    x = (x-xAxis[i])/(xAxis[i+1]-xAxis[i]);

    // calculating linear interpolation
    double f0=yAxis[i];
    double f1=yAxis[i+1];

    double y = -1e50;

    y = f0*(1.0 - x) + f1*x;
    return y;
}

/*
 * Partial derivatives of the bilinear interpolant at (rho, T):
 *
 *   z(x,y) = (1-x)(1-y) f00 + x(1-y) f10 + (1-x) y f01 + x y f11
 *
 * with x = (rho - rho_i)/(rho_{i+1} - rho_i), y = (T - T_j)/(T_{j+1} - T_j), so
 *
 *   dz/drho = [ (1-y)(f10 - f00) + y (f11 - f01) ] / (rho_{i+1} - rho_i)
 *   dz/dT   = [ (1-x)(f01 - f00) + x (f11 - f10) ] / (T_{j+1} - T_j)
 *
 * These are exact for the interpolant and cost one cell lookup. The
 * interpolant is only C0, so on a cell boundary the derivative is the one of
 * the cell that findIndex selects (the cell to the right/above).
 *
 * Returns 0 on success, -1 if (rho, T) is outside the table (outputs untouched).
 */
int interpolateGradientBilinear(double rho, double T, int nT, int nRho, double* rhoAxis,
        double* TAxis, double** zArray, double *dzdrho, double *dzdT)
{
    if (rho < rhoAxis[0] || rho >= rhoAxis[nRho-1] || T < TAxis[0] || T >= TAxis[nT-1])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS interpolateGradientBilinear failed, rho = %.15e, T = %.15e out of bounds\n", rho, T);
#endif
        return -1;
    }

    int i = findIndex(rho, rhoAxis, nRho);
    int j = findIndex(T, TAxis, nT);

    double drho = rhoAxis[i+1] - rhoAxis[i];
    double dT   = TAxis[j+1] - TAxis[j];
    double x = (rho - rhoAxis[i]) / drho;
    double y = (T - TAxis[j]) / dT;

    double f00 = zArray[j][i];
    double f01 = zArray[j+1][i];
    double f10 = zArray[j][i+1];
    double f11 = zArray[j+1][i+1];

    *dzdrho = ((1.0 - y) * (f10 - f00) + y * (f11 - f01)) / drho;
    *dzdT   = ((1.0 - x) * (f01 - f00) + x * (f11 - f10)) / dT;
    return 0;
}