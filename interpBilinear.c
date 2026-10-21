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
 * Exhaustive scan over the T cells of column i for the lowest-T cell that
 * brackets z. Returns 0 and fills j, y, T on success, -1 if none does.
 */
static int scanCellT(double rho, double z, int i, double x, int nT, const double* rhoAxis, const double* TAxis, double* const* zArray, int *pj, double *py, double *pT)
{
    for (int j=0; j<(nT-1); j++)
    {
        double f00=zArray[i][j];
        double f01=zArray[i][j+1];
        double f10=zArray[i+1][j];
        double f11=zArray[i+1][j+1];

        double y = -(z - f10*x + f00*(x - 1.0))/(f10*x - f11*x - f00*(x - 1.0) + f01*(x - 1.0));

        if (y>= -0.0001 && y<=1.0001)
        {
            *pj = j;
            *py = y;
            *pT = (TAxis[j+1]-TAxis[j])*y+TAxis[j];
            return 0;
        }
    }
#ifdef EOSLIB_VERBOSE
    {
        // nothing found, find out why
        // this may not always give the correct answer
        // here we assume that z is somewhat proportional to T for a given rho
        // and that the extrema lie on the grid boundaries
        double zbottom = interpolateCellBilinear(zArray, i, 0, x, 0.0);
        double ztop    = interpolateCellBilinear(zArray, i, nT-2, x, 1.0);
        if (z < zbottom)
            fprintf(stderr,"ANEOS backwardInterpolateTemperatureBilinear failed, z = %.15e is smaller than minz = %.15e\n", z, zbottom);
        if (z > ztop)
            fprintf(stderr,"ANEOS backwardInterpolateTemperatureBilinear failed, z = %.15e is bigger than maxz = %.15e\n", z, ztop);
    }
#endif
    return -1;
}

/*
 * Backward interpolate the temperature: find T with z(rho, T) = z at fixed rho,
 * and report the cell (i, j) and local coordinates (x, y) of the solution.
 *
 * Fast path: binary search over the T cells, which assumes z increases
 * monotonically with T at fixed rho (true for u and s everywhere, and for p
 * outside two-phase regions). If the search lands in a cell that does not
 * bracket z, or the cell is degenerate (z constant in T, so the inversion is
 * undefined), it falls back to the exhaustive scan, which returns the
 * lowest-T cell that brackets z.
 * The search probes zArray[i][c] for varying c, which is contiguous in the
 * rho-major layout.
 *
 * Returns 0 on success, -1 if rho is outside the table or no cell brackets z.
 */
int backwardInterpolateTemperatureCell(double rho, double z, int nT, int nRho, const double* rhoAxis, const double* TAxis, double* const* zArray, int *pi, int *pj, double *px, double *py, double *pT)
{
    if (rho < rhoAxis[0])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS backwardInterpolateTemperatureBilinear failed, rho = %.15e is smaller than minRho = %.15e\n", rho, rhoAxis[0]);
#endif
        return -1;
    }
    if (rho >= rhoAxis[nRho-1])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS backwardInterpolateTemperatureBilinear failed, rho = %.15e is larger than maxRho = %.15e\n", rho, rhoAxis[nRho-1]);
#endif
        return -1;
    }

    int i=findIndex(rho, rhoAxis, nRho);
    double x=(rho-rhoAxis[i])/(rhoAxis[i+1]-rhoAxis[i]);
    *pi = i;
    *px = x;

    int a = 0;
    int b = nT - 2;
    int c = (a + b) / 2;
    int dontexit = 1;
    while (dontexit) {
        if ((b - a) < 2) dontexit = 0;
        double f00=zArray[i][c];
        double f01=zArray[i][c+1];
        double f10=zArray[i+1][c];
        double f11=zArray[i+1][c+1];
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
            *pj = c;
            *py = y;
            *pT = (TAxis[c+1]-TAxis[c])*y+TAxis[c];
            return 0;
        } else {
            break;
        }
    }
    return scanCellT(rho, z, i, x, nT, rhoAxis, TAxis, zArray, pj, py, pT);
}

/*
 * Backward interpolate the temperature (value only). See backwardInterpolateTemperatureCell.
 */
double backwardInterpolateTemperatureBilinear(double rho, double z, int nT, int nRho, const double* rhoAxis, const double* TAxis, double* const* zArray)
{
    int i, j;
    double x, y, T;
    if (backwardInterpolateTemperatureCell(rho, z, nT, nRho, rhoAxis, TAxis, zArray, &i, &j, &x, &y, &T) != 0)
        return -1e50;
    return T;
}

/*
 * Exhaustive-scan backward interpolation (value only). Kept for callers that
 * want the slow path explicitly.
 */
double backwardInterpolateTemperatureBilinear_Scan(double rho, double z, int nT, int nRho, const double* rhoAxis, const double* TAxis, double* const* zArray)
{
    if (rho < rhoAxis[0] || rho >= rhoAxis[nRho-1]) return -1e50;
    int i=findIndex(rho, rhoAxis, nRho);
    double x=(rho-rhoAxis[i])/(rhoAxis[i+1]-rhoAxis[i]);
    int j;
    double y, T;
    if (scanCellT(rho, z, i, x, nT, rhoAxis, TAxis, zArray, &j, &y, &T) != 0) return -1e50;
    return T;
}

/*
 * Column value g_i = z(rhoAxis[i], T) for a fixed T cell j with local
 * coordinate y. The interpolant is linear in rho between consecutive g_i.
 */
static inline double densityColumnValue(double* const* zArray, int i, int j, double y)
{
    return (1.0 - y) * zArray[i][j] + y * zArray[i][j+1];
}

/*
 * Linear scan for the lowest-density cell [i, i+1] whose column values
 * bracket z. Returns 0 and fills lo, glo, ghi on success, -1 if no cell does.
 */
static int scanCellRho(double* const* zArray, int nRho, int j, double y, double z,
        int *plo, double *pglo, double *pghi)
{
    double ghi = densityColumnValue(zArray, 0, j, y);
    for (int i = 0; i < nRho - 1; i++) {
        double glo = ghi;
        ghi = densityColumnValue(zArray, i+1, j, y);
        if ((glo <= z && z <= ghi) || (ghi <= z && z <= glo)) {
            *plo = i;
            *pglo = glo;
            *pghi = ghi;
            return 0;
        }
    }
    return -1;
}

/*
 * Backward interpolate the density: find rho with z(rho, T) = z at fixed T,
 * where z is one of the tabulated quantities (P or u).
 *
 * At fixed T the interpolant is piecewise linear in rho between the column
 * values g_i = (1-y) z[i][j] + y z[i][j+1], so the inversion is exact: find
 * the cell (i, i+1) whose g values bracket z, then solve the linear function
 * for x.
 *
 * Fast path: if the end-point values bracket z, a binary search over i finds
 * the cell, which assumes g is monotonic in rho (true for P; u is not
 * monotonic for solids since the cold curve has a minimum near rho0). If the
 * end points do not bracket z, or the search lands in a cell that does not,
 * a linear scan finds the lowest-density bracketing cell. When several roots
 * exist the lowest-density one is returned; that is the library convention.
 *
 * Returns rho in [rhoAxis[0], rhoAxis[nRho-1]), or -1e50 if T is outside the
 * table or z is not attained anywhere along the row.
 */
double backwardInterpolateDensityBilinear(double T, double z, int nT, int nRho, const double* rhoAxis, const double* TAxis, double* const* zArray)
{
    if (T < TAxis[0] || T >= TAxis[nT-1])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS backwardInterpolateDensityBilinear failed, T = %.15e is outside [%.15e, %.15e)\n", T, TAxis[0], TAxis[nT-1]);
#endif
        return -1e50;
    }

    int j = findIndex(T, TAxis, nT);
    double y = (T - TAxis[j]) / (TAxis[j+1] - TAxis[j]);

    double g0 = densityColumnValue(zArray, 0, j, y);
    double g1 = densityColumnValue(zArray, nRho-1, j, y);

    /* Exact hits at the end points. */
    if (g0 == z) return rhoAxis[0];
    if (g1 == z) return nextafter(rhoAxis[nRho-1], 0.0);   /* stay below the exclusive upper bound */

    int lo, hi;
    double glo, ghi;

    if ((g0 < z) == (g1 < z)) {
        /* End points do not bracket z: g may be non-monotonic. */
        if (scanCellRho(zArray, nRho, j, y, z, &lo, &glo, &ghi) != 0) {
#ifdef EOSLIB_VERBOSE
            fprintf(stderr,"ANEOS backwardInterpolateDensityBilinear failed, z = %.15e is not attained at T = %.15e (z(rhoMin) = %.15e, z(rhoMax) = %.15e)\n", z, T, g0, g1);
#endif
            return -1e50;
        }
        hi = lo + 1;
    } else {
        /* Binary search for the bracketing cell, assuming monotonic g. */
        lo = 0; hi = nRho - 1; glo = g0; ghi = g1;
        while (hi - lo > 1) {
            int mid = lo + (hi - lo) / 2;
            double gmid = densityColumnValue(zArray, mid, j, y);
            if ((gmid < z) == (glo < z)) { lo = mid; glo = gmid; }
            else                         { hi = mid; ghi = gmid; }
        }
        /* Non-monotonic g can still fool the search: verify, else fall back. */
        if ((glo < z) == (ghi < z) && glo != z && ghi != z) {
            if (scanCellRho(zArray, nRho, j, y, z, &lo, &glo, &ghi) != 0) return -1e50;
            hi = lo + 1;
        }
    }

    /* Solve the linear function in the cell. */
    double x;
    if (ghi == glo) x = 0.0;                 /* flat cell: any x works, take the lower node */
    else            x = (z - glo) / (ghi - glo);
    if (x < 0.0) x = 0.0;
    if (x > 1.0) x = 1.0;

    double rho = rhoAxis[lo] + x * (rhoAxis[hi] - rhoAxis[lo]);
    if (rho >= rhoAxis[nRho-1]) rho = nextafter(rhoAxis[nRho-1], 0.0);
    return rho;
}

/*
 * Interpolate a value using nearest neighbor interpolation
 */
int interpolateValueNearest(double rho, double T, int nT, int nRho, const double* rhoAxis, const double* TAxis, int* const* zArray)
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

    // locate the cell (i, j) containing the point; tables are rho-major, zArray[i][j]
    int i=findIndex(rho, rhoAxis, nRho);

    int j=findIndex(T, TAxis, nT);

    // scaling grid rectangle to unit square
    double x=(rho-rhoAxis[i])/(rhoAxis[i+1]-rhoAxis[i]);
    double y=(T-TAxis[j])/(TAxis[j+1]-TAxis[j]);

    // calculating bilinear interpolation
    int f00=zArray[i][j];
    int f01=zArray[i][j+1];
    int f10=zArray[i+1][j];
    int f11=zArray[i+1][j+1];
    
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
double interpolateValueBilinear(double rho, double T, int nT, int nRho, const double* rhoAxis, const double* TAxis, double* const* zArray)
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

    // locate the cell (i, j) containing the point; tables are rho-major, zArray[i][j]
    int i=findIndex(rho, rhoAxis, nRho);

    int j=findIndex(T, TAxis, nT);

    // scaling grid rectangle to unit square
    double x=(rho-rhoAxis[i])/(rhoAxis[i+1]-rhoAxis[i]);
    double y=(T-TAxis[j])/(TAxis[j+1]-TAxis[j]);

    // calculating bilinear interpolation
    double f00=zArray[i][j];
    double f01=zArray[i][j+1];
    double f10=zArray[i+1][j];
    double f11=zArray[i+1][j+1];

    double z = -1e50;

    z = y*(f11*x - f01*(x - 1.0)) - (f10*x - f00*(x - 1.0))*(y - 1.0);
    return z;
}

/*
 * Interpolate a value using linear interpolation
 */
double interpolateValueLinear(double x, int nX, const double* xAxis, const double* yAxis)
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
int interpolateGradientBilinear(double rho, double T, int nT, int nRho, const double* rhoAxis, const double* TAxis, double* const* zArray, double *dzdrho, double *dzdT)
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

    double f00 = zArray[i][j];
    double f01 = zArray[i][j+1];
    double f10 = zArray[i+1][j];
    double f11 = zArray[i+1][j+1];

    *dzdrho = ((1.0 - y) * (f10 - f00) + y * (f11 - f01)) / drho;
    *dzdT   = ((1.0 - x) * (f01 - f00) + x * (f11 - f10)) / dT;
    return 0;
}