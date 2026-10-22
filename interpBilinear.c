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

#define INTERP_EDGE_TOL 1e-4
#define INTERP_REL_TOL  1e-12

/*
 * ---------------------------------------------------------------------------
 * Backward interpolation (inversion along one table axis)
 *
 * Fixing one coordinate reduces the bilinear interpolant to a piecewise
 * linear function of the other, with node values
 *
 *   along T   (fixed rho cell i, local x):  h_j = (1-x) z[i][j] + x z[i+1][j]
 *   along rho (fixed T   cell j, local y):  g_i = (1-y) z[i][j] + y z[i][j+1]
 *
 * Inverting either is: find the cell whose two node values bracket z, then
 * solve the linear function for the local coordinate. invertLine() does this
 * for both axes; the two public functions only translate between table
 * coordinates and physical (rho, T).
 * ---------------------------------------------------------------------------
 */

enum { INTERP_ALONG_T = 0, INTERP_ALONG_RHO = 1 };

/*
 * Node value k of the line obtained by fixing the other axis at cell `fixed`
 * with local coordinate w.
 */
static inline double lineValue(double* const* zArray, int axis, int fixed, double w, int k)
{
    if (axis == INTERP_ALONG_T)
        return (1.0 - w) * zArray[fixed][k] + w * zArray[fixed+1][k];
    else
        return (1.0 - w) * zArray[k][fixed] + w * zArray[k][fixed+1];
}

/*
 * Linear scan for the lowest-index cell [k, k+1] whose node values bracket z
 * (inclusive). Returns 0 and fills k, glo, ghi on success, -1 if none does.
 */
static int scanLine(double* const* zArray, int axis, int fixed, double w, int n, double z,
        int *pk, double *pglo, double *pghi)
{
    double ghi = lineValue(zArray, axis, fixed, w, 0);
    for (int k = 0; k < n - 1; k++) {
        double glo = ghi;
        ghi = lineValue(zArray, axis, fixed, w, k+1);
        if ((glo <= z && z <= ghi) || (ghi <= z && z <= glo)) {
            *pk = k;
            *pglo = glo;
            *pghi = ghi;
            return 0;
        }
    }
    return -1;
}

/*
 * Relaxed fallback once no cell brackets z exactly:
 *   1. z overshoots the first or last node by at most edgeTol times the span
 *      of the end cell -> clamp to that node;
 *   2. z is within INTERP_REL_TOL of some node value (roundoff on a plateau,
 *      e.g. u independent of rho in the vapour region) -> that node, lowest
 *      index first.
 * Returns 0 and fills k, t on success, -1 otherwise.
 */
static int relaxedLocate(double* const* zArray, int axis, int fixed, double w, int n, double z,
        double edgeTol, int *pk, double *pt)
{
    double g0 = lineValue(zArray, axis, fixed, w, 0);
    double g1 = lineValue(zArray, axis, fixed, w, n-1);
    double gb = lineValue(zArray, axis, fixed, w, 1);
    double gt = lineValue(zArray, axis, fixed, w, n-2);
    double zbelow = g0 - edgeTol * (gb - g0);
    double zabove = g1 + edgeTol * (g1 - gt);
    if ((zbelow <= z && z <= g0) || (g0 <= z && z <= zbelow)) { *pk = 0;   *pt = 0.0; return 0; }
    if ((g1 <= z && z <= zabove) || (zabove <= z && z <= g1)) { *pk = n-2; *pt = 1.0; return 0; }
    for (int k = 0; k < n - 1; k++) {
        double gk = lineValue(zArray, axis, fixed, w, k);
        if (fabs(z - gk) <= INTERP_REL_TOL * fmax(fabs(z), fabs(gk))) { *pk = k; *pt = 0.0; return 0; }
    }
    return -1;
}

/*
 * Invert the piecewise linear line: find k and t in [0,1] with
 *
 *   z = (1-t) node[k] + t node[k+1].
 *
 * Fast path: if the end-point node values bracket z, a binary search over k
 * finds the cell, which assumes the line is monotonic. If the end points do
 * not bracket z, or the search lands in a cell that does not (non-monotonic
 * line), a linear scan returns the lowest-index bracketing cell; if no cell
 * brackets z exactly, relaxedLocate() applies the edge and roundoff
 * tolerances.
 *
 * When several roots exist the result depends on the path: the binary search
 * returns the root it bisects onto (normally the one in the monotonic bulk of
 * the line); the scan fallback returns the lowest-index one.
 *
 * Returns 0 on success, -1 if z is not attained anywhere along the line.
 */
static int invertLine(double* const* zArray, int axis, int fixed, double w, int n, double z,
        double edgeTol, int *pk, double *pt)
{
    double g0 = lineValue(zArray, axis, fixed, w, 0);
    double g1 = lineValue(zArray, axis, fixed, w, n-1);
    int lo, hi;
    double glo, ghi;
    int bracketed = 0;

    if ((g0 < z) != (g1 < z) || g0 == z || g1 == z) {
        /* Bracketed by the end points: binary search over node indices. */
        lo = 0; hi = n - 1; glo = g0; ghi = g1;
        while (hi - lo > 1) {
            int mid = lo + (hi - lo) / 2;
            double gmid = lineValue(zArray, axis, fixed, w, mid);
            if ((gmid < z) == (glo < z) && gmid != z) { lo = mid; glo = gmid; }
            else                                      { hi = mid; ghi = gmid; }
        }
        bracketed = (glo <= z && z <= ghi) || (ghi <= z && z <= glo);
    }
    if (!bracketed) {
        /* Non-monotonic line or end points not bracketing: scan. */
        if (scanLine(zArray, axis, fixed, w, n, z, &lo, &glo, &ghi) != 0)
            return relaxedLocate(zArray, axis, fixed, w, n, z, edgeTol, pk, pt);
        hi = lo + 1;
    }

    /* Solve the linear function in cell [lo, hi]. */
    double t;
    if (ghi == glo) t = 0.0;                 /* flat cell: any t works, take the lower node */
    else            t = (z - glo) / (ghi - glo);
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    *pk = lo;
    *pt = t;
    return 0;
}

/*
 * Backward interpolate the temperature: find T with z(rho, T) = z at fixed rho,
 * and report the cell (i, j) and local coordinates (x, y) of the solution so
 * the caller can evaluate other tables in the same cell without searching
 * again. See invertLine() for the search strategy and root convention.
 *
 * Returns 0 on success, -1 if rho is outside the table or z is not attained
 * at this rho.
 */
int backwardInterpolateTemperatureCell(double rho, double z, int nT, int nRho, const double* rhoAxis, const double* TAxis, double* const* zArray, int *pi, int *pj, double *px, double *py, double *pT)
{
    if (rho < rhoAxis[0] || rho >= rhoAxis[nRho-1])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS backwardInterpolateTemperature failed, rho = %.15e is outside [%.15e, %.15e)\n", rho, rhoAxis[0], rhoAxis[nRho-1]);
#endif
        return -1;
    }

    int i = findIndex(rho, rhoAxis, nRho);
    double x = (rho - rhoAxis[i]) / (rhoAxis[i+1] - rhoAxis[i]);
    *pi = i;
    *px = x;

    int j;
    double y;
    if (invertLine(zArray, INTERP_ALONG_T, i, x, nT, z, INTERP_EDGE_TOL, &j, &y) != 0)
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS backwardInterpolateTemperature failed, z = %.15e is not attained at rho = %.15e (z(Tmin) = %.15e, z(Tmax) = %.15e)\n",
                z, rho, lineValue(zArray, INTERP_ALONG_T, i, x, 0), lineValue(zArray, INTERP_ALONG_T, i, x, nT-1));
#endif
        return -1;
    }

    *pj = j;
    *py = y;
    *pT = TAxis[j] + y * (TAxis[j+1] - TAxis[j]);
    return 0;
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
 * Backward interpolate the density: find rho with z(rho, T) = z at fixed T.
 * See invertLine() for the search strategy and root convention.
 *
 * Returns rho in [rhoAxis[0], rhoAxis[nRho-1]), or -1e50 if T is outside the
 * table or z is not attained at this T.
 */
double backwardInterpolateDensityBilinear(double T, double z, int nT, int nRho, const double* rhoAxis, const double* TAxis, double* const* zArray)
{
    if (T < TAxis[0] || T >= TAxis[nT-1])
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS backwardInterpolateDensity failed, T = %.15e is outside [%.15e, %.15e)\n", T, TAxis[0], TAxis[nT-1]);
#endif
        return -1e50;
    }

    int j = findIndex(T, TAxis, nT);
    double y = (T - TAxis[j]) / (TAxis[j+1] - TAxis[j]);

    int i;
    double x;
    if (invertLine(zArray, INTERP_ALONG_RHO, j, y, nRho, z, INTERP_EDGE_TOL, &i, &x) != 0)
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS backwardInterpolateDensity failed, z = %.15e is not attained at T = %.15e (z(rhoMin) = %.15e, z(rhoMax) = %.15e)\n",
                z, T, lineValue(zArray, INTERP_ALONG_RHO, j, y, 0), lineValue(zArray, INTERP_ALONG_RHO, j, y, nRho-1));
#endif
        return -1e50;
    }

    double rho = rhoAxis[i] + x * (rhoAxis[i+1] - rhoAxis[i]);
    if (rho >= rhoAxis[nRho-1]) rho = nextafter(rhoAxis[nRho-1], 0.0);   /* stay below the exclusive upper bound */
    return rho;
}

/*
 * Nearest-neighbour lookup in an integer table (used for the phase table).
 * Returns the value at the table node closest to (rho, T), or -1 if (rho, T)
 * is outside the table or not a number.
 */
int interpolateValueNearest(double rho, double T, int nT, int nRho, const double* rhoAxis, const double* TAxis, int* const* zArray)
{
    // check if (rho,T) is out of bounds (a NaN fails every comparison below)
    if (!(rho >= rhoAxis[0]))
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS interpolateValueNearest failed, rho = %.15e is smaller than minRho = %.15e\n", rho, rhoAxis[0]);
#endif
        return -1;
    }
    if (!(rho < rhoAxis[nRho-1]))
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS interpolateValueNearest failed, rho = %.15e is larger than maxRho = %.15e\n", rho, rhoAxis[nRho-1]);
#endif
        return -1;
    }
    if (!(T >= TAxis[0]))
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS interpolateValueNearest failed, T = %.15e is smaller than minT = %.15e\n", T, TAxis[0]);
#endif
        return -1;
    }
    if (!(T < TAxis[nT-1]))
    {
#ifdef EOSLIB_VERBOSE
        fprintf(stderr,"ANEOS interpolateValueNearest failed, T = %.15e is larger than maxT = %.15e\n", T, TAxis[nT-1]);
#endif
        return -1;
    }

    // locate the cell (i, j) containing the point; tables are rho-major, zArray[i][j]
    int i=findIndex(rho, rhoAxis, nRho);
    int j=findIndex(T, TAxis, nT);

    // scaling grid rectangle to unit square
    double x=(rho-rhoAxis[i])/(rhoAxis[i+1]-rhoAxis[i]);
    double y=(T-TAxis[j])/(TAxis[j+1]-TAxis[j]);

    // pick the nearest of the four cell corners
    if (x > 0.5) i++;
    if (y > 0.5) j++;
    return zArray[i][j];
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