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
        double zlow = ylow*(f11*x - f01*(x - 1)) - (f10*x - f00*(x - 1))*(ylow - 1);
        double zhigh = yhigh*(f11*x - f01*(x - 1)) - (f10*x - f00*(x - 1))*(yhigh - 1);
        if (zhigh > zlow && z > zhigh) {
            // Normal arrangement, we have to look above
            a = c;
            c = (a + b) / 2;
        } else if (zlow >= zhigh && z < zhigh) {
            // Inverted arrangement, we have to look above
            a = c;
            c = (a + b) / 2;
        } else if (zhigh > zlow && z < zlow) {
            // Normal arrangement, we have to look below
            b = c;
            c = (a + b) / 2;
        } else if (zlow >= zhigh && z > zlow){
            // Inverted arrangement, we have to look below
            b = c;
            c = (a + b) / 2;
        } else if (z <= zhigh && z >= zlow) {
            // Normal arrangement, its inside the cell
            double y = -(z - f10*x + f00*(x - 1))/(f10*x - f11*x - f00*(x - 1) + f01*(x - 1));
            T = (TAxis[c+1]-TAxis[c])*y+TAxis[c];
            break;
        } else if (z >= zhigh && z <= zlow) {
            // Inverted arrangement, its inside the cell
            double y = -(z - f10*x + f00*(x - 1))/(f10*x - f11*x - f00*(x - 1) + f01*(x - 1));
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

        double y = -(z - f10*x + f00*(x - 1))/(f10*x - f11*x - f00*(x - 1) + f01*(x - 1));

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
        double zbottom = y*(f11*x - f01*(x - 1)) - (f10*x - f00*(x - 1))*(y - 1);

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
        double ztop = y*(f11*x - f01*(x - 1)) - (f10*x - f00*(x - 1))*(y - 1);
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

    z = y*(f11*x - f01*(x - 1)) - (f10*x - f00*(x - 1))*(y - 1);
    return z;
}

int findIndex(double x, double* xAxis, int nX)
{
    // This code assumes that nX - 1 is is a multiple of 50
    int startIndex = 0;
    int n = 50;
    for (int testind = 1; testind < n + 1; testind++)
    {
        if (xAxis[(nX-1)/n*testind]>x)
        {
            startIndex = (nX-1)/n*(testind-1);
            break;
        }
    }
    int i = 0;
    for (int k=startIndex; k<nX; k++)
    {
        if (xAxis[k]>x)
        {
            i = k-1;
            break;
        }
    }
    return i;
}
