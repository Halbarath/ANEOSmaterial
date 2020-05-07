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

    // searching the rho interval containing the rho value
    int i=findIndex(rho, rhoAxis, nRho);

    // calculating the inverted bilinear interpolation for each of the selected rectangles
    // until the value is found
    double T = -1e50;

    for (int j=0; j<(nT-1); j++)
    {
        //if (indices[j]==0)
        //{
        //	continue;
        //}
        double x=(rho-rhoAxis[i])/(rhoAxis[i+1]-rhoAxis[i]);

        double f00=zArray[j][i];
        double f01=zArray[j+1][i];
        double f10=zArray[j][i+1];
        double f11=zArray[j+1][i+1];

        double y = -(z - f10*x + f00*(x - 1))/(f10*x - f11*x - f00*(x - 1) + f01*(x - 1));

        if (y>= 0.0 && y<=1.0)
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
    //free(indices);
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

    /*// searching the T interval containing the T value
      int i=0;
      for (int k=0; k<nT; k++)
      {
      if (TAxis[k]>T)
      {
      i = k-1;
      break;
      }
      }

    // selecting the rectangles that may contain the correct value
    int *indices = (int *)malloc((nRho-1) * sizeof(int));
    int qq =0;
    for (int k=0; k<nRho-1; k++)
    {
    // calculate minimum of the four points
    double minimum = fmin(fmin(fmin(zArray[i][k],zArray[i+1][k]),zArray[i][k+1]),zArray[i+1][k+1]);
    // calculate maximum of the four points
    double maximum = fmax(fmax(fmax(zArray[i][k],zArray[i+1][k]),zArray[i][k+1]),zArray[i+1][k+1]);

    if (z <= maximum && z >= minimum)
    {
    qq=qq+1;
    indices[k] = 1;
    } else {
    indices[k] = 0;
    }
    }

    // calculating the inverted bilinear interpolation for each of the selected rectangles
    // until the value is found
    double rho = -1e50;

    for (int j=0; j<(nRho-1); j++)
    {
    if (indices[j]==0)
    {
    continue;
    }
    double y=(T-TAxis[i])/(TAxis[i+1]-TAxis[i]);

    double f00=zArray[i][j];
    double f01=zArray[i+1][j];
    double f10=zArray[i][j+1];
    double f11=zArray[i+1][j+1];

    double x = -(z - f01*y + f00*(y - 1))/(y*(f01 - f11) - (f00 - f10)*(y - 1));

    if (x>= 0 && x<=1)
    {
    rho = (rhoAxis[j+1]-rhoAxis[j])*x+rhoAxis[j];
    break;
    }
    }

    if (rho < -1e40)
    {
    // nothing found, find out why
    // this may not always give the correct answer
    // here we assume that z is somewhat proportional to rho for a given T
    // and that the extrema lie on the grid boundaries
    int j = 0;
    double y=(T-TAxis[i])/(TAxis[i+1]-TAxis[i]);
    double x = 0.0;
    double f00=zArray[i][j];
    double f01=zArray[i+1][j];
    double f10=zArray[i][j+1];
    double f11=zArray[i+1][j+1];
    double zleft = y*(f11*x - f01*(x - 1)) - (f10*x - f00*(x - 1))*(y - 1);
    if (z < zleft)
    {
        // left of the grid
        fprintf(stderr,"ANEOS backwardInterpolateDensityBilinear failed, z = %.15e is smaller than minz = %.15e\n", z, zleft);
    }
    j = nT-2;
    y = 1.0;
    f00=zArray[i][j];
    f01=zArray[i+1][j];
    f10=zArray[i][j+1];
    f11=zArray[i+1][j+1];
    double zright = y*(f11*x - f01*(x - 1)) - (f10*x - f00*(x - 1))*(y - 1);
    if (z > zright)
    {
        // right of the grid
        fprintf(stderr,"ANEOS backwardInterpolateDensityBilinear failed, z = %.15e is bigger than maxz = %.15e\n", z, zright);
    }
}

free(indices);
return rho;*/
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
    // This code assumes that nX - 1 is is a multiple of 10
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
