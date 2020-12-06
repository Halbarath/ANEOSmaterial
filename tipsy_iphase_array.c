/*
 * This program reads a Tipsy binary file and generates a tipsy array file that contains the phase
 * of each particle calculated from ANEOS.
 *
 * Author:   Christian Reinhardt
 * Created:  05.12.2020
 * Modified:  
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tipsy.h"
#include "aneos.h"
#include "ANEOSmaterial.h"


double aneosUofRhoT(int iMat, double rho, double T) {
    double p;
    double u;
    double s;
    double cv;
    double dPdT;
    double dPdrho;
    double fkros;
    double cs;
    int iPhase;
    double rhoL;
    double rhoH;
    double ion;

    callaneos_cgs(T, rho, iMat, &p, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &cs, &iPhase, &rhoL, &rhoH,
                  &ion);

    return u;
}

int aneosPhase(int iMat, double rho, double T) {
    double p;
    double u;
    double s;
    double cv;
    double dPdT;
    double dPdrho;
    double fkros;
    double cs;
    int iPhase;
    double rhoL;
    double rhoH;
    double ion;

    callaneos_cgs(T, rho, iMat, &p, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &cs, &iPhase, &rhoL, &rhoH,
                  &ion);

    return iPhase;
}

double aneosTofRhoU(int iMat, double rho, double u) {
    double T_min = 1.0;
    double T_max = 1e6;
    double T_mean = 0.0;
    double u_min, u_max, u_mean;

    u_min = aneosUofRhoT(iMat, rho, T_min);
    u_max = aneosUofRhoT(iMat, rho, T_max);
    u_mean = 0.0;

    /* Fail if the root is not bracketed. */
    if ((u_min > u) || (u_max < u)) {
        fprintf(stderr, "Could not bracket root.\n");
        fprintf(stderr, "rho = %g u_min= %g u_max= %g u= %g\n", rho, u_min, u_max, u);
        assert(0);
    }

    if (u_min >= u) return T_min;
    if (u_max <= u) return T_max;

    /* Bisection. */
    while (2*(T_max-T_min)/(T_max+T_min) > 1e-12) {
        T_mean = 0.5*(T_min + T_max);
        u_mean = aneosUofRhoT(iMat, rho, T_mean);
        
        if (u_mean < u) {
            T_min = T_mean;
            u_min = u_mean;
        } else {
            T_max = T_mean;
            u_max = u_mean;
        }
    }

    return T_mean;
}

int TestBisection(int iMat) {
    double rho_min = 1e-4;
    double rho_max = 1e2;
    double T_min = 1.0;
    double T_max = 1e6;
    int nRho = 200;
    int nT= 100;
    double rho, T, u;
    FILE *fp;
    int i, j;

    fp = fopen("out.txt", "w");
    assert(fp != NULL);

    for (i=0; i<nT; i++) {
        T = T_min + i*(T_max-T_min)/(nT-1);
        for (j=0; j<nRho; j++) {
            rho = rho_min + j*(rho_max-rho_min)/(nRho-1);
            
            printf("T=%g rho= %g\n", T, rho);
            u = aneosUofRhoT(iMat, rho, T);

            fprintf(fp, "%15.7E", (T-aneosTofRhoU(iMat, rho, u))/T);
        }
        fprintf(fp, "\n");
    }
    
    fclose(fp); 
    return 0;
}

int main(int argc, char **argv) {
    /* Tipsy library. */
    TCTX in;
    int type, N;
    double dSoft;
    double time;
    /* ANEOS library. */
    char matFilename[256] = "aneos.input";
    // Uncomment below to use M-ANEOS
    //char matFilename[256] = "maneos.in";
    double p;
    double u;
    double s;
    double cv;
    double dPdT;
    double dPdrho;
    double fkros;
    double cs;
    int iPhase;
    double rhoL;
    double rhoH;
    double ion;
    /* Units. */ 
    const double KBOLTZ = 1.38e-16;      /* bolzman constant in cgs */
    const double MHYDR = 1.67e-24;       /* mass of hydrogen atom in grams */
    const double MSOLG = 1.99e33;        /* solar mass in grams */
    const double GCGS = 6.67e-8;         /* G in cgs */
    const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */
    const double NA = 6.022e23;          /* Avogadro's number */
    double dErgPerGmUnit;
    double dGmPerCcUnit;
    double dSecUnit;	
    double dKpcUnit = 2.06701e-13;
    double dMsolUnit = 4.80438e-08;
    double rho, T;
    int iMat;
    int i;

#if 0
    iMat = 4;
    initaneos(matFilename);
    TestBisection(iMat);
    
    exit(1);
#endif

    if (argc != 1) {
        fprintf(stderr,"Usage: tipsy_iphase_array <tipsy.std\n");
        exit(1);
    }

    /* Units. */
	dErgPerGmUnit = GCGS*dMsolUnit*MSOLG/(dKpcUnit*KPCCM);
	dGmPerCcUnit = (dMsolUnit*MSOLG)/pow(dKpcUnit*KPCCM, 3.0);
    
    fprintf(stderr, "dErgPerGmUnit= %g\n", dErgPerGmUnit);
    fprintf(stderr, "dGmPerCcUnit = %g\n", dGmPerCcUnit);
    	
    /* Initialize ANEOS. */
    fprintf(stderr, "ANEOS: Initializing material...\n");
    initaneos(matFilename);

    /* Initialize tipsy library */
    TipsyInitialize(&in, 0, "stdin");

    N = iTipsyNumParticles(in);
    time = dTipsyTime(in);

    /* Read all particles. */
	TipsyReadAll(in);

    /* Print header for the array file. */
    printf("%i\n",N);

    for (i = 0; i < N; i++) {
        iMat = (int) in->gp[i].metals;
 
        /* Convert iMat to the corresponding material id in the aneos input file. */
        switch(iMat) {
            case 51:
                iMat = 1;
                break;
            case 52:
                iMat = 2;
                break;
            case 54:
                iMat = 4;
                break;
            case 55:
                iMat = 5;
                break;
            default:
                fprintf(stderr, "Unknown material id: %i\n", iMat);
                exit(1);
        }
        
        rho = in->gp[i].rho;
        u = in->gp[i].temp;

        /* Convert to cgs. */
        rho *= dGmPerCcUnit;
        u *= dErgPerGmUnit;
         
        T = aneosTofRhoU(iMat, rho, u);

        callaneos_cgs(T, rho, iMat, &p, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &cs, &iPhase, &rhoL,
                      &rhoH, &ion);

        printf("%i\n", iPhase);
    }

    TipsyFinish(in);

    return 0;
}
