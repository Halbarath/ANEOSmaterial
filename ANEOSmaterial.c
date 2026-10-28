/*  This file is part of ANEOSmaterial.
 *  Copyright (c) 2020-2026 Thomas Meier & Christian Reinhardt
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
 * ANEOS material model
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "ANEOSmaterial.h"

/*
 * The core's inversion policy for ANEOS tables: tolerant mode, which copes with
 * the non-monotonic u(T) columns of the (M-)ANEOS tables (see
 * INTERPINVERTOPT in interpBilinear.h). The tolerances are the ones the
 * library has always used.
 */
static const INTERPINVERTOPT ANEOSInvertOpt = INTERP_INVERT_TOLERANT;

/*
 * Initializes ANEOSMATERIAL given the material number iMat
 */
ANEOSMATERIAL *ANEOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit)
{
    char inputfile[256];

    switch(iMat)
    {
        case MAT_ANEOS_ICE:
            strcpy(inputfile, "ANEOStable_ice.in");
            break;
        case MAT_ANEOS_DUNITE:
            strcpy(inputfile, "ANEOStable_dunite.in");
            break;
        case MAT_ANEOS_IRON:
            strcpy(inputfile, "ANEOStable_iron.in");
            break;
        case MAT_MANEOS_QUARTZ:
            strcpy(inputfile, "MANEOStable_quartz.in");
            break;
        case MAT_MANEOS_FORSTERITE:
            strcpy(inputfile, "MANEOStable_forsterite.in");
            break;
        case MAT_MANEOS_IRON:
            strcpy(inputfile, "MANEOStable_iron.in");
            break;
        case MAT_MANEOS_FE85SI15:
            strcpy(inputfile, "MANEOStable_iron_Fe85Si15.in");
            break;
        case MAT_MANEOS_SERPENTINE:
            strcpy(inputfile, "MANEOStable_serpentine.in");
            break;
        case MAT_MANEOS_WATER:
            strcpy(inputfile, "MANEOStable_water.in");
            break;
        default:
            fprintf(stderr, "ANEOSinitMaterial: unknown ANEOS material %d.\n", iMat);
            assert(0);
            exit(1);
    }

    ANEOSMATERIAL *material;

    material = ANEOSinitMaterialFromFile(iMat, inputfile, dKpcUnit, dMsolUnit);

    return material;
}

/*
 * Read exactly n items or abort. A short read of a binary table is always a
 * corrupt or truncated file.
 */
static void ANEOSReadOrDie(void *ptr, size_t size, size_t n, FILE *file, const char *what, const char *inputfile)
{
    size_t nRead = fread(ptr, size, n, file);
    if (nRead != n) {
        fprintf(stderr, "ANEOSmaterial: failed to read %s from file %s (got %zu of %zu items)\n",
                what, inputfile, nRead, n);
        assert(0);
        exit(1);
    }
}

/*
 * Validation of a table at load (fatal, like a short read): the interpBilinear core needs
 * strictly increasing axes, a NaN or inf anywhere in a table would be interpolated into the
 * neighbouring cells, and rho0, rho and T must be positive (the code <-> cgs conversion and the
 * T_min fallback of T(rho, u)).
 */
static void ANEOSCheckAxis(const double *axis, int n, const char *what, const char *inputfile)
{
    for (int i = 0; i < n; i++) {
        if (!isfinite(axis[i]) || !(axis[i] > 0.0)) {
            fprintf(stderr, "ANEOSmaterial: %s[%d] = %.15e in %s is not positive and finite.\n", what, i, axis[i], inputfile);
            assert(0);
            exit(1);
        }
        if (i > 0 && !(axis[i] > axis[i-1])) {
            fprintf(stderr, "ANEOSmaterial: %s in %s is not strictly increasing at %d (%.15e <= %.15e).\n",
                    what, inputfile, i, axis[i], axis[i-1]);
            assert(0);
            exit(1);
        }
    }
}

static void ANEOSCheckFinite(const double *data, int nRho, int nT, const char *what, const char *inputfile)
{
    for (size_t k = 0; k < (size_t)nRho * nT; k++) {
        if (!isfinite(data[k])) {
            fprintf(stderr, "ANEOSmaterial: %s in %s has a non-finite entry %g at iRho= %zu iT= %zu.\n",
                    what, inputfile, data[k], k / nT, k % nT);
            assert(0);
            exit(1);
        }
    }
}

/*
 * Allocate a rho-major table: one contiguous block of nRho*nT values plus
 * nRho column pointers into it. Returns the column pointers; the block is
 * returned through *pData so it can be freed.
 */
static double **ANEOSAllocTable(int nRho, int nT, double **pData)
{
    double *data = (double *)malloc((size_t)nRho * nT * sizeof(double));
    assert(data != NULL);
    double **cols = (double **)malloc((size_t)nRho * sizeof(double *));
    assert(cols != NULL);
    for (int i = 0; i < nRho; i++) cols[i] = data + (size_t)i * nT;
    *pData = data;
    return cols;
}

/*
 * Read one table from the file, which stores it T-major (nT rows of nRho
 * values), and transpose it into the rho-major layout. tmp must hold
 * nRho*nT doubles.
 */
static void ANEOSReadTableTransposed(double **cols, int nRho, int nT, double *tmp,
                                     FILE *file, const char *what, const char *inputfile)
{
    ANEOSReadOrDie(tmp, sizeof(double), (size_t)nRho * nT, file, what, inputfile);
    for (int j = 0; j < nT; j++) {
        const double *row = tmp + (size_t)j * nRho;
        for (int i = 0; i < nRho; i++) cols[i][j] = row[i];
    }
}

ANEOSMATERIAL *ANEOSinitMaterialFromFile(int iMat, char *inputfile, double dKpcUnit, double dMsolUnit)
{
    double rho0;
    int nRho;
    int nT;

    FILE *file;

    /*
     * Unit system (the same rule in EOSlib and every sublibrary): both units > 0 and finite are code
     * units, both <= 0 is cgs; anything else (exactly one of them <= 0, NaN, inf) is a caller error.
     */
    if (!(((dKpcUnit > 0.0) && (dMsolUnit > 0.0) && isfinite(dKpcUnit) && isfinite(dMsolUnit)) || ((dKpcUnit <= 0.0) && (dMsolUnit <= 0.0)))) {
        fprintf(stderr, "ANEOSinitMaterialFromFile: invalid unit system dKpcUnit= %g dMsolUnit= %g (give both > 0 for code units or both <= 0 for cgs).\n", dKpcUnit, dMsolUnit);
        assert(0);
        exit(1);
    }

    if ((file= fopen(inputfile, "rb")) == NULL)
    {
        fprintf(stderr, "Could not open file %s\n", inputfile);
        assert(0);
        exit(1);
    }
    ANEOSReadOrDie(&rho0, sizeof(rho0), 1, file, "rho0", inputfile);
    ANEOSReadOrDie(&nRho, sizeof(nRho), 1, file, "nRho", inputfile);
    ANEOSReadOrDie(&nT,   sizeof(nT),   1, file, "nT",   inputfile);

    if (nRho < 2 || nT < 2) {
        fprintf(stderr, "ANEOSinitMaterialFromFile: invalid table size nRho = %d, nT = %d in %s\n", nRho, nT, inputfile);
        assert(0);
        exit(1);
    }
    if (!isfinite(rho0) || !(rho0 > 0.0)) {
        fprintf(stderr, "ANEOSinitMaterialFromFile: reference density rho0 = %g in %s is not positive and finite.\n", rho0, inputfile);
        assert(0);
        exit(1);
    }

    ANEOSMATERIAL *material;
    material = (ANEOSMATERIAL *) calloc(1, sizeof(ANEOSMATERIAL));
    assert(material != NULL);
    material->iMat = iMat;
    material->rho0 = rho0;
    material->nRho = nRho;
    material->nT = nT;
    snprintf(material->matName, sizeof(material->matName), "%s", inputfile);

    material->CodeUnitstoCGSforU = 1.0;
    material->CodeUnitstoCGSforP = 1.0;
    material->CodeUnitstoCGSforRho = 1.0;
    material->CodeUnitstoCGSforC = 1.0;

    if (dKpcUnit > 0.0 && dMsolUnit > 0.0) {
        /* Gasoline/pkdgrav3 unit constants. Must match the hydro code and the other EOS libraries exactly. */
        const double MSOLG = 1.99e33;        /* solar mass in grams */
        const double GCGS = 6.67e-8;         /* G in cgs */
        const double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */

        // unit of u is erg/g = cm^2/s^2
        material->CodeUnitstoCGSforU = GCGS*dMsolUnit*MSOLG/(dKpcUnit*KPCCM);
        // unit of rho is g/cm^3
        material->CodeUnitstoCGSforRho = (dMsolUnit*MSOLG)/pow(dKpcUnit*KPCCM,3.0);
        // unit of p is erg/cm^3=1/(cm*s^2)
        // combine the upper two
        material->CodeUnitstoCGSforP = material->CodeUnitstoCGSforU*material->CodeUnitstoCGSforRho;
        // unit of c is cm/s
        material->CodeUnitstoCGSforC = dKpcUnit*KPCCM*sqrt((material->CodeUnitstoCGSforRho*GCGS));
    }

    material->CGStoCodeUnitsforU = 1.0 / material->CodeUnitstoCGSforU;
    material->CGStoCodeUnitsforP = 1.0 / material->CodeUnitstoCGSforP;
    material->CGStoCodeUnitsforRho = 1.0 / material->CodeUnitstoCGSforRho;
    material->CGStoCodeUnitsforC = 1.0 / material->CodeUnitstoCGSforC;

    double *rhoAxis = (double *)malloc(sizeof(double)*nRho);
    assert(rhoAxis != NULL);
    double *TAxis = (double *)malloc(sizeof(double)*nT);
    assert(TAxis != NULL);

    material->rhoAxis = rhoAxis;
    material->TAxis = TAxis;
    material->grid.nx = nRho;
    material->grid.ny = nT;
    material->grid.x = rhoAxis;
    material->grid.y = TAxis;
    material->uArray = ANEOSAllocTable(nRho, nT, &material->uData);
    material->pArray = ANEOSAllocTable(nRho, nT, &material->pData);
    material->cArray = ANEOSAllocTable(nRho, nT, &material->cData);
    material->sArray = ANEOSAllocTable(nRho, nT, &material->sData);

    ANEOSReadOrDie(rhoAxis, sizeof(rhoAxis[0]), nRho, file, "rhoAxis", inputfile);
    ANEOSReadOrDie(TAxis,   sizeof(TAxis[0]),   nT,   file, "TAxis",   inputfile);

    /* Tables are stored T-major in the file; read each with one fread and transpose. */
    double *tmp = (double *)malloc((size_t)nRho * nT * sizeof(double));
    assert(tmp != NULL);
    ANEOSReadTableTransposed(material->pArray, nRho, nT, tmp, file, "pArray", inputfile);
    ANEOSReadTableTransposed(material->uArray, nRho, nT, tmp, file, "uArray", inputfile);
    ANEOSReadTableTransposed(material->sArray, nRho, nT, tmp, file, "sArray", inputfile);
    ANEOSReadTableTransposed(material->cArray, nRho, nT, tmp, file, "cArray", inputfile);
    free(tmp);

    /* Validate before anything is evaluated (see ANEOSCheckAxis). */
    ANEOSCheckAxis(rhoAxis, nRho, "rhoAxis", inputfile);
    ANEOSCheckAxis(TAxis, nT, "TAxis", inputfile);
    ANEOSCheckFinite(material->pData, nRho, nT, "pArray", inputfile);
    ANEOSCheckFinite(material->uData, nRho, nT, "uArray", inputfile);
    ANEOSCheckFinite(material->sData, nRho, nT, "sArray", inputfile);
    ANEOSCheckFinite(material->cData, nRho, nT, "cArray", inputfile);

    /* The material string is optional: older tables end here. */
    size_t nStr = fread(material->matstring, 1, sizeof(material->matstring) - 1, file);
    if (nStr == 0) {
        snprintf(material->matstring, sizeof(material->matstring), "No material string given");
    } else {
        material->matstring[nStr] = '\0';
    }

    fclose(file);

    /* Optional companion files: <inputfile>_ex, _tmelt, _strength. */
    char exinputfile[sizeof(material->matName) + 16];

    snprintf(exinputfile, sizeof(exinputfile), "%s_ex", inputfile);
    if (ANEOSReadExtendedTable(material, exinputfile)) {
        //fprintf(stderr, "No extended EOS tables found.\n");
    }

    snprintf(exinputfile, sizeof(exinputfile), "%s_tmelt", inputfile);
    if (ANEOSReadMeltCurve(material, exinputfile)) {
        //fprintf(stderr, "No melting curve table found.\n");
    }

    snprintf(exinputfile, sizeof(exinputfile), "%s_strength", inputfile);
    if (ANEOSReadStrengthParameters(material, exinputfile)) {
        //fprintf(stderr, "No strength parameters found.\n");
    }

    /*
     * Models 1 and 2 soften the strength with T/T_melt, so a strength file that asks for them
     * without a melt curve is a misconfigured material: stop here, like a bad strength file
     * (ANEOSReadStrengthParameters), instead of running it without strength. Up to 1.2.0 the
     * model was set to -1 silently, for model 0 too, which needs no melt curve and keeps its
     * meaning (no limiter, not "no yield stress"). Without a strength file the model is -1.
     */
    if (material->yieldStrengthModel >= 1 && material->T_melt == NULL) {
        fprintf(stderr, "ANEOSinitMaterialFromFile: %s_strength sets yieldStrengthModel = %d, which needs "
                "the melt curve %s_tmelt, but that file was not found.\n",
                inputfile, material->yieldStrengthModel, inputfile);
        assert(0);
        exit(1);
    }

    return material;
}

/*
 * Releases the memory of ANEOSMATERIAL
 */
void ANEOSfinalizeMaterial(ANEOSMATERIAL *material)
{
    free(material->uData); free(material->uArray);
    free(material->pData); free(material->pArray);
    free(material->cData); free(material->cArray);
    free(material->sData); free(material->sArray);

    if (material->PhaseArray != NULL) {
        free(material->phaseData);
        free(material->PhaseArray);
    }
    if (material->T_melt != NULL) {
        free(material->T_melt);
    }

    free(material->rhoAxis);
    free(material->TAxis);
    free(material);
}

/*
 * Read an extended EOS table.
 */
int ANEOSReadExtendedTable(ANEOSMATERIAL *material, char *inputfile)
{
    FILE *file;
    int nRho;
    int nT;

    if (material == NULL) {
        fprintf(stderr, "ANEOSReadExtendedTable: ANEOSmaterial not initialized.\n");
        assert(0);
        exit(1);
    }

    if ((file= fopen(inputfile, "rb")) == NULL)
    {
        return 1;   /* optional file, absence is not an error */
    }

    /* If the file exists we assume it should be readable. */
    ANEOSReadOrDie(&nRho, sizeof(nRho), 1, file, "nRho", inputfile);
    ANEOSReadOrDie(&nT,   sizeof(nT),   1, file, "nT",   inputfile);

    if ((material->nRho != nRho) || (material->nT != nT)) {
        fprintf(stderr, "ANEOSReadExtendedTable: Inconsistent number of grid points in %s "
                "(table: nRho = %d, nT = %d; main: nRho = %d, nT = %d).\n",
                inputfile, nRho, nT, material->nRho, material->nT);
        assert(0);
        exit(1);
    }

    /*
     * Read the axes and verify that they match the main table instead of
     * skipping them. An extended table built from a different grid would
     * otherwise be silently misaligned.
     */
    double *axisBuf = (double *)malloc(sizeof(double) * (nRho > nT ? nRho : nT));
    assert(axisBuf != NULL);

    ANEOSReadOrDie(axisBuf, sizeof(double), nRho, file, "rhoAxis", inputfile);
    for (int i = 0; i < nRho; i++) {
        if (axisBuf[i] != material->rhoAxis[i]) {
            fprintf(stderr, "ANEOSReadExtendedTable: rhoAxis[%d] = %.15e in %s does not match main table (%.15e).\n",
                    i, axisBuf[i], inputfile, material->rhoAxis[i]);
            assert(0);
            exit(1);
        }
    }
    ANEOSReadOrDie(axisBuf, sizeof(double), nT, file, "TAxis", inputfile);
    for (int j = 0; j < nT; j++) {
        if (axisBuf[j] != material->TAxis[j]) {
            fprintf(stderr, "ANEOSReadExtendedTable: TAxis[%d] = %.15e in %s does not match main table (%.15e).\n",
                    j, axisBuf[j], inputfile, material->TAxis[j]);
            assert(0);
            exit(1);
        }
    }
    free(axisBuf);

    int *phaseData = (int *)malloc((size_t)nRho * nT * sizeof(int));
    assert(phaseData != NULL);
    int **PhaseArray = (int **)malloc((size_t)nRho * sizeof(int *));
    assert(PhaseArray != NULL);
    for (int i = 0; i < nRho; i++) PhaseArray[i] = phaseData + (size_t)i * nT;

    int *tmp = (int *)malloc((size_t)nRho * nT * sizeof(int));
    assert(tmp != NULL);
    ANEOSReadOrDie(tmp, sizeof(int), (size_t)nRho * nT, file, "PhaseArray", inputfile);
    for (int j = 0; j < nT; j++)
        for (int i = 0; i < nRho; i++)
            PhaseArray[i][j] = tmp[(size_t)j * nRho + i];
    free(tmp);

    /* Only publish the array once it is complete. */
    material->phaseData = phaseData;
    material->PhaseArray = PhaseArray;
    fclose(file);

    return 0;
}

/*
 * Read the melt curve.
 */
int ANEOSReadMeltCurve(ANEOSMATERIAL *material, char *inputfile)
{
    FILE *file;
    int nRho;

    if (material == NULL) {
        fprintf(stderr, "ANEOSReadMeltCurve: ANEOSmaterial not initialized.\n");
        assert(0);
        exit(1);
    }

    if ((file= fopen(inputfile, "rb")) == NULL)
    {
        return 1;   /* optional file */
    }

    ANEOSReadOrDie(&nRho, sizeof(nRho), 1, file, "nRho", inputfile);

    if (material->nRho != nRho) {
        fprintf(stderr, "ANEOSReadMeltCurve: Inconsistent number of grid points in %s "
                "(table: %d, main: %d).\n", inputfile, nRho, material->nRho);
        assert(0);
        exit(1);
    }

    double *T_melt = (double *)malloc(sizeof(double)*nRho);
    assert(T_melt != NULL);

    ANEOSReadOrDie(T_melt, sizeof(T_melt[0]), nRho, file, "T_melt", inputfile);

    /* The melt curve is 0 where there is no condensed phase; negative or non-finite is a defect. */
    for (int i = 0; i < nRho; i++) {
        if (!isfinite(T_melt[i]) || !(T_melt[i] >= 0.0)) {
            fprintf(stderr, "ANEOSReadMeltCurve: T_melt[%d] = %g in %s is negative or not finite.\n", i, T_melt[i], inputfile);
            assert(0);
            exit(1);
        }
    }

    fclose(file);
    material->T_melt = T_melt;

    return 0;
}

#define ANEOS_STRENGTH_LINE_MAX 256

/*
 * Read the strength parameters from a text file of the form
 *
 *   # comment
 *   yieldStrengthModel = 2
 *   Y0  = 1e8
 *   YM  = 3.5e10
 *   mui = 1.5
 *   mud = 0.6
 *   xi  = 1.2
 *   Gamma = -1      # -1: compute shear modulus from the sound speed
 *
 * Whitespace around '=' is ignored; '#' starts a comment anywhere on a line.
 * Returns 1 if the file does not exist (no strength available), 0 on success.
 * Malformed files are fatal, and so are parameters the yield model cannot use,
 * so that a bad file stops the run at initialisation and not at the first
 * strength evaluation:
 *   - every value must be finite (strtod accepts "nan" and "inf");
 *   - yieldStrengthModel must be 0, 1 or 2 (the models EOSYieldStrength knows);
 *   - Y0 >= 0, YM >= 0, mui >= 0, xi >= 0, Gamma >= 0 (or the sentinel -1);
 *   - for models 1 and 2: YM > Y0 (otherwise the intact strength
 *     Y0 + mui P/(1 + mui P/(YM - Y0)) has a pole at P = (Y0 - YM)/mui and changes
 *     sign), mui > 0, xi > 0 (xi = 0 makes the thermal softening tanh(0) = 0, i.e.
 *     no strength at all) and a fixed Gamma > 0;
 *   - for model 2: mud >= 0 must be given (the damaged strength mud*P).
 */
int ANEOSReadStrengthParameters(ANEOSMATERIAL *material, char *inputfile)
{
    FILE *file;
    char buffer[ANEOS_STRENGTH_LINE_MAX];
    char key[ANEOS_STRENGTH_LINE_MAX];
    char value[ANEOS_STRENGTH_LINE_MAX];
    int iLine = 0;

    if (material == NULL) {
        fprintf(stderr, "ANEOSReadStrengthParameters: ANEOSmaterial not initialized.\n");
        assert(0);
        exit(1);
    }

    if ((file = fopen(inputfile, "r")) == NULL)
    {
        material->yieldStrengthModel = -1; /* No yield strength available */
        return 1;
    }

    material->yieldStrengthModel = -1;
    material->Y0 = -1.0;
    material->YM = -1.0;
    material->mui = -1.0;
    material->mud = -1.0;
    material->xi = -1.0;
    material->Gamma = -2.0;
    material->fixedGamma = 1;

    while (fgets(buffer, sizeof(buffer), file)) {
        iLine++;

        /* Reject lines longer than the buffer rather than silently splitting them. */
        if (strchr(buffer, '\n') == NULL && !feof(file)) {
            fprintf(stderr, "ANEOSReadStrengthParameters: %s line %d exceeds %d characters.\n",
                    inputfile, iLine, ANEOS_STRENGTH_LINE_MAX - 1);
            fclose(file);
            assert(0);
            exit(1);
        }

        /* Strip comments. */
        char *hash = strchr(buffer, '#');
        if (hash) *hash = '\0';

        /* Skip blank lines (only whitespace left). */
        if (sscanf(buffer, " %255s", key) != 1) continue;

        /* key = value, whitespace around '=' allowed. */
        if (sscanf(buffer, " %255[^= \t\r\n] = %255s", key, value) != 2) {
            fprintf(stderr, "ANEOSReadStrengthParameters: cannot parse %s line %d: '%s'\n",
                    inputfile, iLine, buffer);
            fclose(file);
            assert(0);
            exit(1);
        }

        char *end = NULL;
        if (strcmp(key, "yieldStrengthModel") == 0) {
            material->yieldStrengthModel = (int) strtol(value, &end, 10);
        } else {
            double v = strtod(value, &end);
            if      (strcmp(key, "Y0")  == 0) material->Y0  = v;
            else if (strcmp(key, "YM")  == 0) material->YM  = v;
            else if (strcmp(key, "mui") == 0) material->mui = v;
            else if (strcmp(key, "mud") == 0) material->mud = v;
            else if (strcmp(key, "xi")  == 0) material->xi  = v;
            else if (strcmp(key, "Gamma") == 0) {
                if (v == -1.0) {
                    /* Sentinel: derive the shear modulus from the sound speed. */
                    material->fixedGamma = 0;
                    material->Gamma = 0.0;
                } else {
                    material->fixedGamma = 1;
                    material->Gamma = v;
                }
            } else {
                fprintf(stderr, "ANEOSReadStrengthParameters: unknown parameter '%s' in %s line %d\n",
                        key, inputfile, iLine);
                fclose(file);
                assert(0);
                exit(1);
            }
        }

        if (end == value || *end != '\0') {
            fprintf(stderr, "ANEOSReadStrengthParameters: invalid value '%s' for '%s' in %s line %d\n",
                    value, key, inputfile, iLine);
            fclose(file);
            assert(0);
            exit(1);
        }
    }

    fclose(file);

    /* Missing parameters keep their negative initial values. NaN fails every comparison. */
    if ((material->yieldStrengthModel < 0) || (material->yieldStrengthModel > 2) ||
        !(material->Y0 >= 0.0) || !(material->YM >= 0.0) || !(material->mui >= 0.0) ||
        !(material->xi >= 0.0) || !(material->Gamma >= 0.0) ||
        !isfinite(material->Y0) || !isfinite(material->YM) || !isfinite(material->mui) ||
        !isfinite(material->xi) || !isfinite(material->Gamma) || !isfinite(material->mud)) {
        fprintf(stderr, "ANEOSReadStrengthParameters: missing or invalid parameter in %s "
                "(yieldStrengthModel = %d (0, 1 or 2) Y0 = %g YM = %g mui = %g xi = %g Gamma = %g mud = %g)\n",
                inputfile, material->yieldStrengthModel, material->Y0, material->YM,
                material->mui, material->xi, material->Gamma, material->mud);
        assert(0);
        exit(1);
    }

    if (material->yieldStrengthModel >= 1) {
        if (!(material->YM > material->Y0) || !(material->mui > 0.0) || !(material->xi > 0.0) ||
            (material->fixedGamma && !(material->Gamma > 0.0))) {
            fprintf(stderr, "ANEOSReadStrengthParameters: yieldStrengthModel = %d in %s needs YM > Y0, "
                    "mui > 0, xi > 0 and a positive fixed Gamma (Y0 = %g YM = %g mui = %g xi = %g Gamma = %g)\n",
                    material->yieldStrengthModel, inputfile, material->Y0, material->YM,
                    material->mui, material->xi, material->Gamma);
            assert(0);
            exit(1);
        }
        if (material->yieldStrengthModel == 2 && !(material->mud >= 0.0)) {
            fprintf(stderr, "ANEOSReadStrengthParameters: yieldStrengthModel = 2 in %s needs mud >= 0 "
                    "(damaged material), got mud = %g (missing?)\n", inputfile, material->mud);
            assert(0);
            exit(1);
        }
    }

    return 0;
}

/*
 * Converts the reference density in cgs to Code Units
 */
double ANEOSgetRho0(ANEOSMATERIAL *material)
{
    return material->rho0 * material->CGStoCodeUnitsforRho;
}

/*
 * The lowest tabulated temperature.
 */
double ANEOSTmin(ANEOSMATERIAL *material)
{
    return material->TAxis[0];
}

/*
 * The cold curve of the table, u(rho, T_min).
 */
double ANEOSUCold(ANEOSMATERIAL *material, double rho)
{
    return ANEOSUofRhoT(material, rho, material->TAxis[0]);
}

/*
 * Table status of (rho, u), see ANEOSmaterial.h. One search on the rho axis; the u bounds are the
 * interpolated first and last rows of that column, i.e. u(rho, T_min) and u(rho, T_max), compared
 * with the relative tolerance ANEOS_EDGE_TOL (the same rule as the T(rho, u) inversion's table edges).
 */
int ANEOSIsInTable(ANEOSMATERIAL *material, double rho, double u)
{
    const double rho_cgs = rho*material->CodeUnitstoCGSforRho;
    const double u_cgs = u*material->CodeUnitstoCGSforU;
    double w, uMin, uMax;
    int i;

    switch (interpLocateAxis(material->rhoAxis, material->nRho, rho_cgs, ANEOS_EDGE_TOL, &i, &w)) {
        case INTERP_SUCCESS: break;
        case INTERP_ABOVE:   return ANEOS_OUTSIDE_RHOMAX;
        default:             return ANEOS_OUTSIDE_RHOMIN;   /* below, or NaN */
    }
    uMin = interpLineValue(&material->grid, material->uData, INTERP_ALONG_Y, i, w, 0);
    uMax = interpLineValue(&material->grid, material->uData, INTERP_ALONG_Y, i, w, material->nT - 1);
    if (!(u_cgs >= uMin - ANEOS_EDGE_TOL*fabs(uMin))) return ANEOS_OUTSIDE_TMIN;   /* also NaN */
    if (u_cgs > uMax + ANEOS_EDGE_TOL*fabs(uMax)) return ANEOS_OUTSIDE_TMAX;
    return ANEOS_SUCCESS;
}

/*
 * Calculates the temperature T(p,u)
 */
double ANEOSTofPU(ANEOSMATERIAL *material, double p, double u)
{
    double rho = ANEOSRhoofPU(material, p, u);
    if (ANEOSFailed(rho)) return ANEOS_FAIL_VALUE;
    return ANEOSTofRhoU(material, rho, u);
}

/*
 * Locate the state (rho, u) in the table: the density cell and the inversion of the u column
 * there (the same two steps as interpInvertY), giving the cell (ix along rho, iy along T) with its
 * local coordinates and T. Every (rho, u) function evaluates its tables in this one cell, as
 * ANEOSPCTofRhoU always did (scvheos has the same helper, scvheosLocateRhoU):
 *
 *   - one search instead of two (T(rho, u), then a second locate of (rho, T) for P, c, s or the
 *     partials), and no roundoff from recomputing the T coordinate of the cell from T;
 *   - one bad state is counted once in nLookupFail. Up to 1.2.0 a density outside the table was
 *     counted twice (the T fallback, then the failed lookup at (rho, T_min)) by ANEOSPofRhoU,
 *     ANEOSCofRhoU, ANEOSSofRhoU, ANEOSPhaseofRhoU and the derivatives, once by ANEOSPCTofRhoU.
 *
 * Returns 0 on success, -1 (counted, message under EOSLIB_VERBOSE) if rho is outside the table or
 * not a number. If u cannot be inverted (outside the column [u(T_min), u(T_max)], or NaN) the state
 * falls back to T_min by contract: cell iy = 0, v = 0, T = TAxis[0], counted, and 0 is returned,
 * so the values are those of (rho, T_min).
 */
static int ANEOSLocateRhoU(ANEOSMATERIAL *material, const char *chFunc, double rho, double u,
                           INTERPCELL *cell, double *pT)
{
    const double rho_cgs = rho*material->CodeUnitstoCGSforRho;
    const double u_cgs = u*material->CodeUnitstoCGSforU;

    if (interpLocateAxis(material->rhoAxis, material->nRho, rho_cgs, ANEOS_EDGE_TOL, &cell->ix, &cell->u) != INTERP_SUCCESS) {
        material->nLookupFail++;
        ANEOS_WARN("%s: rho = %.15e outside the table (u = %.15e)\n", chFunc, rho, u);
        *pT = material->TAxis[0];
        return -1;
    }

    if (interpInvertLine(&material->grid, material->uData, INTERP_ALONG_Y, cell->ix, cell->u, u_cgs,
                         &ANEOSInvertOpt, &cell->iy, &cell->v) != INTERP_SUCCESS) {
        /* Fallback, by contract: the lowest tabulated temperature. */
        material->nLookupFail++;
        ANEOS_WARN("%s: T(rho,u) failed for rho = %.15e, u = %.15e, falling back to T = %.15e\n",
                   chFunc, rho, u, material->TAxis[0]);
        cell->iy = 0;
        cell->v = 0.0;
    }

    *pT = interpLineCoord(&material->grid, INTERP_ALONG_Y, cell->iy, cell->v);
    return 0;
}

/*
 * Calculates temperature T(rho,u). Never fails: T_min if rho is outside the table or u cannot be
 * inverted (counted).
 */
double ANEOSTofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
    INTERPCELL cell;
    double T;
    (void) ANEOSLocateRhoU(material, __func__, rho, u, &cell, &T);   /* T = T_min on failure */
    return T;
}

/*
 * Calculates temperature T(rho,p)
 */
double ANEOSTofRhoP(ANEOSMATERIAL *material, double rho, double p)
{
    INTERPCELL cell;
    double T;
    if (interpInvertY(&material->grid, material->pData, rho*material->CodeUnitstoCGSforRho, p*material->CodeUnitstoCGSforP,
                      ANEOS_EDGE_TOL, &ANEOSInvertOpt, &cell, &T) != INTERP_SUCCESS) {
        material->nLookupFail++;
        ANEOS_WARN("ANEOSTofRhoP failed for rho = %.15e, p = %.15e\n", rho, p);
        return ANEOS_FAIL_VALUE;
    }
    return T;
}

/*
 * Calculates internal energy u(p,T).
 *
 * The chained lookups (u(p, T), u(rho, p), p(u, T)) stop at a failed intermediate: they return
 * ANEOS_FAIL_VALUE, counted once by the call that failed, instead of passing the sentinel on to
 * the next lookup (which failed again and counted it a second time).
 */
double ANEOSUofPT(ANEOSMATERIAL *material, double p, double T)
{
    double rho = ANEOSRhoofPT(material,p,T);
    if (ANEOSFailed(rho)) return ANEOS_FAIL_VALUE;
    return ANEOSUofRhoT(material,rho,T);
}

/*
 * Calculates internal energy u(rho,T)
 */
double ANEOSUofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
    double u;
    if (interpValue(&material->grid, material->uData, rho*material->CodeUnitstoCGSforRho, T, ANEOS_EDGE_TOL, &u) != INTERP_SUCCESS) {
        material->nLookupFail++;
        ANEOS_WARN("ANEOSUofRhoT failed for rho = %.15e, T = %.15e\n", rho, T);
        return ANEOS_FAIL_VALUE;
    }
    return u * material->CGStoCodeUnitsforU;
}

/*
 * Calculates internal energy u(rho,p)
 */
double ANEOSUofRhoP(ANEOSMATERIAL *material, double rho, double p)
{
    double T = ANEOSTofRhoP(material,rho,p);
    if (ANEOSFailed(T)) return ANEOS_FAIL_VALUE;   /* see ANEOSUofPT */
    return ANEOSUofRhoT(material,rho,T);
}

/*
 * Calculates pressure p(u,T)
 */
double ANEOSPofUT(ANEOSMATERIAL *material, double u, double T)
{
    double rho = ANEOSRhoofUT(material,u,T);
    if (ANEOSFailed(rho)) return ANEOS_FAIL_VALUE;   /* see ANEOSUofPT */
    return ANEOSPofRhoT(material,rho,T);
}

/*
 * Calculates pressure p(rho,T)
 */
double ANEOSPofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
    double p;
    if (interpValue(&material->grid, material->pData, rho*material->CodeUnitstoCGSforRho, T, ANEOS_EDGE_TOL, &p) != INTERP_SUCCESS) {
        material->nLookupFail++;
        ANEOS_WARN("ANEOSPofRhoT failed for rho = %.15e, T = %.15e\n", rho, T);
        return ANEOS_FAIL_VALUE;
    }
    return p * material->CGStoCodeUnitsforP;
}

/*
 * Calculates pressure p(rho,u)
 */
double ANEOSPofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
    INTERPCELL cell;
    double T;
    if (ANEOSLocateRhoU(material, __func__, rho, u, &cell, &T)) return ANEOS_FAIL_VALUE;
    return interpCellValue(&material->grid, material->pData, &cell) * material->CGStoCodeUnitsforP;
}

/*
 * Calculates density rho(u,T)
 */
double ANEOSRhoofUT(ANEOSMATERIAL *material, double u, double T)
{
    INTERPCELL cell;
    double rho;
    if (interpInvertX(&material->grid, material->uData, T, u*material->CodeUnitstoCGSforU,
                      ANEOS_EDGE_TOL, &ANEOSInvertOpt, &cell, &rho) != INTERP_SUCCESS) {
        material->nLookupFail++;
        ANEOS_WARN("ANEOSRhoofUT failed for u = %.15e, T = %.15e\n", u, T);
        return ANEOS_FAIL_VALUE;
    }
    return rho * material->CGStoCodeUnitsforRho;
}

/*
 * Calculates density rho(p,T)
 */
double ANEOSRhoofPT(ANEOSMATERIAL *material, double p, double T)
{
    INTERPCELL cell;
    double rho;
    if (interpInvertX(&material->grid, material->pData, T, p*material->CodeUnitstoCGSforP,
                      ANEOS_EDGE_TOL, &ANEOSInvertOpt, &cell, &rho) != INTERP_SUCCESS) {
        material->nLookupFail++;
        ANEOS_WARN("ANEOSRhoofPT failed for p = %.15e, T = %.15e\n", p, T);
        return ANEOS_FAIL_VALUE;
    }
    return rho * material->CGStoCodeUnitsforRho;
}

/*
 * g(rho) = P(rho, u) - P in table (cgs) units.
 *
 * *outside reports whether the state (rho, u) exists in the table: 0 if it
 * does, +-1 if u lies below the cold curve u(rho, TAxis[0]) (the state would
 * be colder than the table), +-2 if u lies above u(rho, TAxis[nT-1]) (hotter
 * than the table). The sign is the sign g must be given for bracketing: + when
 * the attainable region lies at lower density (the cold curve increases, or
 * the top row decreases, with rho), - when it lies at higher density. In both
 * cases 0.0 is returned as value; the caller uses only *outside.
 *
 * Returns -1e50 (with *outside = 0) when u is unattainable for any other
 * reason: rho outside the table, or u in a gap of a non-monotonic u(T) line.
 */
static double ANEOSRhoofPUResidual(ANEOSMATERIAL *material, double rho_cgs, double p_cgs, double u_cgs, int *outside)
{
    INTERPCELL cell;
    int i, iret;
    double x, T;

    *outside = 0;
    iret = interpInvertY(&material->grid, material->uData, rho_cgs, u_cgs, ANEOS_EDGE_TOL, &ANEOSInvertOpt, &cell, &T);
    if (iret == INTERP_SUCCESS)
        return interpCellValue(&material->grid, material->pData, &cell) - p_cgs;

    /* rho outside the table (or NaN): no state at all. Otherwise u is not attained in this column. */
    if (interpLocateAxis(material->rhoAxis, material->nRho, rho_cgs, ANEOS_EDGE_TOL, &i, &x) != INTERP_SUCCESS) return ANEOS_FAIL_VALUE;
    int jn = material->nT - 1;
    double u0 = interpLineValue(&material->grid, material->uData, INTERP_ALONG_Y, i, x, 0);
    double un = interpLineValue(&material->grid, material->uData, INTERP_ALONG_Y, i, x, jn);
    if (u_cgs < u0)      *outside = (material->uArray[i+1][0]  > material->uArray[i][0])  ? +1 : -1;
    else if (u_cgs > un) *outside = (material->uArray[i+1][jn] > material->uArray[i][jn]) ? -2 : +2;
    else return ANEOS_FAIL_VALUE;
    return 0.0;
}

/*
 * Refine a sign change of g between the densities a and b by bisection.
 * Returns 1 and the root in *prho (cgs) on success, 0 if the bracket turned
 * out to be an artefact of the table edge (the caller keeps scanning), and
 * -1 if the residual became undefined (gap of a non-monotonic u(T) line).
 */
static int ANEOSRhoofPURefine(ANEOSMATERIAL *material, double p_cgs, double u_cgs,
        double a, double ga, int outA, double b, double gb, int outB, double *prho)
{
    const double tol = 1e-10;
    const int maxIter = 100;

    /* Reference for the final acceptance test: an attainable end of the
     * initial bracket (or the first attainable midpoint) and the bracket width. */
    double w0 = b - a, g0 = 0.0;
    int haveRef = 0;
    if (outA == 0)      { g0 = ga; haveRef = 1; }
    else if (outB == 0) { g0 = gb; haveRef = 1; }
    int negA = outA ? (outA < 0) : (ga < 0.0);

    for (int n = 0; n < maxIter && 2.0 * (b - a) / (a + b) > tol; n++) {
        int out;
        double c = 0.5 * (a + b);
        double gc = ANEOSRhoofPUResidual(material, c, p_cgs, u_cgs, &out);
        if (ANEOSFailed(gc)) return -1;
        if (out == 0 && gc == 0.0) { *prho = c; return 1; }
        if (out == 0 && !haveRef) { g0 = gc; haveRef = 1; }
        int neg = out ? (out < 0) : (gc < 0.0);
        if (neg == negA) { a = c; ga = gc; outA = out; }
        else             { b = c; gb = gc; outB = out; }
    }

    if (outA == 0 && outB == 0) { *prho = 0.5 * (a + b); return 1; }
    if (outA != 0 && outB != 0) return 0;   /* attainable band narrower than tol, or none */
    /* One end on the edge of the attainable region: the root is genuine only
     * if the residual there has shrunk with the bracket. */
    double gx = outA ? gb : ga;
    if (fabs(gx) > 8.0 * fabs(g0) * (b - a) / w0) return 0;
    *prho = outA ? b : a;
    return 1;
}

/*
 * Calculates density rho(p,u): the density at which P(rho, u) = p.
 *
 * The residual is evaluated at every tabulated density, scanning from low to
 * high density; the first interval where it changes sign is refined by
 * bisection to a relative tolerance. If several roots exist the lowest-density
 * one is returned (library convention).
 *
 * For a cold compressed state the next tabulated density above the root is
 * already too cold for u (u lies below the cold curve there), so the root
 * could not be bracketed between two attainable densities. Densities at
 * which u is outside the tabulated T range therefore take part with the
 * sign of the side on which the attainable region lies (see
 * ANEOSRhoofPUResidual). A sign change between two densities that are both
 * too cold, or both too hot, is an artefact and is ignored; between a too
 * cold and a too hot one the attainable band lies in between and is found
 * by the bisection. A bracket whose bisection ends on the boundary of the
 * attainable region without the residual vanishing there is an artefact of
 * the table edge as well, and the scan continues past it.
 *
 * Close to the cold curve the tabulated P(rho, u) may decrease with rho
 * (large effective Gamma at low T), so rho(P, u) can be double-valued
 * there; the bisection then finds the lower root only if it does not step
 * over it, and otherwise reports failure.
 *
 * Returns rho in code units, or -1e50 if no root exists.
 */
double ANEOSRhoofPU(ANEOSMATERIAL *material, double p, double u)
{
    double p_cgs = p * material->CodeUnitstoCGSforP;
    double u_cgs = u * material->CodeUnitstoCGSforU;
    int nRho = material->nRho;

    double gprev = 0.0;
    int outPrev = 0, iprev = -1;
    /*
     * Every node, the last one included: the table is closed at rhoAxis[nRho-1], so the residual
     * is defined there, and without it a root in the last cell is never bracketed (up to 1.2.0
     * the loop stopped at nRho-2 and reported "no root found" for such a state).
     */
    for (int i = 0; i < nRho; i++) {
        int out;
        double g = ANEOSRhoofPUResidual(material, material->rhoAxis[i], p_cgs, u_cgs, &out);
        if (g < -1e40) continue;
        if (out == 0 && g == 0.0) return material->rhoAxis[i] * material->CGStoCodeUnitsforRho;
        int neg = out ? (out < 0) : (g < 0.0);
        if (iprev >= 0) {
            int negPrev = outPrev ? (outPrev < 0) : (gprev < 0.0);
            if (neg != negPrev && !(out != 0 && abs(out) == abs(outPrev))) {
                double rho;
                int r = ANEOSRhoofPURefine(material, p_cgs, u_cgs,
                            material->rhoAxis[iprev], gprev, outPrev, material->rhoAxis[i], g, out, &rho);
                if (r > 0) return rho * material->CGStoCodeUnitsforRho;
                if (r < 0) {
                    material->nLookupFail++;
                    ANEOS_WARN("ANEOSRhoofPU: residual undefined between rho = %.15e and %.15e for p = %.15e, u = %.15e\n",
                               material->rhoAxis[iprev] * material->CGStoCodeUnitsforRho,
                               material->rhoAxis[i] * material->CGStoCodeUnitsforRho, p, u);
                    return ANEOS_FAIL_VALUE;
                }
            }
        }
        gprev = g; outPrev = out; iprev = i;
    }
    material->nLookupFail++;
    ANEOS_WARN("ANEOSRhoofPU: no root found for p = %.15e, u = %.15e\n", p, u);
    return ANEOS_FAIL_VALUE;
}

/*
 * Calculates the sound speed c(rho,u)
 */
double ANEOSCofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
    INTERPCELL cell;
    double T;
    if (ANEOSLocateRhoU(material, __func__, rho, u, &cell, &T)) return ANEOS_FAIL_VALUE;
    return interpCellValue(&material->grid, material->cData, &cell) * material->CGStoCodeUnitsforC;
}

/*
 * Calculates the sound speed c(rho,T)
 */
double ANEOSCofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
    double c;
    if (interpValue(&material->grid, material->cData, rho*material->CodeUnitstoCGSforRho, T, ANEOS_EDGE_TOL, &c) != INTERP_SUCCESS) {
        material->nLookupFail++;
        ANEOS_WARN("ANEOSCofRhoT failed for rho = %.15e, T = %.15e\n", rho, T);
        return ANEOS_FAIL_VALUE;
    }
    return c * material->CGStoCodeUnitsforC;
}

/*
 * P, c and T of (rho, u) from one lookup (ANEOSLocateRhoU). rho outside the table: *T = T_min,
 * *c and the return value ANEOS_FAIL_VALUE, counted once; the separate calls give the same values.
 */
double ANEOSPCTofRhoU(ANEOSMATERIAL *material, double rho, double u, double *c, double *T)
{
    INTERPCELL cell;

    if (ANEOSLocateRhoU(material, __func__, rho, u, &cell, T)) {
        *c = ANEOS_FAIL_VALUE;
        return ANEOS_FAIL_VALUE;
    }
    *c = interpCellValue(&material->grid, material->cData, &cell) * material->CGStoCodeUnitsforC;
    return interpCellValue(&material->grid, material->pData, &cell) * material->CGStoCodeUnitsforP;
}

/*
 * Calculates the phase(rho,u) from the extended EOS table. See ANEOSPhaseofRhoT.
 */
int ANEOSPhaseofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
    INTERPCELL cell;
    double T;

    if (material->PhaseArray == NULL) {
        material->nLookupFail++;
        ANEOS_WARN("ANEOSPhaseofRhoU: no extended EOS table loaded for material %d.\n", material->iMat);
        return -1;
    }
    if (ANEOSLocateRhoU(material, __func__, rho, u, &cell, &T)) return -1;

    /* The node nearest to the state in its cell; exactly halfway the lower node wins (interpNearest). */
    return material->phaseData[(size_t) (cell.ix + (cell.u > 0.5)) * material->nT + (cell.iy + (cell.v > 0.5))];
}

/*
 * Calculates the phase(rho,T) from the extended EOS table (nearest neighbour).
 *
 * Returns -1 if no extended table was loaded for this material or if (rho, T)
 * lies outside the table.
 */
int ANEOSPhaseofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
    if (material->PhaseArray == NULL) {
        material->nLookupFail++;
        ANEOS_WARN("ANEOSPhaseofRhoT: no extended EOS table loaded for material %d.\n", material->iMat);
        return -1;
    }

    int phase;
    if (interpNearest(&material->grid, material->phaseData, rho*material->CodeUnitstoCGSforRho, T, ANEOS_EDGE_TOL, &phase) != INTERP_SUCCESS) {
        material->nLookupFail++;
        ANEOS_WARN("ANEOSPhaseofRhoT failed for rho = %.15e, T = %.15e\n", rho, T);
        return -1;
    }
    return phase;
}

/*
 * Calculates the isentropic evolution u2(rho1,u1,rho2).
 *
 * If the entropy S(rho1,u1) is not attained at rho2, u2 is clamped to the
 * table edge on the side the isentrope leaves it: to the cold curve
 * u(rho2, TAxis[0]) when S lies below the s(T) line at rho2 (the same
 * convention ANEOSTofRhoU uses for sub-table states), to the top row
 * u(rho2, TAxis[nT-1]) when it lies above (a hot state compressed beyond the
 * table; scvheosIsentropicU does the same). If rho1 is outside the table
 * the cold curve at rho2 is used. A call counts at most once in nLookupFail,
 * however many of these steps fall back. Returns -1e50 only if rho2 itself
 * is outside the table.
 */
double ANEOSisentropicU(ANEOSMATERIAL *material, double rho1, double u1, double rho2)
{
    const double rho2_cgs = rho2*material->CodeUnitstoCGSforRho;
    const long nLookupFail0 = material->nLookupFail;
    double oldT, oldS, newT, u2;
    int bClamp = 0;
    INTERPCELL cell;

    if (ANEOSLocateRhoU(material, __func__, rho1, u1, &cell, &oldT)) {
        /* rho1 outside the table (counted): no entropy, so the cold curve at rho2. */
        newT = material->TAxis[0];
        oldS = 0.0;
    } else {
        /* The entropy in the cell of (rho1, u1), then the T of that entropy at rho2. */
        oldS = interpCellValue(&material->grid, material->sData, &cell);
        bClamp = (interpInvertY(&material->grid, material->sData, rho2_cgs, oldS,
                                ANEOS_EDGE_TOL, &ANEOSInvertOpt, &cell, &newT) != INTERP_SUCCESS);
    }
    if (bClamp) {
        /*
         * S is not attained on the s(T) line at rho2. The tolerant inversion scans every segment,
         * so S lies below the line's minimum or above its maximum. The top row is not above the
         * maximum, so S > s(rho2, T_max) tells the two apart even on a defective (non-monotonic)
         * line. rho2 outside the table: no line, u(rho2, T) below returns the sentinel anyway.
         */
        int i;
        double w;
        material->nLookupFail++;
        newT = material->TAxis[0];
        if (interpLocateAxis(material->rhoAxis, material->nRho, rho2_cgs, ANEOS_EDGE_TOL, &i, &w) == INTERP_SUCCESS &&
            oldS > interpLineValue(&material->grid, material->sData, INTERP_ALONG_Y, i, w, material->nT - 1)) {
            newT = material->TAxis[material->nT - 1];
        }
        ANEOS_WARN("ANEOSisentropicU newT failed for rho1 = %.15e, u1 = %.15e, rho2 = %.15e, using %s at rho2\n",
                   rho1, u1, rho2, (newT == material->TAxis[0]) ? "the cold curve" : "the top row of the table");
    }

    /* ANEOSUofRhoT returns -1e50 if rho2 is outside the table. */
    u2 = ANEOSUofRhoT(material,rho2,newT);

    /*
     * One call, one count: a fallback at rho1, a clamp at rho2 and a failure at rho2 are steps of
     * the same evaluation (up to 1.2.0 rho1 and rho2 outside the table counted 3).
     */
    if (material->nLookupFail > nLookupFail0) material->nLookupFail = nLookupFail0 + 1;
    return u2;
}

/*
 * Calculates the entropy S(rho,u)
 */
double ANEOSSofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
    INTERPCELL cell;
    double T;
    if (ANEOSLocateRhoU(material, __func__, rho, u, &cell, &T)) return ANEOS_FAIL_VALUE;
    return interpCellValue(&material->grid, material->sData, &cell);   /* cgs, see ANEOSSofRhoT */
}

/*
 * Calculates the entropy S(rho,T).
 * Returned in cgs (erg/g/K) directly from the table; no unit conversion by design.
 */
double ANEOSSofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
    double S;
    if (interpValue(&material->grid, material->sData, rho*material->CodeUnitstoCGSforRho, T, ANEOS_EDGE_TOL, &S) != INTERP_SUCCESS) {
        material->nLookupFail++;
        ANEOS_WARN("ANEOSSofRhoT failed for rho = %.15e, T = %.15e\n", rho, T);
        return ANEOS_FAIL_VALUE;
    }
    return S;
}

/*
 * Calculates the melting temperature T_melt(rho).
 *
 * Outside the tabulated density range the melt curve is held constant at the
 * nearest edge value (conservative for the strength model: T_melt increases
 * with compression, so clamping above the table underestimates it).
 *
 * A missing melt curve is a configuration error (a strength model was
 * requested without its melt curve) and fatal. Returns ANEOS_FAIL_VALUE
 * (counted in nLookupFail) if rho is not a valid density.
 */
double ANEOSTmeltofRho(ANEOSMATERIAL *material, double rho)
{
    if (material->T_melt == NULL) {
        fprintf(stderr, "ANEOSTmeltofRho: Melt curve was not loaded for material %d.\n", material->iMat);
        assert(0);
        exit(1);
    }

    double rho_cgs = rho * material->CodeUnitstoCGSforRho;
    if (!(rho_cgs > 0.0) || !isfinite(rho_cgs)) {
        material->nLookupFail++;
        ANEOS_WARN("ANEOSTmeltofRho: invalid rho = %.15e for material %d.\n", rho, material->iMat);
        return ANEOS_FAIL_VALUE;
    }

    /* Clamp to the tabulated range. */
    double rhoMin = material->rhoAxis[0];
    double rhoMax = material->rhoAxis[material->nRho-1];
    if (rho_cgs <= rhoMin) return material->T_melt[0];
    if (rho_cgs >= rhoMax) return material->T_melt[material->nRho-1];

    /*
     * Linear interpolation along the rho axis. rho_cgs lies strictly inside the axis after the
     * clamp, so the locate cannot fail; its status is ignored on purpose (it returned -1.0 here,
     * outside the failure contract, on a path that cannot be reached).
     */
    int i = 0;
    double x = 0.0;
    (void) interpLocateAxis(material->rhoAxis, material->nRho, rho_cgs, 0.0, &i, &x);
    return material->T_melt[i] + x * (material->T_melt[i+1] - material->T_melt[i]);
}

/*
 * Partial derivatives of P(rho,T) and u(rho,T) at fixed T and fixed rho,
 * in cgs. Shared by the public derivative functions below.
 * Returns 0 on success, -1 if (rho, T) is outside the table.
 */
static int ANEOSPartialsCGS(ANEOSMATERIAL *material, double rho, double T,
                            double *dPdrho_T, double *dPdT_rho,
                            double *dudrho_T, double *dudT_rho)
{
    double rho_cgs = rho * material->CodeUnitstoCGSforRho;
    int iret = 0;

    iret |= interpGradient(&material->grid, material->pData, rho_cgs, T, ANEOS_EDGE_TOL, dPdrho_T, dPdT_rho);
    if (dudrho_T && dudT_rho) {
        iret |= interpGradient(&material->grid, material->uData, rho_cgs, T, ANEOS_EDGE_TOL, dudrho_T, dudT_rho);
    }
    if (iret) {
        material->nLookupFail++;
        ANEOS_WARN("ANEOSPartialsCGS failed for rho = %.15e, T = %.15e\n", rho, T);
        return -1;
    }
    return 0;
}

/*
 * dP/drho at constant u (code units).
 *
 *   dP/drho|_u = dP/drho|_T - dP/dT|_rho * (du/drho|_T) / (du/dT|_rho)
 */
double ANEOSdPdRhoofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
    double dPdrho_T, dPdT_rho, dudrho_T, dudT_rho;
    const long nLookupFail0 = material->nLookupFail;
    double T;
    INTERPCELL cell;

    /* Partials of the interpolant in the cell of the state (ANEOSLocateRhoU). */
    if (ANEOSLocateRhoU(material, __func__, rho, u, &cell, &T)) return ANEOS_FAIL_VALUE;
    interpCellGradient(&material->grid, material->pData, &cell, &dPdrho_T, &dPdT_rho);
    interpCellGradient(&material->grid, material->uData, &cell, &dudrho_T, &dudT_rho);
    if (!(dudT_rho > 0.0)) {
        /* Heat capacity must be positive; a non-positive du/dT is a table artefact. */
        /* Counted once per call: not again after a T fallback (the cold rows of M-ANEOS tables have
         * du/dT <= 0, so a state above or below the column fails here as well). */
        material->nLookupFail = nLookupFail0 + 1;
        ANEOS_WARN("ANEOSdPdRhoofRhoU: du/dT|rho = %.15e <= 0 at rho = %.15e, T = %.15e\n", dudT_rho, rho, T);
        return ANEOS_FAIL_VALUE;
    }

    double dPdrho_u = dPdrho_T - dPdT_rho * dudrho_T / dudT_rho;

    /* cgs: erg cm^-3 / (g cm^-3) -> code units. */
    return dPdrho_u * material->CGStoCodeUnitsforP * material->CodeUnitstoCGSforRho;
}

/*
 * dP/du at constant rho (code units).
 *
 *   dP/du|_rho = dP/dT|_rho / du/dT|_rho
 */
double ANEOSdPdUofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
    double dPdrho_T, dPdT_rho, dudrho_T, dudT_rho;
    const long nLookupFail0 = material->nLookupFail;
    double T;
    INTERPCELL cell;

    /* Partials of the interpolant in the cell of the state (ANEOSLocateRhoU). */
    if (ANEOSLocateRhoU(material, __func__, rho, u, &cell, &T)) return ANEOS_FAIL_VALUE;
    interpCellGradient(&material->grid, material->pData, &cell, &dPdrho_T, &dPdT_rho);
    interpCellGradient(&material->grid, material->uData, &cell, &dudrho_T, &dudT_rho);
    if (!(dudT_rho > 0.0)) {
        /* Counted once per call: not again after a T fallback (the cold rows of M-ANEOS tables have
         * du/dT <= 0, so a state above or below the column fails here as well). */
        material->nLookupFail = nLookupFail0 + 1;
        ANEOS_WARN("ANEOSdPdUofRhoU: du/dT|rho = %.15e <= 0 at rho = %.15e, T = %.15e\n", dudT_rho, rho, T);
        return ANEOS_FAIL_VALUE;
    }

    /* cgs: erg cm^-3 / (erg g^-1) -> code units. */
    return (dPdT_rho / dudT_rho) * material->CGStoCodeUnitsforP * material->CodeUnitstoCGSforU;
}

/*
 * du/drho along an isentrope (code units): du/drho|_S = P/rho^2.
 * Returns ANEOS_FAIL_VALUE if P(rho, u) fails (counted there, not again). Up to 1.2.0 the
 * sentinel was divided by rho^2: -1e50/rho^2 is not the sentinel, and for rho > 1e5 (code
 * units) not even ANEOSFailed.
 */
double ANEOSdUdRhoofRhoU(ANEOSMATERIAL *material, double rho, double u)
{
    double P = ANEOSPofRhoU(material, rho, u);
    if (ANEOSFailed(P)) return ANEOS_FAIL_VALUE;
    return P/(rho*rho);
}

/*
 * dP/drho at constant T (code units).
 */
double ANEOSdPdRhoofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
    double dPdrho_T, dPdT_rho;

    if (ANEOSPartialsCGS(material, rho, T, &dPdrho_T, &dPdT_rho, NULL, NULL)) {
        return ANEOS_FAIL_VALUE;
    }
    return dPdrho_T * material->CGStoCodeUnitsforP * material->CodeUnitstoCGSforRho;
}

/*
 * dP/dT at constant rho (code units; T is in K in both unit systems).
 */
double ANEOSdPdTofRhoT(ANEOSMATERIAL *material, double rho, double T)
{
    double dPdrho_T, dPdT_rho;

    if (ANEOSPartialsCGS(material, rho, T, &dPdrho_T, &dPdT_rho, NULL, NULL)) {
        return ANEOS_FAIL_VALUE;
    }
    return dPdT_rho * material->CGStoCodeUnitsforP;
}

/*
 * Return yield strength parameters in code units
 */
int ANEOSYieldParameters(ANEOSMATERIAL *material, double *Y0, double *YM, double *mui, double *mud, double *xi) {
    if (Y0) *Y0 = material->Y0 * material->CGStoCodeUnitsforP;
    if (YM) *YM = material->YM * material->CGStoCodeUnitsforP;
    if (mui) *mui = material->mui;
    if (mud) *mud = material->mud;
    if (xi) *xi = material->xi;
    return material->yieldStrengthModel;
}

/*
 * Return shear modulus in code units
 */
double ANEOSGammaofRhoT(ANEOSMATERIAL *material, double rho, double T) {
    if (material->fixedGamma) {
        return material->Gamma * material->CGStoCodeUnitsforP;
    } else {
        // We guess nu = 0.25 as "normal" materials have nu \in [0.2,0.3]
        // and the difference in the resulting shear modulus is less than 25%
        double nu = 0.25;
        double c = ANEOSCofRhoT(material, rho, T);
        if (c < -1e40) return ANEOS_FAIL_VALUE;
        double bulkModulus = rho * c * c;
        return bulkModulus * 3.0 * (1.0 - 2.0 * nu) / (2.0 * (1.0 + nu));
    }
}

void ANEOSPrintMat(ANEOSMATERIAL *material, FILE *fp)
{
    fprintf(fp,"# Material: %i (%s)\n", material->iMat, material->matName);
    fprintf(fp,"# Reference density rho0: %g\n", material->rho0 * material->CGStoCodeUnitsforRho);
    fprintf(fp,"# Table size: nRho = %d, nT = %d\n", material->nRho, material->nT);
    fprintf(fp,"# Table boundaries (cgs units): minRho = %g, maxRho = %g, minT = %g, maxT = %g\n",material->rhoAxis[0],material->rhoAxis[material->nRho-1],material->TAxis[0],material->TAxis[material->nT-1]);
    fprintf(fp,"# %s\n",material->matstring);
    fprintf(fp,"# Failed lookups and fallbacks (nLookupFail): %ld\n", material->nLookupFail);
}