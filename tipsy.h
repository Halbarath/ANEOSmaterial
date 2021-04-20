/*  This file is part of ANEOSmaterial.
 *  Copyright (c) 2000 Joachim Stadel
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

#ifndef TIPSY_INCLUDED
#define TIPSY_INCLUDED

#include <rpc/types.h>
#include <rpc/xdr.h>

#define TIPSY_TYPE_GAS	1
#define TIPSY_TYPE_DARK	2
#define TIPSY_TYPE_STAR	3

#define TIPSY_BLOCK 10000

struct gas_particle {
    float mass;
#ifdef DOUBLE_POS
    double pos[3];
#else
    float pos[3];
#endif
    float vel[3];
    float rho;
    float temp;
    float hsmooth;
    float metals;
    float phi;
    };

struct dark_particle {
    float mass;
#ifdef DOUBLE_POS
    double pos[3];
#else
    float pos[3];
#endif
    float vel[3];
    float eps;
    float phi;
    };

struct star_particle {
    float mass;
#ifdef DOUBLE_POS
    double pos[3];
#else
    float pos[3];
#endif
    float vel[3];
    float metals;
    float tform;
    float eps;
    float phi;
    };

struct base_particle {
    float mass;
#ifdef DOUBLE_POS
    double pos[3];
#else
    float pos[3];
#endif
    float vel[3];
    };

struct dump {
    double time;
    unsigned nbodies;
    unsigned ndim;
    unsigned nsph;
    unsigned ndark;
    unsigned nstar;
    };


typedef struct TipsyContext {
    int bNative;
    unsigned iCurr;
    int iCurrType;
    unsigned nRemain;
    unsigned nBlock;
    unsigned iBlock;
    struct gas_particle gpCurr;
    struct dark_particle dpCurr;
    struct star_particle spCurr;
    unsigned nMaxGas,nGas;
    unsigned nMaxDark,nDark;
    unsigned nMaxStar,nStar;
    struct gas_particle *gp;
    struct dark_particle *dp;
    struct star_particle *sp;
    FILE *fp;
    XDR xdr;
    struct dump head;
    char *pchInFile;
    } * TCTX;


int xdr_header(XDR *xdrs, struct dump *header);
int xdr_gas(XDR *xdrs,struct gas_particle *p);
int xdr_dark(XDR *xdrs,struct dark_particle *p);
int xdr_star(XDR *xdrs,struct star_particle *p);

void TipsyInitialize(TCTX *pctx,int bNative,char *pchInFile);
void TipsyFinish(TCTX ctx);
struct base_particle *pTipsyReadNative(TCTX ctx,int *piType,double *pdSoft);
struct base_particle *pTipsyRead(TCTX ctx,int *piType,double *pdSoft);
struct base_particle *pTipsyParticle(TCTX ctx,unsigned iIndex,int *piType,double *pdSoft);
void TipsyReadAll(TCTX ctx);
void TipsyAddGas(TCTX ctx,struct gas_particle *pGas);
void TipsyAddDark(TCTX ctx,struct dark_particle *pDark);
void TipsyAddStar(TCTX ctx,struct star_particle *pStar);
void TipsyWriteAll(TCTX ctx,double dTime,char *pchFileName);
unsigned iTipsyNumParticles(TCTX ctx);
double dTipsyTime(TCTX ctx);


#endif








