/* ------------------------------------------------------------------------ */
/* Copyright (c) 2021-2025 Cadence Design Systems, Inc. ALL RIGHTS RESERVED.*/
/* These coded instructions, statements, and computer programs ('Cadence    */
/* Libraries') are the copyrighted works of Cadence Design Systems Inc.     */
/* Cadence IP is licensed for use with Cadence processor cores only and     */
/* must not be used for any other processors and platforms. Your use of the */
/* Cadence Libraries is subject to the terms of the license agreement you   */
/* have entered into with Cadence Design Systems, or a sublicense granted   */
/* to you by a direct Cadence license.                                      */
/* ------------------------------------------------------------------------ */
/*  IntegrIT, Ltd.   www.integrIT.com, info@integrIT.com                    */
/*                                                                          */
/* NatureDSP Signal Library for HiFi 3 and 3z DSP                                  */
/*                                                                          */
/* This library contains copyrighted materials, trade secrets and other     */
/* proprietary information of IntegrIT, Ltd. This software is licensed for  */
/* use with Cadence processor cores only and must not be used for any other */
/* processors and platforms. The license to use these sources was given to  */
/* Cadence, Inc. under Terms and Condition of a Software License Agreement  */
/* between Cadence, Inc. and IntegrIT, Ltd.                                 */
/* ------------------------------------------------------------------------ */
/*          Copyright (c) 2009-2021 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
    Real FIR Filter with decimation
    C code optimized for HiFi3
 */

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
#include "fir_decimaf_Dx.h"
#include "common_fpu.h"

#if (HAVE_VFPU)

/*-------------------------------------------------------------------------
    universal decimator:
    Input/output:
    delay[M] - circular delay line
    Input:
    p        - pointer to delay line
    x[N*D]   - input signal
    h[M]     - impulse response
    D          decimation factor
    N        - number of output samples
    Output:
    z[N]     - output samples
    Restrictions:
    N>0, M>0
    D>4
    M multiple of 2
    N multiple of 8
    delay should be aligned on 8 byte boundary

    Returns:
    updated pointer to delay line
-------------------------------------------------------------------------*/
float32_t * fir_decimaf_Dx(float32_t * restrict z, float32_t * restrict delay, float32_t * restrict p, const float32_t * restrict h, const float32_t * restrict x,  int M, int D, int N)
{
    const xtfloat   *px;
          xtfloat   *pz;
    const xtfloatx2 *ph;
          xtfloat   *pp;
    const xtfloatx2 *pd;
    xtfloatx2 A0,A1,A2,A3,XX,X0,X1,X2,X3,H0,H1,H2,H3;
    ae_valign ah,ad;

    int m,n,d;
    NASSERT(N>0 && M>0);
    NASSERT(M%2==0);
    NASSERT(N%8==0);
    NASSERT(D>4);
    NASSERT_ALIGN(delay,8);
    WUR_AE_CBEGIN0( (uintptr_t)( delay + 0 ) );
    WUR_AE_CEND0  ( (uintptr_t)( delay + M ) );
    px=(const xtfloat*)x;
    pz=(      xtfloat*)z;
    pp=(xtfloat*)p;
    for (n=0; n<N; n+=1)
    {
        ph=(const xtfloatx2*)h;
        ah=AE_LA64_PP(ph);

        A0=A1=A2=A3=(xtfloatx2)0.0f;
        pd=(const xtfloatx2*)pp;
        XX=XT_LSI(px,0);
        XT_LASX2NEGPC(ad,pd);
        XT_LASX2RIC(X0,ad,pd);
        XT_LASX2RIC(X1,ad,pd);
        XT_LASX2RIC(X2,ad,pd);
        XT_LASX2RIC(X3,ad,pd);
        X0=XT_SEL32_HL_SX2(XX,X0);
        for (m=0; m<(M&~7); m+=8)
        {
            XT_LASX2IP(H0,ah,ph);
            XT_LASX2IP(H1,ah,ph);
            XT_LASX2IP(H2,ah,ph);
            XT_LASX2IP(H3,ah,ph);
            XT_MADD_SX2(A0,H0,X0);
            XT_MADD_SX2(A1,H1,X1);
            XT_MADD_SX2(A2,H2,X2);
            XT_MADD_SX2(A3,H3,X3);
            XT_LASX2RIC(X0,ad,pd);
            XT_LASX2RIC(X1,ad,pd);
            XT_LASX2RIC(X2,ad,pd);
            XT_LASX2RIC(X3,ad,pd);
        }
        if(M&4)
        {
            XT_LASX2IP(H0,ah,ph);
            XT_LASX2IP(H1,ah,ph);
            XT_MADD_SX2(A0,H0,X0);
            XT_MADD_SX2(A1,H1,X1);
            X0=X2;
        }
        if(M&2)
        {
            XT_LASX2IP(H0,ah,ph);
            XT_MADD_SX2(A0,H0,X0);
        }
        A0=A0+A1;
        A2=A2+A3;
        A0=A0+A2;
        A0=A0+XT_SEL32_LH_SX2(A0,A0);
        __Pragma("loop_count min=4");
        for(d=0; d<D; d++) 
        {
            xtfloat t;
            XT_LSIP(t,px,sizeof(float32_t));
            XT_SSXC(t,pp,sizeof(float32_t));
        }
        XT_SSIP(A0,pz,sizeof(float32_t));
    }
    return (float32_t*)pp;
} /* fir_decimaf_Dx() */
#elif (HAVE_FPU)
// for scalar FPU
float32_t * fir_decimaf_Dx(float32_t * restrict z, float32_t * restrict delay, float32_t * restrict p, const float32_t * restrict h, const float32_t * restrict x,  int M, int D, int N)
{
    const xtfloat * restrict pX=(const xtfloat *)x;
          xtfloat * restrict pZ=(      xtfloat *)z;
          xtfloat * restrict pD;
          xtfloat * restrict pP=(      xtfloat *)p;
    int k,n,m,j;

    NASSERT_ALIGN(x,8);
    NASSERT(N%8==0);
    NASSERT(x);
    NASSERT(z);
    if(N<=0) return p;
    NASSERT(N>0);
    NASSERT(M>0);
    NASSERT(D>1);
    NASSERT(M%2==0);
        /* set circular buffer boundaries */
    WUR_AE_CBEGIN0((uintptr_t)(delay + 0));
    WUR_AE_CEND0  ((uintptr_t)(delay + M));
    for (k=0,n=0; n<N*D; n+=D,k++)
    {
        xtfloat a0,a1,xn;
        a0=a1=XT_CONST_S(0);
        pD=(xtfloat *)pP;
        XT_LSIP(xn,pX,sizeof(xtfloat));
        XT_SSXC(xn,pP,sizeof(xtfloat));
        for (m=0; m<M; m+=2)
        {
            xtfloat d0,d1;
            XT_LSXC(d0,pD,-(int)sizeof(xtfloat));
            XT_LSXC(d1,pD,-(int)sizeof(xtfloat));
            XT_MADD_S(a0,h[m],d0);
            XT_MADD_S(a1,h[m+1],d1);
        }
        __Pragma("no_reorder")
        for(j=0; j<D-1; j++) 
        {
            XT_LSIP(xn,pX,sizeof(xtfloat));
            XT_SSXC(xn,pP,sizeof(xtfloat));
        }
        a0=XT_ADD_S(a0,a1);
        XT_SSIP(a0,pZ,sizeof(xtfloat));
    }
    return (float32_t*)pP;
}
#endif
