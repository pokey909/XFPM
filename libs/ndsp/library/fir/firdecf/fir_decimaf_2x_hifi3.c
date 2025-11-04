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
    Real FIR Filter with decimation (2x)
    C code optimized for HiFi3
 */

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Filters and transformations */
#include "NatureDSP_Signal_fir.h"
#include "fir_decimaf_2x.h"
#include "common_fpu.h"

#if (HAVE_VFPU)

/*-------------------------------------------------------------------------
    2x decimator:
    Input/output:
    delay[M] - circular delay line
    Input:
    p        - pointer to delay line
    x[N*2]   - input signal
    h[M]     - impulse response
    N        - number of output samples
    Output:
    z[N]     - output samples
    Restrictions:
    N>0, M>0
    M multiple of 2
    N multiple of 8
    delay should be aligned on 8 byte boundary

    Returns:
    updated pointer to delay line
-------------------------------------------------------------------------*/
float32_t * fir_decimaf_2x(float32_t * restrict z, float32_t * restrict delay, float32_t * restrict p, const float32_t * restrict h, const float32_t * restrict x,  int M, int N)
{
    xtfloat* restrict pz;
    const xtfloatx2 *restrict ph;
    const xtfloatx2 *restrict px;
          xtfloatx2 * pp;
    const xtfloatx2 * pd;
    ae_valign ax,ad,ah;
    int n,m;
    NASSERT(x);
    NASSERT(z);
    NASSERT(h);
    NASSERT(delay);
    NASSERT(p);
    NASSERT_ALIGN(delay,8);
    NASSERT(N>0);
    NASSERT(M>0);
    NASSERT(M%2==0);
    NASSERT(N%8==0);

    /* set circular buffer boundaries */
    WUR_AE_CBEGIN0( (uintptr_t)( delay + 0 ) );
    WUR_AE_CEND0  ( (uintptr_t)( delay + M ) );

    pp=(xtfloatx2*)p;
    pz=(xtfloat*)z;
    /* process by 8 input samples */
    for (n=0; n<(N>>2); n++,x+=8)
    {
        xtfloatx2 A0,A1,A2,A3,XX,X0,X1,X2,X3,H0;

        ph=(const xtfloatx2*)h;
        ah=AE_LA64_PP(ph);

        A0=A1=A2=A3=(xtfloatx2)0.0f;
        px=(xtfloatx2*)(x);
        ax=AE_LA64_PP(px);
        XT_LASX2IP(X0,ax,px);
        XT_LASX2IP(X1,ax,px);
        XT_LASX2IP(X2,ax,px);
        XT_LASX2IP(X3,ax,px);
        pd=(const xtfloatx2*)(((const float32_t*)pp));
        XT_LASX2NEGPC(ad,pd);
        XT_LASX2RIC(XX,ad,pd);

        X3=XT_SEL32_HL_SX2(X3,X2);
        X2=XT_SEL32_HL_SX2(X2,X1);
        X1=XT_SEL32_HL_SX2(X1,X0);
        X0=XT_SEL32_HL_SX2(X0,XX);
        for (m=0; m<M; m+=2)
        {
            XT_LASX2IP(H0,ah,ph);
            XT_MADD_SX2(A0,H0,X0);
            XT_MADD_SX2(A1,H0,X1);
            XT_MADD_SX2(A2,H0,X2);
            XT_MADD_SX2(A3,H0,X3);
            X3=X2;
            X2=X1;
            X1=X0;
            XT_LASX2RIC(X0,ad,pd);
        }
        A0=A0+XT_SEL32_LH_SX2(A0,A0);
        A1=A1+XT_SEL32_LH_SX2(A1,A1);
        A2=A2+XT_SEL32_LH_SX2(A2,A2);
        A3=A3+XT_SEL32_LH_SX2(A3,A3);
        px=(xtfloatx2*)(x);
        ax=AE_LA64_PP(px);
        XT_LASX2IP(X0,ax,px);
        XT_LASX2IP(X1,ax,px);
        XT_LASX2IP(X2,ax,px);
        XT_LASX2IP(X3,ax,px);
        XT_SSX2XC(X0,pp,sizeof(X0));
        XT_SSX2XC(X1,pp,sizeof(X1));
        XT_SSX2XC(X2,pp,sizeof(X2));
        XT_SSX2XC(X3,pp,sizeof(X3));

        XT_SSIP(A0,pz,sizeof(float32_t));
        XT_SSIP(A1,pz,sizeof(float32_t));
        XT_SSIP(A2,pz,sizeof(float32_t));
        XT_SSIP(A3,pz,sizeof(float32_t));
    }

    return (float32_t*)pp;
} /* fir_decimaf_2x() */
#elif (HAVE_FPU)
// for scalar FPU
float32_t * fir_decimaf_2x(float32_t * restrict z, float32_t * restrict delay, float32_t * restrict p, const float32_t * restrict h, const float32_t * restrict x,  int M, int N)
{
  xtfloat* restrict pZ;
  const xtfloat *restrict pH;
  const xtfloat *restrict pX;
  xtfloat * pp;
  xtfloat* restrict pD;
  xtfloat x0, x1, x2, x3;
  xtfloat s0, s1, s2, s3;
  xtfloat h0, h1;
  int n, m;
  NASSERT(x);
  NASSERT(z);
  NASSERT(h);
  NASSERT(delay);
  NASSERT(p);
  NASSERT_ALIGN(delay, 8);
  NASSERT(N>0);
  NASSERT(M>0);
  NASSERT(M % 2 == 0);
  NASSERT(N % 8 == 0);
  pD = (xtfloat*)p;
  /* set circular buffer boundaries */
  WUR_AE_CBEGIN0((uintptr_t)(delay + 0));
  WUR_AE_CEND0((uintptr_t)(delay + M));

  pp = (xtfloat*)p;
  pZ = (xtfloat*)z;
  pX = (const xtfloat*)x;
  
  /* process by 2 input samples */
  for (n=0; n<(N>>1); n++)
  {
    xtfloat A0, A1, A2, A3, xx;
    pH = (const xtfloat*)h;
    pD = (xtfloat*)pp;
    A0=A1=A2=A3=XT_CONST_S(0);
    XT_LSIP(x0, pX, 4);
    XT_LSIP(x1, pX, 4);
    XT_LSIP(x2, pX, 4);
    XT_LSIP(x3, pX, 4);
    s0 = x0;
    s1 = x1;
    s2 = x2;
    s3 = x3;
    { xtfloat dummy;  XT_LSXC(dummy, pD, -4); xx = dummy;}
    { xtfloat dummy;  XT_LSXC(dummy, pD, -4); xx = dummy; }
    for (m=0; m<M; m+=2)
    {
      h0 = pH[m+0];
      h1 = pH[m+1];

      XT_MADD_S(A0, x0, h0);
      XT_MADD_S(A1, x2, h0);
      XT_MADD_S(A2, xx, h1);
      XT_MADD_S(A3, x1, h1);

      x3 = x2;
      x2 = x0;
      x1 = xx;
      XT_LSXC(x0, pD, -4);
      XT_LSXC(xx, pD, -4);
    }
    A0 = XT_ADD_S(A0, A2);
    A1 = XT_ADD_S(A1, A3);
    { xtfloat dummy;  XT_LSXC(dummy, pD, 4); }
    XT_SSXC(s0, pp, 4);
    XT_SSXC(s1, pp, 4);
    XT_SSXC(s2, pp, 4);
    XT_SSXC(s3, pp, 4);

    XT_SSIP(A0, pZ, 4);
    XT_SSIP(A1, pZ, 4);
  }
  return (float32_t*)pp;
}

#endif
