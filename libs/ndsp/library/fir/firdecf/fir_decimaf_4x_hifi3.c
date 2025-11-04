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
    Real FIR Filter with decimation (4x)
    C code optimized for HiFi3
 */

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Filters and transformations */
#include "NatureDSP_Signal_fir.h"
#include "fir_decimaf_4x.h"
#include "common_fpu.h"

#if (HAVE_VFPU)

/*-------------------------------------------------------------------------
    4x decimator:
    Input/output:
    delay[M] - circular delay line
    Input:
    p        - pointer to delay line
    x[N*4]   - input signal
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
float32_t * fir_decimaf_4x(float32_t * restrict z, float32_t * restrict delay, float32_t * restrict p, const float32_t * restrict h, const float32_t * restrict x,  int M, int N)
{
          xtfloat   *restrict pz;
    const xtfloatx2 *restrict px;
    const xtfloatx2 *restrict ph;
          xtfloatx2 *pp;
    const xtfloatx2 *pd;
    ae_valign ax,ah,ad;
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

    pz=(xtfloat*)z;
    pp=(xtfloatx2*)p;
    for (n=0; n<(N>>1); n++,x+=8)
    {
        xtfloatx2 A0,B0,A1,B1,X0,X1,X2,X3,H0,H1;
        xtfloatx2 Y0,Y1,Y2,Y3,Y4;
        A0=B0=A1=B1=(xtfloatx2)0.0f;

        ph=(const xtfloatx2*)h;
        ah=AE_LA64_PP(ph);

        px=(xtfloatx2*)(x);
        ax=AE_LA64_PP(px);
        XT_LASX2IP(Y0,ax,px); 
        XT_LASX2IP(Y1,ax,px); 
        XT_LASX2IP(Y2,ax,px); 
        pd=(const xtfloatx2*)(((const float32_t*)pp));
        XT_LASX2NEGPC(ad,pd);
        XT_LASX2RIC(Y3,ad,pd);
        XT_LASX2RIC(Y4,ad,pd);
        X3=XT_SEL32_HL_SX2(Y1,Y0);
        X2=XT_SEL32_HL_SX2(Y2,Y1);
        X0=XT_SEL32_HL_SX2(Y0,Y3);
        X1=Y4;

        for (m=0; m<(M>>2); m++)
        {
            XT_LASX2IP(H0,ah,ph);
            XT_LASX2IP(H1,ah,ph);
            XT_MADD_SX2(A0,H0,X0);
            XT_MADD_SX2(A1,H0,X2);
            XT_MADD_SX2(B0,H1,X1);
            XT_MADD_SX2(B1,H1,X3);
            X2=X0;
            X3=X1;
            XT_LASX2RIC(X0,ad,pd);
            XT_LASX2RIC(X1,ad,pd);
        }
        if(M&2)
        {
            XT_LASX2IP(H0,ah,ph);
            XT_MADD_SX2(A0,H0,X0);
            XT_MADD_SX2(A1,H0,X2);
        }
        A0=A0+B0;
        A1=A1+B1;
        A0=A0+XT_SEL32_LH_SX2(A0,A0);
        A1=A1+XT_SEL32_LH_SX2(A1,A1);

        px=(xtfloatx2*)(x);
        ax=AE_LA64_PP(px);
        XT_LASX2IP(Y0,ax,px); 
        XT_LASX2IP(Y1,ax,px); 
        XT_LASX2IP(Y2,ax,px); 
        XT_LASX2IP(Y3,ax,px); 
        XT_SSX2XC(Y0,pp,sizeof(Y0));
        XT_SSX2XC(Y1,pp,sizeof(Y1));
        XT_SSX2XC(Y2,pp,sizeof(Y2));
        XT_SSX2XC(Y3,pp,sizeof(Y3));
        XT_SSIP(A0,pz,sizeof(float32_t));
        XT_SSIP(A1,pz,sizeof(float32_t));
    }
    return (float32_t*)pp;
} /* fir_decimaf_4x() */
#elif (HAVE_FPU)
// for scalar FPU
float32_t * fir_decimaf_4x(float32_t * restrict z, float32_t * restrict delay, float32_t * restrict p, const float32_t * restrict h, const float32_t * restrict x,  int M, int N)
{
  xtfloat* restrict pZ;
  const xtfloat *restrict pH;
  const xtfloat *restrict pX;
  xtfloat * pp;
  xtfloat* restrict pD;
  xtfloat x0, x1, x2, x3, x4, x5, x6, x7;
  xtfloat s0, s1, s2, s3, s4, s5, s6, s7;
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
  NASSERT(M%2==0);
  NASSERT(N%8==0);
  pD = (xtfloat*)p;
  /* set circular buffer boundaries */
  WUR_AE_CBEGIN0((uintptr_t)(delay + 0));
  WUR_AE_CEND0((uintptr_t)(delay + M));

  pp = (xtfloat*)p;
  pZ = (xtfloat*)z;
  pX = (const xtfloat*)x;

  /* process by 2 input samples */
  for (n = 0; n<(N >> 1); n++)
  {
    xtfloat A0, A1, A2, A3, xx;
    pH = (const xtfloat*)h;
    pD = (xtfloat*)pp;
    A0 = A1 = A2 = A3 = XT_CONST_S(0);
    XT_LSIP(x0, pX, 4);
    XT_LSIP(x1, pX, 4);
    XT_LSIP(x2, pX, 4);
    XT_LSIP(x3, pX, 4);
    XT_LSIP(x4, pX, 4);
    XT_LSIP(x5, pX, 4);
    XT_LSIP(x6, pX, 4);
    XT_LSIP(x7, pX, 4);
    s0 = x0;
    s1 = x1;
    s2 = x2;
    s3 = x3;
    s4 = x4;
    s5 = x5;
    s6 = x6;
    s7 = x7;
    { xtfloat dummy;  XT_LSXC(dummy, pD, -4); }
    { xtfloat dummy;  XT_LSXC(dummy, pD, -4); xx = dummy; }
    for (m = 0; m<M; m += 2)
    {
      h0 = pH[m + 0];
      h1 = pH[m + 1];

      XT_MADD_S(A0, x0, h0);
      XT_MADD_S(A1, x4, h0);
      XT_MADD_S(A2, xx, h1);
      XT_MADD_S(A3, x3, h1);

      x4 = x2;
      x3 = x1;
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
    XT_SSXC(s4, pp, 4);
    XT_SSXC(s5, pp, 4);
    XT_SSXC(s6, pp, 4);
    XT_SSXC(s7, pp, 4);

    XT_SSIP(A0, pZ, 4);
    XT_SSIP(A1, pZ, 4);
  }
  return (float32_t*)pp; 
}
#endif
