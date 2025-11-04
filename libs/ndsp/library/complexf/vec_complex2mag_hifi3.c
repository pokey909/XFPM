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
  NatureDSP Signal Processing Library. Complex Math functions
    Complex magnitude
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/
/* Signal Processing Library API. */
#include "NatureDSP_Signal_complex.h"
#include "common.h"
#include "common_fpu.h"
#if (HAVE_VFPU==0 && HAVE_FPU==0)
DISCARD_FUN(void,vec_complex2mag,(float32_t  * restrict y, const complex_float  * restrict x, int N))
#else
/*-------------------------------------------------------------------------
  Complex magnitude
  Routines compute complex magnitude or its reciprocal

  Precision: 
  f     single precision floating point

  Input:
  x[N]  input complex data
  N     length of vector
  Output:
  y[N]  output data

  Restriction:
  none
-------------------------------------------------------------------------*/
void       vec_complex2mag (float32_t  * restrict y, const complex_float  * restrict x, int N)
#if HAVE_VFPU
{
  /*
    union ufloat32uint32 R, I, T, X, Z, TMP, U0, U1, T2;
    float32_t mnt_re, mnt_im;
    T.u = 0x7f800000;
    T2.u = 0x00400000;
    TMP.u = 0x7f000000;

    R.f = fabsf( crealf(x) );
    I.f = fabsf( cimagf(x) );

    if (isnan(R.f) || isnan(I.f)) return NAN;
    if (isinf(R.f) || isinf(I.f)) return T.f;


    U0.u = R.u&T.u;
    U1.u = I.u&T.u;
    X.u = MAX(U0.u, U1.u);
    T.u = (T.u - X.u );
    T.f = fminf(T.f, TMP.f); 
    mnt_re = R.f*T.f;
    mnt_im = I.f*T.f;
    Z.f = sqrtf(mnt_re*mnt_re + mnt_im*mnt_im);
    T.u = 0x00800000;
    X.u = (X.u - T.u);
    X.f = fmaxf(X.f, T2.f);
    Z.f = Z.f*X.f;
    return Z.f;
  */
  const xtfloatx2 * restrict X = (const xtfloatx2 *)x;
  xtfloatx2 * restrict Y = (xtfloatx2 *)y;
  int n;
  xtfloatx2 x0, x1, y0, z0, xre, xim, MAXEXP;
  xtfloat l;
  ae_int32x2 u0, u1, I0, QNAN;
  xtbool2 bnan, binf;
  ae_valign Y_va;
  if (N <= 0) return;
  Y_va = AE_ZALIGN64();
 
  for (n = 0; n<(N >> 1); n++)
  {
    MAXEXP = XT_AE_MOVXTFLOATX2_FROMINT32X2(0x7f000000);/*254 << 23*/
    QNAN = 0xffc00000;/* Unordored */
    I0 = 0x7f800000;

    XT_LSX2IP(x0, X, sizeof(complex_float));
    XT_LSX2IP(x1, X, sizeof(complex_float));

    x0 = XT_ABS_SX2(x0);
    x1 = XT_ABS_SX2(x1);
    xre = XT_SEL32_HH_SX2(x0, x1);
    xim = XT_SEL32_LL_SX2(x0, x1);
    u0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xre);
    u1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xim);
    bnan = XT_UN_SX2(xre, xim); /* is NaN */

    u0 = AE_AND32(u0, I0);
    u1 = AE_AND32(u1, I0);
    u1 = AE_MAX32(u0, u1);
    binf = AE_EQ32(u1, I0); /* is Inf */

    u0 = AE_SUB32(I0, u1);
    y0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u0);
    y0 = XT_MIN_SX2(y0, MAXEXP);
    xre = XT_MUL_SX2(xre, y0);
    xim = XT_MUL_SX2(xim, y0);

    x0 = XT_MUL_SX2(xre, xre);
    x1 = XT_MUL_SX2(xim, xim);
    x0 = XT_ADD_SX2(x0, x1);


    z0 = XT_SQRT_SX2(x0);
    
    u0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(z0);
    AE_MOVT32X2(u0, I0, binf);
    AE_MOVT32X2(u0, QNAN, bnan);
    z0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u0);
    u1 = AE_SUB32(u1, 0x00800000);
    y0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u1);
    y0 = XT_MAX_SX2(y0, XT_AE_MOVXTFLOATX2_FROMINT32X2(0x00400000));
    y0 = XT_MUL_SX2(y0, z0);
    XT_SASX2IP(y0, Y_va, Y);
  }
  XT_SASX2POSFP(Y_va, Y);
  if (N&1)
  {
    MAXEXP = XT_AE_MOVXTFLOATX2_FROMINT32X2(0x7f000000);/*254 << 23*/
    QNAN = 0xffc00000;/* Unordored */
    I0 = 0x7f800000;
    XT_LSX2IP(x0, X, sizeof(complex_float));
    
    x0 = XT_ABS_SX2(x0);
    xre = XT_SEL32_HH_SX2(x0, x0);
    xim = XT_SEL32_LL_SX2(x0, x0);
    u0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xre);
    u1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xim);
    bnan = XT_UN_SX2(xre, xim); /* is NaN */

    u0 = AE_AND32(u0, I0);
    u1 = AE_AND32(u1, I0);
    u1 = AE_MAX32(u0, u1);
    binf = AE_EQ32(u1, I0); /* is Inf */

    u0 = AE_SUB32(I0, u1);
    y0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u0);
    y0 = XT_MIN_SX2(y0, MAXEXP);
    xre = XT_MUL_SX2(xre, y0);
    xim = XT_MUL_SX2(xim, y0);

    x0 = XT_MUL_SX2(xre, xre);
    x1 = XT_MUL_SX2(xim, xim);
    x0 = XT_ADD_SX2(x0, x1);

    z0 = XT_SQRT_SX2(x0);

    u0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(z0);
    AE_MOVT32X2(u0, I0, binf);
    AE_MOVT32X2(u0, QNAN, bnan);
    z0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u0);
    u1 = AE_SUB32(u1, 0x00800000);
    y0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u1);
    y0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u1);
    y0 = XT_MAX_SX2(y0, XT_AE_MOVXTFLOATX2_FROMINT32X2(0x00400000));
    y0 = XT_MUL_SX2(y0, z0);
    l = XT_LOW_S(y0);
    XT_SSX(l, (xtfloat *)Y, 0);
  }
}
#else // scalar FPU
{
    const int blkSize = (MAX_ALLOCA_SZ/(sizeof(xtfloat)))/20;
    xtfloat scr[blkSize];
    int n,NN;
    const xtfloat* restrict pX=(const xtfloat*)x;
    const xtfloat* restrict pZrd;
          xtfloat* restrict pZwr;
    const xtfloat* restrict pYrd=(const xtfloat*)y;
          xtfloat* restrict pYwr=(      xtfloat*)y;
    NN=N;
    while(NN>0)
    {
        __Pragma("no_reorder")
        N=XT_MIN(NN,blkSize);
        pZwr=scr;
        pYwr=(      xtfloat*)y;
        for (n=0; n<N; n++)
        {
            xtfloat x0, x1, y0;
            xtfloat xre, xim;
            int t0, t1, nsa;
            int e0;
            int nsa0;

            XT_LSIP(xre,pX,1*sizeof(xtfloat));
            XT_LSIP(xim,pX,1*sizeof(xtfloat));
            xre = XT_ABS_S(xre);
            xim = XT_ABS_S(xim);
            t0 = XT_RFR(xre);
            t1 = XT_RFR(xim);
            nsa = XT_MAX(t0, t1);
            nsa = ((uint32_t)nsa)>> 23;
            nsa = (nsa-127);
            nsa = XT_MIN(nsa, 127);
            e0 = (127-nsa);
            nsa0 = (e0<<23);
            XT_MOVEQZ(nsa0,0x00400000,e0);
            y0 = XT_WFR(nsa0);

            xre = XT_MUL_S(xre, y0);
            xim = XT_MUL_S(xim, y0);

            x0 = XT_MUL_S(xre, xre);
            x1 = XT_MUL_S(xim, xim);

            x0 = XT_ADD_S(x0, x1);
            XT_SSIP(x0,pYwr,sizeof(xtfloat));

            e0 = (127+nsa);
            nsa0 = (e0<<23);
            XT_MOVEQZ(nsa0, 0x00400000, e0);
            x0 = XT_WFR(nsa0);
            XT_SSIP(x0,pZwr,sizeof(xtfloat));
        }
        __Pragma("no_reorder")
        pZrd=scr;
        pYrd=(const xtfloat*)y;
        pYwr=(      xtfloat*)y;
        for (n=0; n<N; n++)
        {
            xtfloat y0,x0;
            XT_LSIP(y0,pYrd,sizeof(xtfloat));
            y0 = XT_SQRT_S(y0);
            XT_LSIP(x0,pZrd,sizeof(xtfloat));
            y0 = XT_MUL_S(y0, x0);    
            XT_SSIP(y0,pYwr,sizeof(xtfloat));
        }
        NN-=N;
        y+=N;
    }
}
#endif
#endif
