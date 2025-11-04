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
#include "inff_tbl.h"
#if (HAVE_VFPU==0 && HAVE_FPU==0)
DISCARD_FUN(void,vec_complex2invmag,(float32_t  * restrict y, const complex_float  * restrict x, int N))
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
#define sz_f32    (int)sizeof(float32_t)
void       vec_complex2invmag (float32_t  * restrict y, const complex_float  * restrict x, int N)
#if HAVE_VFPU
{
  /*
  union ufloat32uint32 R, I, T, X, Z, TMP, U0, U1;
  float32_t mnt_re, mnt_im;
  T.u = 0x7f800000;
  TMP.u = 0x7f000000;

  R.f = fabsf( crealf(x) );
  I.f = fabsf( cimagf(x) );

  if (isnan(R.f) || isnan(I.f)) return NAN;
  if (isinf(R.f) || isinf(I.f)) return 0;


  U0.u = R.u&T.u;
  U1.u = I.u&T.u;
  X.u = MAX(U0.u, U1.u);
  T.u = (T.u - X.u );
  T.f = fminf(T.f, TMP.f); 
  mnt_re = R.f*T.f;
  mnt_im = I.f*T.f;
  Z.f = 1.f / sqrtf(mnt_re*mnt_re + mnt_im*mnt_im);
  Z.f = Z.f*T.f;
  return Z.f;
  */
  #if 1
  const xtfloatx2 * restrict X = (const xtfloatx2 *)x;
        xtfloatx2 * restrict Y = (xtfloatx2 *)y;
  int n;
  xtfloatx2 x0, x1, y0, z0, xre, xim, MAXEXP, _0;
  xtfloat l;
  ae_int32x2 u0, u1, I0, QNAN;
  xtbool2 bnan, binf, bzero;
  ae_valign Y_va;
  if (N <= 0) return;

  Y_va = AE_ZALIGN64();
  MAXEXP = XT_AE_MOVXTFLOATX2_FROMINT32X2(0x7f000000);/*254 << 23*/
  QNAN = 0xffc00000;/* Unordored */
  I0 = 0x7f800000;
  _0 = (xtfloatx2)0.0f;
  for (n = 0; n<(N >> 2); n++) 
  {
    XT_LSX2IP(x0, X, sizeof(complex_float));
    XT_LSX2IP(x1, X, sizeof(complex_float));

    MAXEXP = XT_AE_MOVXTFLOATX2_FROMINT32X2(0x7f000000);/*254 << 23*/
    QNAN = 0xffc00000;/* Unordored */
    I0 = 0x7f800000;
    _0 = (xtfloatx2)0.0f;

    x0 = XT_ABS_SX2(x0); 
    x1 = XT_ABS_SX2(x1);
    xre = XT_SEL32_HH_SX2(x0, x1);
    xim = XT_SEL32_LL_SX2(x0, x1);
    u0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xre);
    u1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xim);
    bnan = XT_UN_SX2(xre, xim); /* is NaN */

    u0 = AE_AND32(u0, I0);
    u1 = AE_AND32(u1, I0);
    u0 = AE_MAX32(u0, u1);
    binf = AE_EQ32(u0, I0); /* is Inf */

    u0 = AE_SUB32(I0, u0);
    y0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u0);
    y0 = XT_MIN_SX2(y0, MAXEXP);
    xre = XT_MUL_SX2(xre, y0);
    xim = XT_MUL_SX2(xim, y0);

    x0 = XT_MUL_SX2(xre, xre);
    x1 = XT_MUL_SX2(xim, xim);
    x0 = XT_ADD_SX2(x0, x1);
   
    z0 = XT_RSQRT_SX2(x0);

    bzero = XT_OEQ_SX2(x0, _0); /* is zero */
    u0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(z0);
    AE_MOVT32X2(u0, I0, bzero);
    AE_MOVT32X2(u0, AE_ZERO32(), binf);
    AE_MOVT32X2(u0, QNAN, bnan);
    z0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u0);
    y0 = XT_MUL_SX2(y0, z0);
    XT_SASX2IP(y0, Y_va, Y);


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
    u0 = AE_MAX32(u0, u1);
    binf = AE_EQ32(u0, I0); /* is Inf */

    u0 = AE_SUB32(I0, u0);
    y0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u0);
    y0 = XT_MIN_SX2(y0, MAXEXP);
    xre = XT_MUL_SX2(xre, y0);
    xim = XT_MUL_SX2(xim, y0);

    x0 = XT_MUL_SX2(xre, xre);
    x1 = XT_MUL_SX2(xim, xim);
    x0 = XT_ADD_SX2(x0, x1);

    z0 = XT_RSQRT_SX2(x0);

    bzero = XT_OEQ_SX2(x0, _0); /* is zero */
    u0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(z0);
    AE_MOVT32X2(u0, I0, bzero);
    AE_MOVT32X2(u0, AE_ZERO32(), binf);
    AE_MOVT32X2(u0, QNAN, bnan);
    z0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u0);
    y0 = XT_MUL_SX2(y0, z0);
    XT_SASX2IP(y0, Y_va, Y);
  }
  if (N&2)
  {
    XT_LSX2IP(x0, X, sizeof(complex_float));
    XT_LSX2IP(x1, X, sizeof(complex_float));

    MAXEXP = XT_AE_MOVXTFLOATX2_FROMINT32X2(0x7f000000);/*254 << 23*/
    QNAN = 0xffc00000;/* Unordored */
    I0 = 0x7f800000;
    _0 = (xtfloatx2)0.0f;

    x0 = XT_ABS_SX2(x0); 
    x1 = XT_ABS_SX2(x1);
    xre = XT_SEL32_HH_SX2(x0, x1);
    xim = XT_SEL32_LL_SX2(x0, x1);
    u0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xre);
    u1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xim);
    bnan = XT_UN_SX2(xre, xim); /* is NaN */

    u0 = AE_AND32(u0, I0);
    u1 = AE_AND32(u1, I0);
    u0 = AE_MAX32(u0, u1);
    binf = AE_EQ32(u0, I0); /* is Inf */

    u0 = AE_SUB32(I0, u0);
    y0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u0);
    y0 = XT_MIN_SX2(y0, MAXEXP);
    xre = XT_MUL_SX2(xre, y0);
    xim = XT_MUL_SX2(xim, y0);

    x0 = XT_MUL_SX2(xre, xre);
    x1 = XT_MUL_SX2(xim, xim);
    x0 = XT_ADD_SX2(x0, x1);

    z0 = XT_RSQRT_SX2(x0);

    bzero = XT_OEQ_SX2(x0, _0); /* is zero */
    u0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(z0);
    AE_MOVT32X2(u0, I0, bzero);
    AE_MOVT32X2(u0, AE_ZERO32(), binf);
    AE_MOVT32X2(u0, QNAN, bnan);
    z0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u0);
    y0 = XT_MUL_SX2(y0, z0);
    XT_SASX2IP(y0, Y_va, Y);

  }
  XT_SASX2POSFP(Y_va, Y);

  if (N&1)
  {
    MAXEXP = XT_AE_MOVXTFLOATX2_FROMINT32X2(0x7f000000);/*254 << 23*/
    QNAN = 0xffc00000;/* Unordored */
    I0 = 0x7f800000;
    _0 = (xtfloatx2)0.0f;
    XT_LSX2IP(x0, X, sizeof(complex_float));

    x0 = XT_ABS_SX2(x0); 
    xre = XT_SEL32_HH_SX2(x0, x0);
    xim = XT_SEL32_LL_SX2(x0, x0);
    u0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xre);
    u1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xim);
    bnan = XT_UN_SX2(xre, xim); /* is NaN */

    u0 = AE_AND32(u0, I0);
    u1 = AE_AND32(u1, I0);
    u0 = AE_MAX32(u0, u1);
    binf = AE_EQ32(u0, I0); /* is Inf */

    u0 = AE_SUB32(I0, u0);
    y0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u0);
    y0 = XT_MIN_SX2(y0, MAXEXP);
    xre = XT_MUL_SX2(xre, y0);
    xim = XT_MUL_SX2(xim, y0);

    x0 = XT_MUL_SX2(xre, xre);
    x1 = XT_MUL_SX2(xim, xim);
    x0 = XT_ADD_SX2(x0, x1);

    z0 = XT_RSQRT_SX2(x0);

    bzero = XT_OEQ_SX2(x0, _0); /* is zero */
    u0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(z0);
    AE_MOVT32X2(u0, I0, bzero);
    AE_MOVT32X2(u0, AE_ZERO32(), binf);
    AE_MOVT32X2(u0, QNAN, bnan);
    z0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u0);
    y0 = XT_MUL_SX2(y0, z0);
    l = XT_LOW_S(y0);
    XT_SSX(l, (xtfloat *)Y, 0);
  }

  #else
  const xtfloatx2 * restrict X = (const xtfloatx2 *)x;
  xtfloatx2 * restrict Y = (xtfloatx2 *)y;
  int n;
  
  xtfloatx2 x0, x1, y0, xre, xim, I0, xabs;
  xtfloat l;
  ae_int32x2 t0, t1, nsa, u0, u1;
  ae_int32x2 e0, e1;
  ae_int32x2 nsa0;
  xtbool2 b0, b1;
  ae_valign Y_va;
  if (N <= 0) return;
  Y_va = AE_ZALIGN64();
  I0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(0x7f800000);
  for (n = 0; n<(N >> 1); n++)
  {
    XT_LSX2IP(x0, X, sizeof(complex_float));
    XT_LSX2IP(x1, X, sizeof(complex_float));
    x0 = XT_ABS_SX2(x0);
    x1 = XT_ABS_SX2(x1);

    xre = XT_SEL32_HH_SX2(x0, x1);
    xim = XT_SEL32_LL_SX2(x0, x1);

    t0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xre);
    t1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xim);

    u0 = AE_SRLI32(t0, 23);
    u1 = AE_SRLI32(t1, 23);
    u0 = AE_AND32(u0, 0xFF);
    u1 = AE_AND32(u1, 0xFF);
    nsa = AE_MAX32(u0, u1);
    nsa = AE_SUB32(nsa,127);

    nsa = AE_MIN32(nsa,127);

    e0 = AE_SUB32(127, nsa);
    b0 = AE_EQ32(e0, 0);
    nsa0 = AE_SLLI32(e0, 23);
    AE_MOVT32X2(nsa0, 0x00400000, b0);
    y0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(nsa0);
    xre = XT_MUL_SX2(xre, y0);
    xim = XT_MUL_SX2(xim, y0);

    x0 = XT_MUL_SX2(xre, xre);
    x1 = XT_MUL_SX2(xim, xim);
    x0 = XT_ADD_SX2(x0, x1);
    y0 = XT_RSQRT_SX2(x0);
    b0 = XT_OEQ_SX2(x0, 0.0f);
    b1 = XT_OEQ_SX2(x0, I0);
    t0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(y0);
    AE_MOVT32X2(t0, 0x7f800000, b0);
    AE_MOVT32X2(t0, 0x0, b1);
    y0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(t0);

    x0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(nsa0);
    y0 = XT_MUL_SX2(y0, x0);

    XT_SASX2IP(y0, Y_va, Y);
  }
  XT_SASX2POSFP(Y_va, Y);
  if (N & 1)
  {
    XT_LSX2IP(x0, X, sizeof(complex_float));
    x0 = XT_ABS_SX2(x0);

    xre = XT_SEL32_HH_SX2(x0, x0);
    xim = XT_SEL32_LL_SX2(x0, x0);

    t0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xre);
    t1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xim);

    u0 = AE_SRLI32(t0, 23);
    u1 = AE_SRLI32(t1, 23);
    u0 = AE_AND32(u0, 0xFF);
    u1 = AE_AND32(u1, 0xFF);
    e0 = AE_SUB32(u0, 127);
    e1 = AE_SUB32(u1, 127);
    nsa = AE_MAX32(e0, e1);
    b0 = AE_LT32(127, nsa);
    AE_MOVT32X2(nsa, 127, b0);

    e0 = AE_SUB32(127, nsa);
    b0 = AE_EQ32(e0, 0);
    nsa0 = AE_SLLI32(e0, 23);
    AE_MOVT32X2(nsa0, 0x00400000, b0);
    y0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(nsa0);
    xre = XT_MUL_SX2(xre, y0);
    xim = XT_MUL_SX2(xim, y0);

    x0 = XT_MUL_SX2(xre, xre);
    x1 = XT_MUL_SX2(xim, xim);
    x0 = XT_ADD_SX2(x0, x1);
    y0 = XT_RSQRT_SX2(x0);
    b0 = XT_OEQ_SX2(x0, 0.0f);
    b1 = XT_OEQ_SX2(x0, I0);
    t0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(y0);
    AE_MOVT32X2(t0, 0x7f800000, b0);
    AE_MOVT32X2(t0, 0x0, b1);
    y0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(t0);
    x0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(nsa0);
    y0 = XT_MUL_SX2(y0, x0);
    l = XT_LOW_S(y0);
    XT_SSX(l, (xtfloat *)Y, 0);
  }
  #endif
}
#else   // scalar FPU
{
    int n;
    const xtfloat* restrict pX=(const xtfloat*)x;
          xtfloat* restrict pY=(      xtfloat*)y;
    if (N<=0) return;
    for (n=0; n<N; n++)
    {
        xtfloat x0, x1, y0,coef;
        xtfloat xre, xim;
        int t0, t1, nsa;
        int e0;
        int nsa0;
        xtbool b0, b1;

        XT_LSIP(xre,pX,1*sizeof(xtfloat));
        XT_LSIP(xim,pX,1*sizeof(xtfloat));
        xre = XT_ABS_S(xre);
        xim = XT_ABS_S(xim);

        t0 = XT_RFR(xre);
        t1 = XT_RFR(xim);

        nsa = XT_MAX(t0, t1);
        nsa = ((uint32_t)nsa)>>23;
        nsa = nsa-127;
        nsa = XT_MIN(nsa, 127);
        e0 = (127-nsa);
        nsa0 = e0<<23;
        XT_MOVEQZ(nsa0,0x00400000,e0);
        coef = XT_WFR(nsa0);

        xre = XT_MUL_S(xre, coef);
        xim = XT_MUL_S(xim, coef);

        x0 = XT_MUL_S(xre, xre);
        x1 = XT_MUL_S(xim, xim);

        x0 = XT_ADD_S(x0, x1);

        y0 = XT_RSQRT_S(x0);
        b0 = XT_OEQ_S(x0, XT_CONST_S(0));
        b1 = XT_OEQ_S(x0, plusInff.f);
        XT_MOVT_S(y0,plusInff.f,b0);
        XT_MOVT_S(y0,XT_CONST_S(0),b1);
        y0 = XT_MUL_S(y0, coef);
        XT_SSIP(y0,pY,sizeof(xtfloat));
    }
}
#endif
#endif
