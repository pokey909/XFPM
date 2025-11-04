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
    Real interpolating FIR Filter
    C code optimized for HiFi3
    IntegrIT, 2006-2018
 */

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
#include "fir_interpf_Dx.h"
#include "common_fpu.h"

#if (HAVE_VFPU)

/*-------------------------------------------------------------------------
    Dx interpolator:
    Input/output:
    delay[M*D] - circular delay line
    Input:
    p        - pointer to delay line
    x[N]     - input signal
    h[M*D]   - impulse response
    N        - number of output samples
    Output:
    z[N*D]     - output samples
    Restrictions:
    N>0, M>0
    N   - multiple of 8
    M   - multiple of 4
    D > 1
    delay should be aligned on 8 byte boundary

    Returns:
    updated pointer to delay line
-------------------------------------------------------------------------*/
float32_t * fir_interpf_Dx(float32_t * restrict z, float32_t * restrict delay, float32_t * restrict p, const float32_t * restrict h, const float32_t * restrict x, int M, int D, int N)
{
#define P 8
  int n, m, j, _N;
  const xtfloatx2*  restrict pX = (const xtfloatx2*)x;
  const xtfloat  *  restrict px = (const xtfloat  *)x;
        xtfloat  * restrict pD = (      xtfloat  *)p;
  const xtfloat*   restrict pH = (const xtfloat*)h;
  xtfloat*          pZ = (xtfloat*)z;
  ae_valign x_align, d_align;
  NASSERT(x);
  NASSERT(z);
  NASSERT(N>0);
  NASSERT(M>0);
  NASSERT(D>1);
  NASSERT(N%8==0);
  NASSERT(M%4==0);

  _N = N&(~7);
  x_align = AE_LA64_PP(pX);
  WUR_AE_CBEGIN0((uintptr_t)(delay));
  WUR_AE_CEND0((uintptr_t)(delay + M));
  for (n = 0; n<_N; n += P)
  {
    xtfloatx2 x0, x1, x2, x3;
    xtfloatx2 acc0, acc1, acc2, acc3;
    xtfloat    s0, s1, s2, s3, s4, s5, s6, s7;
    __Pragma("loop_count min=5")
    for (j = 0; j<D; j++)
    {
      xtfloatx2 temp;
      pX = (xtfloatx2*)((uintptr_t)(x + n));
      x_align = AE_LA64_PP(pX);
      XT_LASX2IP(x0, x_align, pX);
      XT_LASX2IP(x1, x_align, pX);
      XT_LASX2IP(x2, x_align, pX);
      XT_LASX2IP(x3, x_align, pX);
      acc0 = 
      acc1 = 
      acc2 = 
      acc3 = (xtfloatx2)0.0f;
      d_align = AE_LA64_PP(pD);
      XT_LASX2RIC(temp, d_align, castxcc(const xtfloatx2,pD));

      x3 = XT_SEL32_LH_SX2(x3, x3);
      x2 = XT_SEL32_LH_SX2(x2, x2);
      x1 = XT_SEL32_LH_SX2(x1, x1);
      x0 = XT_SEL32_LH_SX2(x0, x0);

      __Pragma("loop_count min=1")
      for (m = 0; m<M; m += 2)
      {
        xtfloatx2 hm;
        xtfloatx2 t0, t1, t2;
        hm = pH[m + j*M + 0];
        XT_MADD_SX2(acc0, hm, x0);
        XT_MADD_SX2(acc1, hm, x1);
        XT_MADD_SX2(acc2, hm, x2);
        XT_MADD_SX2(acc3, hm, x3);
        /* select and load next sample */
        x3 = XT_SEL32_LH_SX2(x3, x2);
        t2 = XT_SEL32_LH_SX2(x2, x1);
        t1 = XT_SEL32_LH_SX2(x1, x0);
        XT_LASX2RIC(temp, d_align, castxcc(const xtfloatx2,pD));
        //temp = XT_SEL32_LH_SX2(temp, temp);
        t0 = XT_SEL32_LH_SX2(x0, temp);

        hm = pH[m + j*M + 1];
        XT_MADD_SX2(acc0, hm, t0);
        XT_MADD_SX2(acc1, hm, t1);
        XT_MADD_SX2(acc2, hm, t2);
        XT_MADD_SX2(acc3, hm, x3);
        /* select and load next sample */
        x3 = x2;
        x2 = x1;
        x1 = x0;
        x0 = temp;
      }
      XT_LASX2IC(temp, d_align, castxcc(const xtfloatx2,pD));
      pZ = (xtfloat*)((uintptr_t*)z + n*D + j);
      XT_SSXP((xtfloat)acc0, pZ, D * 4);
      acc0 = XT_SEL32_LH_SX2(acc0, acc0);
      XT_SSXP((xtfloat)acc0, pZ, D * 4);
      XT_SSXP((xtfloat)acc1, pZ, D * 4);
      acc1 = XT_SEL32_LH_SX2(acc1, acc1);
      XT_SSXP((xtfloat)acc1, pZ, D * 4);
      XT_SSXP((xtfloat)acc2, pZ, D * 4);
      acc2 = XT_SEL32_LH_SX2(acc2, acc2);
      XT_SSXP((xtfloat)acc2, pZ, D * 4);
      XT_SSXP((xtfloat)acc3, pZ, D * 4);
      acc3 = XT_SEL32_LH_SX2(acc3, acc3);
      XT_SSXP((xtfloat)acc3, pZ, -7 * D * 4 + 4);
    }
    XT_LSIP(s0, px, 4);
    XT_LSIP(s1, px, 4);
    XT_LSIP(s2, px, 4);
    XT_LSIP(s3, px, 4);
    XT_LSIP(s4, px, 4);
    XT_LSIP(s5, px, 4);
    XT_LSIP(s6, px, 4);
    XT_LSIP(s7, px, 4);
    XT_SSXC(s0, pD, 4);
    XT_SSXC(s1, pD, 4);
    XT_SSXC(s2, pD, 4);
    XT_SSXC(s3, pD, 4);
    XT_SSXC(s4, pD, 4);
    XT_SSXC(s5, pD, 4);
    XT_SSXC(s6, pD, 4);
    XT_SSXC(s7, pD, 4);
  }
  p = (float32_t*)pD;

  return (float32_t*)pD;
}/* fir_interpf_Dx() */
#elif HAVE_FPU
// for scalar FPU
float32_t * fir_interpf_Dx(float32_t * restrict z, float32_t * restrict delay, float32_t * restrict p, const float32_t * restrict h, const float32_t * restrict x, int M, int D, int N)
{
#define P 4
  int n, m, j, _N;
  const xtfloat*  restrict pX = (const xtfloat*)x;
  const xtfloat*  restrict px = (const xtfloat*)x;
  const xtfloat* restrict pD  = (const xtfloat*)p;
  const xtfloat*   restrict pH = (const xtfloat*)h;
  xtfloat*          pZ = (xtfloat*)z;
  ae_valign ax, az, ad;
  NASSERT(x);
  NASSERT(z);
  NASSERT(N>0);
  NASSERT(M>0);
  NASSERT(D>1);
  NASSERT(M % 4 == 0);
  NASSERT(N % 8 == 0);
  _N = N&(~3);
  ax = AE_LA64_PP(pX);
  az = AE_ZALIGN64();
  WUR_AE_CBEGIN0((uintptr_t)(delay));
  WUR_AE_CEND0((uintptr_t)(delay + M));
  for (n = 0; n<_N; n +=P)
  {
    xtfloat x0, x1, x2, x3, H0;
    xtfloat A0, A1, A2, A3;
    xtfloat s0, s1, s2, s3;
    pH = (const xtfloat*)h;
    __Pragma("loop_count min=5")
      for (j = 0; j<D; j++)
      {
        xtfloat temp;
        pX = (xtfloat*)((uintptr_t)(x + n));
        ax = AE_LA64_PP(pX);
        XT_LSIP(x0, castxcc(xtfloat, pX), 4);
        XT_LSIP(x1, castxcc(xtfloat, pX), 4);
        XT_LSIP(x2, castxcc(xtfloat, pX), 4);
        XT_LSIP(x3, castxcc(xtfloat, pX), 4);
        A0 = 
        A1 = 
        A2 = 
        A3 = XT_CONST_S(0);
        ad = AE_LA64_PP(pD);
        XT_LSXC(temp, castxcc(xtfloat, pD), -4);

        __Pragma("loop_count min=1")
        for (m = 0; m<M; m ++)
        {
          XT_LSIP(H0, pH, 4);
          XT_MADD_S(A0, H0, x0);
          XT_MADD_S(A1, H0, x1);
          XT_MADD_S(A2, H0, x2);
          XT_MADD_S(A3, H0, x3);
          x3=x2;x2=x1;x1=x0;
          XT_LSXC(x0, castxcc(xtfloat, pD), -4);
         
        }
        XT_LSXC(temp, castxcc(xtfloat, pD), 4);
        pZ = (xtfloat*)((uintptr_t*)z + n*D + j);
        XT_SSXP(A0, pZ, 4 * D);
        XT_SSXP(A1, pZ, 4 * D);
        XT_SSXP(A2, pZ, 4 * D);
        XT_SSXP(A3, pZ, -3 * 4 * D + 4);
      }
    
      XT_LSIP(s0, castxcc(xtfloat, px), 4);
      XT_LSIP(s1, castxcc(xtfloat, px), 4);
      XT_LSIP(s2, castxcc(xtfloat, px), 4);
      XT_LSIP(s3, castxcc(xtfloat, px), 4);
    
      XT_SSXC(s0, castxcc(xtfloat, pD), 4);
      XT_SSXC(s1, castxcc(xtfloat, pD), 4);
      XT_SSXC(s2, castxcc(xtfloat, pD), 4);
      XT_SSXC(s3, castxcc(xtfloat, pD), 4);
    }
  p = (float32_t*)pD;
  return (float32_t*)pD;
} 
#endif
