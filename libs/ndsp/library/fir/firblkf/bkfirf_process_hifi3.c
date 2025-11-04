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
  NatureDSP Signal Processing Library. FIR part
    Real block FIR filter, floating point
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Real FIR filter.
  Computes a real FIR filter (direct-form) using IR stored in vector h. The 
  real data input is stored in vector x. The filter output result is stored 
  in vector y. The filter calculates N output samples using M coefficients 
  and requires last M-1 samples in the delay line.
  NOTE: 
  1. User application is not responsible for management of delay lines
  2. User has an option to set IR externally or copy from original location 
     (i.e. from the slower constant memory). In the first case, user is 
     responsible for right alignment, ordering and zero padding of filter 
     coefficients - usually array is composed from zeroes (left padding), 
     reverted IR and right zero padding.


  Precision: 
  16x16    16-bit data, 16-bit coefficients, 16-bit outputs
  24x24    24-bit data, 24-bit coefficients, 24-bit outputs
  24x24p   use 24-bit data packing for internal delay line buffer
           and internal coefficients storage
  32x16    32-bit data, 16-bit coefficients, 32-bit outputs
  32x32    32-bit data, 32-bit coefficients, 32-bit outputs
  f        floating point

  Input:
  x[N]     input samples, Q31, Q15, floating point
  h[M]     filter coefficients in normal order, Q31, Q15, floating point
  N        length of sample block, should be a multiple of 4
  M        length of filter, should be a multiple of 4
  extIR    if zero, IR is copied from original location, otherwise not
           but user should keep alignment, order of coefficients 
           and zero padding requirements shown below
  Output:
  y[N]     output samples, Q31, Q15, floating point

  Alignment, ordering and zero padding for external IR  (extIR!=0)
  ------------------------+----------+--------------+--------------+----------------
  Function                |Alignment,|Left zero     |   Coefficient| Right zero 
                          | bytes    |padding, bytes|   order      | padding, bytes
  ------------------------+----------+--------------+--------------+----------------
  bkfir16x16_init         |     8    |      2       |  inverted    |  6
  bkfir24x24_init         |     8    |      4       |  inverted    |  12
  bkfir24x24p_init        |     8    |((-M&4)+1)*3  |  inverted    |  7
  bkfir32x16_init (M>32)  |     8    |     10       |  inverted    |  6
  bkfir32x16_init (M<=32) |     8    |      2       |  inverted    |  6
  bkfir32x32_init         |     8    |      4       |  inverted    |  12
  bkfirf_init             |     8    |      0       |  direct      |  0
  ------------------------+----------+--------------+--------------+----------------

  Restrictions:
  x, y     should not be overlapping
  x, h     aligned on a 8-bytes boundary
  N, M     multiples of 4 
-------------------------------------------------------------------------*/
/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "bkfirf.h"
#include "common_fpu.h"

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void,bkfirf_process,( bkfirf_handle_t _bkfir, 
                         float32_t * restrict  y,
                   const float32_t * restrict  x, int N ))
#elif HAVE_VFPU
/* process block of samples */
void bkfirf_process( bkfirf_handle_t _bkfir, 
                         float32_t * restrict  y,
                   const float32_t * restrict  x, int N )
{
    ae_valign yalign;
    bkfirf_t *state;
    int n, m, M, _N;
    xtfloatx2 h0;
    xtfloatx2 x0, x1, x2, x3;
    xtfloatx2 s0, s1, s2, s3;
    xtfloatx2 acc0, acc1, acc2, acc3;
    const xtfloatx2* restrict pX = (const xtfloatx2*)x;
          xtfloatx2* restrict pD;
    const xtfloat*   restrict pH;
          xtfloatx2*          pZ = (      xtfloatx2*)y;
    NASSERT(_bkfir);
    state=(bkfirf_t*)_bkfir;
    NASSERT(state->h);
    NASSERT(state->d);
    NASSERT(state->p);
    NASSERT_ALIGN(state->h,8);
    NASSERT_ALIGN(state->d,8);
    NASSERT_ALIGN(state->p,8);
    NASSERT(N%4==0);
    NASSERT_ALIGN(x,8);
    NASSERT((state->M%4)==0);
    NASSERT(x);
    NASSERT(y);
    if(N<=0) return;
    pD = (      xtfloatx2*)state->p;
    pH = (const xtfloat*  )state->h;

    M=state->M;
    NASSERT(N>0);
    NASSERT(M>0);
#define P1 8    /* unroll factor for the first loop */
#define P2 2    /* unroll factor for the last loop */

    _N = N&(~(P1 - 1));
    WUR_AE_CBEGIN0((uintptr_t)(state->d));
    WUR_AE_CEND0((uintptr_t)(state->d + M));
    yalign=AE_ZALIGN64();
    for (n = 0; n<_N; n += P1)
    {
      XT_LSX2IP(x0, pX, 8);
      XT_LSX2IP(x1, pX, 8);
      XT_LSX2IP(x2, pX, 8);
      XT_LSX2IP(x3, pX, 8);
      s0 = x0;
      s1 = x1;
      s2 = x2;
      s3 = x3;
      acc0 = 
      acc1 = 
      acc2 = 
      acc3 = (xtfloatx2)0.0f;
      { xtfloatx2 dummy;  XT_LSX2XC(dummy, pD, -8); }
      __Pragma("loop_count min=2,factor=2")
      for (m = 0; m<M; m += 2)
      {
        xtfloatx2 t0, t1, t2;
        xtfloatx2 tmp0;
        /* SIMD multiply */
        h0 = pH[m];
        XT_MADD_SX2(acc0, h0, x0);
        XT_MADD_SX2(acc1, h0, x1);
        XT_MADD_SX2(acc2, h0, x2);
        XT_MADD_SX2(acc3, h0, x3);

        /* select and load next sample */
        x3 = XT_SEL32_LH_SX2(x2, x3);
        t2 = XT_SEL32_LH_SX2(x1, x2);
        t1 = XT_SEL32_LH_SX2(x0, x1);
        XT_LSX2XC(tmp0, pD, -8);
        t0 = XT_SEL32_LH_SX2(tmp0, x0);

        /* SIMD multiply */
        h0 = pH[m+1];
        XT_MADD_SX2(acc0, h0, t0);
        XT_MADD_SX2(acc1, h0, t1);
        XT_MADD_SX2(acc2, h0, t2);
        XT_MADD_SX2(acc3, h0, x3);
        /* select and load next sample */
        x3 = x2;
        x2 = x1;
        x1 = x0;
        x0 = tmp0;
      }
      { xtfloatx2 dummy;  XT_LSX2XC(dummy, pD, 8); }

      XT_SASX2IP(acc0, yalign, pZ);
      XT_SASX2IP(acc1, yalign, pZ);
      XT_SASX2IP(acc2, yalign, pZ);
      XT_SASX2IP(acc3, yalign, pZ);

      NASSERT_ALIGN(pD, 8);
      XT_SSX2XC(s0, pD, 8);
      XT_SSX2XC(s1, pD, 8);
      XT_SSX2XC(s2, pD, 8);
      XT_SSX2XC(s3, pD, 8);
    }

    for (; n<N; n += P2)
    {
      XT_LSX2IP(x0, pX, 8);
      s0 = x0;
      acc0 = (xtfloatx2)0.0f;
      { xtfloatx2 dummy;  XT_LSX2XC(dummy, pD, -8); }
      __Pragma("loop_count min=2,factor=2")
      for (m = 0; m<M; m += 2)
      {
        xtfloatx2 tmp0;
        /* SIMD multiply */
        h0 = pH[m];
        XT_MADD_SX2(acc0, h0, x0);
        /* select and load next sample */
        XT_LSX2XC(tmp0, pD, -8);
        x0 = XT_SEL32_LH_SX2(tmp0, x0);
        h0 = pH[m + 1];
        XT_MADD_SX2(acc0, h0, x0);
        /* select and load next sample */
        x0 = tmp0;
      }
      { xtfloatx2 dummy;  XT_LSX2XC(dummy, pD, 8); }
      XT_SASX2IP(acc0, yalign, pZ);
      NASSERT_ALIGN(pD, 8);
      XT_SSX2XC(s0, pD, 8);
    }
    state->p = (float32_t*)pD;
    AE_SA64POS_FP(yalign,pZ);
} /* bkfirf_process() */
#else

/* process block of samples */
void bkfirf_process( bkfirf_handle_t _bkfir, 
                         float32_t * restrict  y,
                   const float32_t * restrict  x, int N )
{
  bkfirf_t *state;
  int n, m, M;
  xtfloat h0;
  xtfloat x0, x1, x2, x3;
  xtfloat s0, s1, s2, s3;
  xtfloat acc0, acc1, acc2, acc3;
  const xtfloat* restrict pX = (const xtfloat*)x;
   xtfloat* restrict pD;
  const xtfloat*   restrict pH;
  xtfloat*          pZ = (      xtfloat*)y;
  NASSERT(_bkfir);
  state=(bkfirf_t*)_bkfir;
  NASSERT(state->h);
  NASSERT(state->d);
  NASSERT(state->p);
  NASSERT_ALIGN(state->h,8);
  NASSERT_ALIGN(state->d,8);
  NASSERT_ALIGN(state->p,8);
  NASSERT(N%4==0);
  NASSERT_ALIGN(x,8);
  NASSERT((state->M%4)==0);
  NASSERT(x);
  NASSERT(y);
  if(N<=0) return;
  pD = ( xtfloat*)state->p;
  pH = (const xtfloat*  )state->h;

  M=state->M;
  NASSERT(N>0);
  NASSERT(M>0);
  WUR_AE_CBEGIN0((uintptr_t)(state->d));
  WUR_AE_CEND0((uintptr_t)(state->d + M));
  for (n = 0; n<(N>>2); n ++)
  {
    XT_LSIP(x0, castxcc(xtfloat, pX), 4);
    XT_LSIP(x1, castxcc(xtfloat, pX), 4);
    XT_LSIP(x2, castxcc(xtfloat, pX), 4);
    XT_LSIP(x3, castxcc(xtfloat, pX), 4);
    acc0 = XT_CONST_S(0);
    acc1 = XT_CONST_S(0);
    acc2 = XT_CONST_S(0);
    acc3 = XT_CONST_S(0);
    s0 = x0;
    s1 = x1;
    s2 = x2;
    s3 = x3;
    { xtfloat dummy;  XT_LSXC(dummy, pD, -4); }
    for (m = 0; m<M; m ++)
    {
      h0 = pH[m];
      XT_MADD_S(acc0, x0, h0);
      XT_MADD_S(acc1, x1, h0);
      XT_MADD_S(acc2, x2, h0);
      XT_MADD_S(acc3, x3, h0);
      x3 = x2;
      x2 = x1;
      x1 = x0;
      XT_LSXC(x0, pD, -4); 
    } 
    { xtfloat dummy;  XT_LSXC(dummy, pD, 4); }
    XT_SSXC(s0, pD, 4);
    XT_SSXC(s1, pD, 4);
    XT_SSXC(s2, pD, 4);
    XT_SSXC(s3, pD, 4);

    XT_SSIP(acc0, pZ, 4);
    XT_SSIP(acc1, pZ, 4);
    XT_SSIP(acc2, pZ, 4);
    XT_SSIP(acc3, pZ, 4);

  }
  state->p = (float32_t*)pD;
}

#endif
