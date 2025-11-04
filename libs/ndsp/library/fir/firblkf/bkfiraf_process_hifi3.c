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
  These functions implement FIR filter with no limitation on size of data
  block, alignment and length of impulse response at the cost of increased
  processing complexity.
  NOTE: 
  User application is not responsible for management of delay lines.

  Precision: 
  16x16    16-bit data, 16-bit coefficients, 16-bit outputs
  24x24    24-bit data, 24-bit coefficients, 24-bit outputs
  32x16    32-bit data, 16-bit coefficients, 32-bit outputs
  32x32    32-bit data, 32-bit coefficients, 32-bit outputs
  f        floating point
  Input:
  x[N]     input samples, Q15, Q31, floating point
  h[M]     filter coefficients in normal order, Q15, Q31, floating point
  N        length of sample block
  M        length of filter
  Output:
  y[N]     input samples, Q15, Q31, floating point 

  Restrictions:
  x,y      should not be overlapping
-------------------------------------------------------------------------*/
/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "bkfiraf.h"
#include "common_fpu.h"

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void,bkfiraf_process,( bkfiraf_handle_t _bkfir, 
                         float32_t * restrict  y,
                   const float32_t * restrict  x, int N ))
#elif HAVE_VFPU
/* process block of samples */
void bkfiraf_process( bkfiraf_handle_t _bkfir, 
                         float32_t * restrict  y,
                   const float32_t * restrict  x, int N )
{
    bkfiraf_t *state;
    int n, m, M, _N;
    xtfloat tmp;
    xtfloatx2 h0;
    xtfloatx2 x0, x1, x2, x3;
    xtfloat   a, s0, s1, s2, s3,s4,s5,s6,s7;
    xtfloatx2 acc0, acc1, acc2, acc3;
    ae_valign x_align, z_align;
    const xtfloatx2*  restrict pX = (const xtfloatx2*)x;
    const xtfloatx2*  restrict px = (const xtfloatx2*)x;
    const xtfloatx2* restrict pD;
    const xtfloat*   restrict pH;
    xtfloatx2*          pZ = (xtfloatx2*)y;
    NASSERT(_bkfir);
    state=(bkfiraf_t*)_bkfir;
    NASSERT(state->h);
    NASSERT(state->d);
    NASSERT(state->p);
    NASSERT_ALIGN(state->h,8);
    NASSERT_ALIGN(state->d,8);
    NASSERT((state->M%4)==0);
    NASSERT(x);
    NASSERT(y);
    if(N<=0) return;
    pD = (const xtfloatx2*)state->p;
    pH = (const xtfloat*)state->h;
    M=state->M;
    NASSERT(N>0);
    NASSERT(M>0);
#define P 8    /* unroll factor for the first loop */

    _N = N&(~(P - 1));
    WUR_AE_CBEGIN0((uintptr_t)(state->d));
    WUR_AE_CEND0((uintptr_t)(state->d + M));
    x_align = AE_LA64_PP(pX);
    z_align = AE_ZALIGN64();
    for (n = 0; n<_N; n += P)
    {
        XT_LASX2IP(x0, x_align, pX);
        XT_LASX2IP(x1, x_align, pX);
        XT_LASX2IP(x2, x_align, pX);
        XT_LASX2IP(x3, x_align, pX);

        acc0 = 
        acc1 = 
        acc2 = 
        acc3 = (xtfloatx2)0.0f;
        XT_LSXC(tmp, castxcc(const xtfloat,pD), -4);
        __Pragma("loop_count min=1,factor=4")
        for (m = 0; m<M; m++)
        {
            /* SIMD multiply */
            h0 = pH[m];
            XT_MADD_SX2(acc0, h0, x0);
            XT_MADD_SX2(acc1, h0, x1);
            XT_MADD_SX2(acc2, h0, x2);
            XT_MADD_SX2(acc3, h0, x3);

            /* select and load next sample */
            x3 = XT_SEL32_LH_SX2(x2, x3);
            x2 = XT_SEL32_LH_SX2(x1, x2);
            x1 = XT_SEL32_LH_SX2(x0, x1);
            XT_LSXC(tmp, castxcc(const xtfloat,pD), -4);
            x0 = XT_SEL32_LH_SX2((xtfloatx2)tmp, x0);
        }
        XT_LSXC(tmp, castxcc(const xtfloat,pD), 4);
        XT_SASX2IP(acc0, z_align, pZ);
        XT_SASX2IP(acc1, z_align, pZ);
        XT_SASX2IP(acc2, z_align, pZ);
        XT_SASX2IP(acc3, z_align, pZ);

        XT_LSIP(s0, castxcc(const xtfloat,px), 4);
        XT_LSIP(s1, castxcc(const xtfloat,px), 4);
        XT_LSIP(s2, castxcc(const xtfloat,px), 4);
        XT_LSIP(s3, castxcc(const xtfloat,px), 4);
        XT_LSIP(s4, castxcc(const xtfloat,px), 4);
        XT_LSIP(s5, castxcc(const xtfloat,px), 4);
        XT_LSIP(s6, castxcc(const xtfloat,px), 4);
        XT_LSIP(s7, castxcc(const xtfloat,px), 4);

        XT_SSXC(s0, castxcc(xtfloat,pD), 4);
        XT_SSXC(s1, castxcc(xtfloat,pD), 4);
        XT_SSXC(s2, castxcc(xtfloat,pD), 4);
        XT_SSXC(s3, castxcc(xtfloat,pD), 4);
        XT_SSXC(s4, castxcc(xtfloat,pD), 4);
        XT_SSXC(s5, castxcc(xtfloat,pD), 4);
        XT_SSXC(s6, castxcc(xtfloat,pD), 4);
        XT_SSXC(s7, castxcc(xtfloat,pD), 4);
    }
    AE_SA64POS_FP(z_align, pZ);
    for (; n<N; n ++)
    {
        XT_LSIP(s0, castxcc(const xtfloat,pX), 4);
        a = (xtfloatx2)0.0f;
        XT_LSXC(tmp, castxcc(const xtfloat,pD), -4);
        __Pragma("loop_count min=1,factor=4")
        for (m = 0; m<M; m++)
        {
            s1 = pH[m];
            XT_MADD_S(a, s1, s0);
            XT_LSXC(s0, castxcc(const xtfloat,pD), -4);
        }
        XT_LSXC(tmp, castxcc(const xtfloat,pD), 4);
        XT_SSIP(a,   castxcc(      xtfloat,pZ), 4);

        XT_LSIP(s0, castxcc(const xtfloat,px), 4);
        XT_SSXC(s0, castxcc(      xtfloat,pD), 4);
    }
    state->p = (float32_t*)pD;
} /* bkfiraf_process() */
#else   
// for SFPU
void bkfiraf_process( bkfiraf_handle_t _bkfir, 
                         float32_t * restrict  y,
                   const float32_t * restrict  x, int N )
{
    bkfiraf_t *state;
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
    state = (bkfiraf_t*)_bkfir;
    NASSERT(state->h);
    NASSERT(state->d);
    NASSERT(state->p);
    NASSERT_ALIGN(state->h,8);
    NASSERT_ALIGN(state->d,8);
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
    for (n = (4*n); n<N; n++)
    {
        XT_LSIP(x0, castxcc(xtfloat, pX), 4);
        acc0 = XT_CONST_S(0);
        s0 = x0;
        { xtfloat dummy;  XT_LSXC(dummy, pD, -4); }
        for (m = 0; m<M; m ++)
        {
            h0 = pH[m];
            XT_MADD_S(acc0, x0, h0);
            XT_LSXC(x0, pD, -4); 
        } 
        { xtfloat dummy;  XT_LSXC(dummy, pD, 4); }
        XT_SSXC(s0, pD, 4);
        XT_SSIP(acc0, pZ, 4);
    }
    state->p = (float32_t*)pD;
}
#endif
