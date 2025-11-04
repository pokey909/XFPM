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
    Complex block FIR filter, floating point
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Complex Block FIR Filter
  Computes a complex FIR filter (direct-form) using complex IR stored in 
  vector h. The complex data input is stored in vector x. The filter output
  result is stored in vector r. The filter calculates N output samples using 
  M coefficients and requires last M-1 samples in the delay line. Real and
  imaginary parts are interleaved and real parts go first (at even indexes).
  NOTE: 
  1. User application is not responsible for management of delay lines
  2. User has an option to set IR externally or copy from original location 
     (i.e. from the slower constant memory). In the first case, user is 
     responsible for right alignment, ordering and zero padding of filter 
     coefficients - usually array is composed from zeroes (left padding), 
     reverted IR and right zero padding.

  Precision: 
  16x16     16-bit data, 16-bit coefficients, 16-bit outputs
  24x24     24-bit data, 24-bit coefficients, 24-bit outputs
  32x16     32-bit data, 16-bit coefficients, 32-bit outputs
  32x32     32-bit data, 32-bit coefficients, 32-bit outputs
  f         floating point

  Input:
  h[M]      complex filter coefficients; h[0] is to be multiplied with the 
            newest sample, Q15, Q31, floating point
  x[N]      input samples, Q15, Q31, floating point
  N         length of sample block (in complex samples) 
  M         length of filter 
  extIR     if zero, IR is copied from original location, otherwise not
            but user should keep alignment, order of coefficients 
            and zero padding requirements shown below
  Output:			
  y[N]      output samples, Q15, Q31, floating point

  Alignment, ordering and zero padding for external IR  (extIR!=0)
  -----------------+----------+--------------+--------------+----------------
  Function	       |Alignment,|Left zero     |   Coefficient| Right zero 
                   | bytes    |padding, bytes|   order      | padding, bytes
  -----------------+----------+--------------+--------------+----------------
  cxfir16x16_init, |     8    |    4         |  inverted    |  4
  hifi3            |          |              |              |
  cxfir16x16_init, |     8    |  2 before    |  *           |  6 after
  hifi3z           |          |  each copy   |              |  each copy
  cxfir24x24_init  |     8    |    0         |  inverted    |  0
  cxfir32x16_init  |     8    |    4         |  inverted    |  4
  cxfir32x32_init  |     8    |    0         |  inv,conj    |  0
  cxfirf_init      |     8    |    0         |  direct      |  0
  -----------------+----------+--------------+--------------+----------------
  * inverted: conjugated copy and (imaginary; real) copy at 4*(M+4) bytes offset

  Restriction:
  x,y       should not overlap
  x,h       aligned on a 8-bytes boundary
  N,M       multiples of 4
-------------------------------------------------------------------------*/
/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "cxfirf.h"
#include "common_fpu.h"

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void,cxfirf_process,( cxfirf_handle_t _cxfir, 
                         complex_float * restrict  y,
                   const complex_float * restrict  x, int N ))
#elif HAVE_VFPU

/* process block of samples */
void cxfirf_process( cxfirf_handle_t _cxfir, 
                         complex_float * restrict  y,
                   const complex_float * restrict  x, int N )
{
    int n, m, M, _N;
    xtfloatx2 x0, x1, x2, x3;
    xtfloatx2 s0, s1, s2, s3;
    xtfloatx2 acc0, acc1, acc2, acc3;
    const xtfloatx2* restrict pX = (const xtfloatx2*)x;
          xtfloatx2* restrict pD;
    const xtfloatx2* restrict pH;
    xtfloatx2*          pZ = (xtfloatx2*)y;
    cxfirf_t* state;
    NASSERT(_cxfir);
    state=(cxfirf_t*)_cxfir;
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
    M=state->M;
    pD = (      xtfloatx2*)state->p;
    pH = (const xtfloatx2*)state->h;
    NASSERT(N>0);
    NASSERT(M>0);
#define P1 4    /* unroll factor for the first loop */

    _N = N&(~(P1 - 1));
    WUR_AE_CBEGIN0((uintptr_t)(state->d));
    WUR_AE_CEND0((uintptr_t)(state->d + M));
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
        pH = (const xtfloatx2 *)state->h;
        { xtfloatx2 dummy;  XT_LSX2XC(dummy, pD, -8); }
        __Pragma("loop_count min=1,factor=4")
        for (m = 0; m<M; m ++)
        {
            xtfloatx2 hm;
            XT_LSX2IP(hm, pH, 8);
            /*Complex SIMD multiply */
            XT_MADDC_S(acc0, hm, x0);
            XT_MADDC_S(acc1, hm, x1);
            XT_MADDC_S(acc2, hm, x2);
            XT_MADDC_S(acc3, hm, x3);
            /*load next sample */
            x3 = x2; x2 = x1; x1 = x0;
            XT_LSX2XC(x0, pD, -8);
        }
        { xtfloatx2 dummy;  XT_LSX2XC(dummy, pD, 8); }

        XT_SSX2IP(acc0, pZ, 8);
        XT_SSX2IP(acc1, pZ, 8);
        XT_SSX2IP(acc2, pZ, 8);
        XT_SSX2IP(acc3, pZ, 8);

        NASSERT_ALIGN(pD, 8);
        XT_SSX2XC(s0, pD, 8);
        XT_SSX2XC(s1, pD, 8);
        XT_SSX2XC(s2, pD, 8);
        XT_SSX2XC(s3, pD, 8);
    }
    __Pragma("loop_count max=3")
    for (; n<N; n++)
    {
        XT_LSX2IP(x0, pX, 8);
        s0 = x0;
        acc0 = 
        acc1 = (xtfloatx2)0.0f;
        pH = (const xtfloatx2 *)state->h;
        { xtfloatx2 dummy;  XT_LSX2XC(dummy, pD, -8); }
        __Pragma("loop_count min=1,factor=4")
        for (m = 0; m<M; m++)
        {
            xtfloatx2 hm;
            XT_LSX2IP(hm, pH, 8);
            XT_MADDMUX_S( acc0, hm, x0, 0);
            XT_MADDMUX_S( acc1, hm, x0, 1);
            XT_LSX2XC(x0, pD, -8);
        }
        { xtfloatx2 dummy;  XT_LSX2XC(dummy, pD, 8); }
        acc0 = XT_ADD_SX2(acc0, acc1);
        XT_SSX2IP(acc0, pZ, 8);
        NASSERT_ALIGN(pD, 8);
        XT_SSX2XC(s0, pD, 8);
    }
    state->p = (complex_float*)pD;
} /* cxfirf_process() */
#else

/* process block of samples */
void cxfirf_process( cxfirf_handle_t _cxfir, 
                         complex_float * restrict  y,
                   const complex_float * restrict  x, int N )
{
  int n, m, M;
  xtfloat x0, x1, x2, x3;
  xtfloat s0, s1, s2, s3;
  xtfloat acc0, acc1, acc2, acc3;
  const xtfloat* restrict pX = (const xtfloat*)x;
        xtfloat* restrict pD;
  const xtfloat* restrict pH;
        xtfloat*          pZ = (xtfloat*)y;
  cxfirf_t* state;
  NASSERT(_cxfir);
  state=(cxfirf_t*)_cxfir;
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
  M=state->M;
  pD = ( xtfloat*)state->p;
  pH = (const xtfloat*)state->h;
  NASSERT(N>0);
  NASSERT(M>0);
  WUR_AE_CBEGIN0((uintptr_t)(state->d));
  WUR_AE_CEND0((uintptr_t)(state->d + M));
  for (n = 0; n<(N>>1); n ++)
  {
    XT_LSIP(x0, pX, 4);
    XT_LSIP(x1, pX, 4);
    XT_LSIP(x2, pX, 4);
    XT_LSIP(x3, pX, 4);
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
      xtfloat h0, h1;
      h0 = pH[2*m + 0];
      h1 = pH[2*m + 1];
      XT_MADD_S(acc0, x0, h0); /*re*re*/
      XT_MADD_S(acc1, x0, h1); /*re*im*/

      XT_MADD_S(acc2, x2, h0);
      XT_MADD_S(acc3, x2, h1);

      XT_MSUB_S(acc0, x1, h1); /*-im*im*/
      XT_MADD_S(acc1, x1, h0); /*im*re*/

      XT_MSUB_S(acc2, x3, h1);
      XT_MADD_S(acc3, x3, h0);
      x3 = x1;
      x2 = x0;
      XT_LSXC(x1, pD, -4); 
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
  state->p = (complex_float*)pD;
} // cxfirf_process()
#endif
