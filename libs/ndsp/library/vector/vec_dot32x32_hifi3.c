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
  NatureDSP Signal Processing Library. Vector Operations
    Vector Dot product
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/
/* Signal Processing Library API. */
#include "NatureDSP_Signal_vector.h"
#include "common.h"

/*-------------------------------------------------------------------------
  Vector Dot product
  These routines take two vectors and calculates their dot product.
  Two versions of routines are available: regular versions (vec_dot64x32,
  vec_dot64x64, vec_dot64x64i, vec_dot24x24, vec_dot32x16, vec_dot32x32,
  vec_dot16x16, vec_dotf) work with arbitrary arguments, faster versions
  (vec_dot64x32_fast, vec_dot64x64_fast, vec_dot64x64i_fast,
  vec_dot24x24_fast, vec_dot32x16_fast, vec_dot32x32_fast,
  vec_dot16x16_fast) apply some restrictions.  
  NOTE:
  vec_dot16x16_fast utilizes 32-bit saturating accumulator, so input data 
  should be scaled properly to avoid erroneous results.

  Precision: 
  64x32  64x32-bit data, 64-bit output (fractional multiply Q63xQ31->Q63)
  64x64  64x64-bit data, 64-bit output (fractional multiply Q63xQ63->Q63)
  64x64i 64x64-bit data, 64-bit output (low 64 bit of integer multiply)
  24x24  24x24-bit data, 64-bit output
  32x32  32x32-bit data, 64-bit output
  32x16  32x16-bit data, 64-bit output
  16x16  16x16-bit data, 64-bit output for regular version and 32-bit for 
                        fast version
  f      single precision floating point

  Input:
  x[N]  input data, Q15, Q31, Q63 or floating point
  y[N]  input data, Q15, Q31, Q63 or floating point
  N	    length of vectors
  Returns:
  dot product of all data pairs, Q31, Q63 or floating point

  Restrictions:
  Regular versions:
    none
  Faster versions:
    x,y - aligned on 8-byte boundary
    N   - multiple of 4
-------------------------------------------------------------------------*/
int64_t vec_dot32x32 (const int32_t * restrict x, const int32_t * restrict y, int N)
{
#ifndef AE_MULAAFD32R_HH_LL
  int n;
  const ae_int32x2 * restrict px = (const ae_int32x2 *)x;
  const ae_int32x2 * restrict py = (const ae_int32x2 *)y;
  ae_int32x2 X, Y;
  ae_valign ax, ay;
  ae_f64 Acc;
  ae_int64 Acci;

  if (N <= 0) return 0;
  NASSERT(x);
  NASSERT(y);

  ax = AE_LA64_PP(px);
  ay = AE_LA64_PP(py);
  Acc = AE_ZERO64();
  for (n = 0; n < (N >> 1); n++)
  {
    AE_LA32X2_IP(X, ax, px);
    AE_LA32X2_IP(Y, ay, py);
    AE_MULAF32R_HH(Acc, X, Y);
    AE_MULAF32R_LL(Acc, X, Y);
  }
  if (N & 1)
  {
    X = AE_L32_I((const ae_int32*)px, 0);
    Y = AE_L32_I((const ae_int32*)py, 0);
    AE_MULAF32R_LL(Acc, X, Y);
  }

  Acci = (Acc);
  Acci = AE_SRAI64(Acci, 16);
  return_int64(Acci);
#else
  int n;
  const ae_int32x2 * restrict px = (const ae_int32x2 *)x;
  const ae_int32x2 * restrict py = (const ae_int32x2 *)y;
  ae_int32x2 X, Y;
  ae_valign ax, ay;
  ae_f64 Acc;
  ae_int64 Acci;

  if (N <= 0) return 0;
  NASSERT(x);
  NASSERT(y);

  ax = AE_LA64_PP(px);
  ay = AE_LA64_PP(py);
  Acc = AE_ZERO64();
  for (n = 0; n < (N >> 1); n++)
  {
    AE_LA32X2_IP(X, ax, px);
    AE_LA32X2_IP(Y, ay, py);
    AE_MULAAFD32R_HH_LL(Acc, X, Y);
  }
  if (N & 1)
  {
    X = AE_L32_I((const ae_int32*)px, 0);
    Y = AE_L32_I((const ae_int32*)py, 0);
    AE_MULAF32R_LL(Acc, X, Y);
  }

  Acci = (Acc);
  Acci = AE_SRAI64(Acci, 16);
  return_int64(Acci);
#endif

} /* vec_dot32x32() */
