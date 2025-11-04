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
#include "common_fpu.h"

#if (HAVE_VFPU==0 && HAVE_FPU==0)
DISCARD_FUN(float32_t,vec_dotf, (const float32_t * restrict x,const float32_t * restrict y,int N))
#elif (HAVE_VFPU)
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
float32_t vec_dotf   (const float32_t * restrict x,const float32_t * restrict y,int N)

{
  int n;

  xtfloatx2 vxf, vyf, vacc;
  xtfloat xf, yf, zf, acc;

  const xtfloatx2 * restrict px = (const xtfloatx2 *)x;
  const xtfloatx2 * restrict py = (const xtfloatx2 *)y;
  ae_valign x_align, y_align;
  NASSERT(x);
  NASSERT(y);
  if (N <= 0) return 0;

  x_align = AE_LA64_PP(px);
  y_align = AE_LA64_PP(py);

  XT_LASX2IP(vxf, x_align, px);
  XT_LASX2IP(vyf, y_align, py);
  vacc = XT_MOV_SX2(0.f);
  zf = XT_MOV_S(0.f);
  for (n = 0; n<N - 1; n += 2)
  {
    XT_MADD_SX2(vacc, vxf, vyf);
    XT_LASX2IP(vxf, x_align, px);
    XT_LASX2IP(vyf, y_align, py);
  }
  if (N & 1)
  {
    xf = XT_LSI((const xtfloat *)px, -8);
    yf = XT_LSI((const xtfloat *)py, -8);
    zf = XT_MUL_S(xf, yf);
  }
  acc=XT_RADD_SX2(vacc);
  acc = XT_ADD_S(acc,zf);

  return acc;
} /* vec_dotf() */
#elif (HAVE_FPU)
float32_t vec_dotf   (const float32_t * restrict x,const float32_t * restrict y,int N)
{
  xtfloat acc0, acc1,x0,y0;
  int n;
  const xtfloat  * restrict pX = (const xtfloat  *)x;
  const xtfloat  * restrict pY = (const xtfloat  *)y;
  if (N <= 0) return 0.f;
  acc0 = acc1 = XT_CONST_S(0);
  for (n = 0; n<(N&~1); n+=2)
  {
    XT_LSIP(x0, pX, sizeof(xtfloat));
    XT_LSIP(y0, pY, sizeof(xtfloat));
    XT_MADD_S(acc0,x0,y0);
    XT_LSIP(x0, pX, sizeof(xtfloat));
    XT_LSIP(y0, pY, sizeof(xtfloat));
    XT_MADD_S(acc1,x0,y0);
  }
  if (N&1)
  {
    XT_LSIP(x0, pX, sizeof(xtfloat));
    XT_LSIP(y0, pY, sizeof(xtfloat));
    XT_MADD_S(acc0,x0,y0);
  }
  return XT_ADD_S(acc0 , acc1 );
}
#endif
