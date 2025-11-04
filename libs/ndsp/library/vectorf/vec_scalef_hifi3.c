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
    Vector Scaling with Saturation
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/
/* Signal Processing Library API. */
#include "NatureDSP_Signal_vector.h"
#include "common.h"
#include "common_fpu.h"

#if (HAVE_VFPU==0 && HAVE_FPU==0)
DISCARD_FUN(void,vec_scalef,(     float32_t * restrict y,
                    const float32_t * restrict x,
                    float32_t s,
                    int N))
#elif (HAVE_VFPU)
/*-------------------------------------------------------------------------
  Vector Scaling with Saturation
  These routines make shift with saturation of data values in the vector 
  by given scale factor (degree of 2).
  Functions vec_scale() make multiplication of a vector to a coefficient 
  which is not a power of 2 forming Q31, Q15 or floating-point result.
  Two versions of routines are available: regular versions (vec_shift24x24, 
  vec_shift32x32, vec_shift16x16, vec_shiftf, vec_scale32x24, vec_scale32x32, 
  vec_scale24x24, vec_scale16x16, vec_scalef, vec_scale_sf) work with 
  arbitrary arguments, faster versions (vec_shift24x24_fast, vec_shift32x32_fast,
  vec_shift16x16_fast, vec_scale32x24_fast, vec_scale32x32_fast, 
  vec_scale24x24_fast, vec_scale16x16_fast) apply some restrictions

  For floating point:
  Fuction vec_shiftf() makes scaling without saturation of data values by given 
  scale factor (degree of 2). 
  Functions vec_scalef() and vec_scale_sf() make multiplication of input vector
  to coefficient which is not a power of 2.
  Two versions of routines are available: 
    without saturation - vec_scalef;
    with saturation - vec_scale_sf; 

  Precision:
  24x24	24-bit input, 24-bit output
  32x24	32-bit input, 32-bit output, 24-bit scale factor
  32x32 32-bit input, 32-bit output
  16x16 16-bit input, 16-bit output
  f     single precision floating point

  Input:
  x[N]    input data, Q31, Q15 or floating point
  t       shift count. If positive, it shifts left with saturation, if
          negative it shifts right
  s       scale factor, Q31, Q15 or floating point
  N       length of vector
  fmin    minimum output value (only for vec_scale_sf)
  fmax    maximum output value (only for vec_scale_sf)
  Output:
  y[N]    output data, Q31, Q15 or floating point

  Restrictions:
  For regular versions (vec_shift24x24, vec_shift32x32, vec_shift16x16, 
  vec_shiftf, vec_scale32x24, vec_scale32x32, vec_scale24x24, vec_scale16x16, 
  vec_scalef, vec_scale_sf):
  x,y should not overlap
  t   should be in range -31...31 for fixed-point functions and -129...146 
      for floating point
  For vec_scale_sf:
  fmin<=fmax;

  For faster versions (vec_shift24x24_fast, vec_shift32x32_fast, 
  vec_shift16x16_fast, vec_scale32x24_fast, vec_scale32x32_fast,
  vec_scale24x24_fast, vec_scale16x16_fast):
  x,y should not overlap
  t should be in range -31...31 
  x,y - aligned on 8-byte boundary
  N   - multiple of 4 
-------------------------------------------------------------------------*/
void vec_scalef     (     float32_t * restrict y,
                    const float32_t * restrict x,
                    float32_t s,
                    int N)
{
  int n;

  xtfloatx2 vxf, vyf, vsf;
  xtfloat xf, yf;

  const xtfloatx2* restrict px = (const xtfloatx2 *)x;
        xtfloatx2* restrict py = (      xtfloatx2 *)y;
  ae_valign x_align, y_align;
  NASSERT(x);
  NASSERT(y);
  if (N <= 0) return;

  x_align = AE_LA64_PP(px);
  y_align = AE_ZALIGN64();

  vsf = XT_MOV_SX2(s);
  for (n = 0; n<N-1; n += 2)
  {
    XT_LASX2IP(vxf, x_align, px);
    vyf = XT_MUL_SX2(vxf, vsf);
    XT_SASX2IP(vyf, y_align, py);
  }
  AE_SA64POS_FP(y_align, py);

  if (N & 1)
  {
    xf = XT_LSI((const xtfloat *)px, 0);
    yf = XT_MUL_S(xf, s);
    XT_SSI(yf, (xtfloat *)py, 0);
  }
} /* vec_scalef() */
#else
void vec_scalef     (     float32_t * restrict y,
                    const float32_t * restrict x,
                    float32_t s,
                    int N)
{
          xtfloat  * restrict pY = (      xtfloat  *)y;
    const xtfloat  * restrict pX = (const xtfloat  *)x;
    int n;
    for (n=0; n<N; n++)
    {
        xtfloat x0;
        XT_LSIP(x0, pX, sizeof(xtfloat));
        x0=XT_MUL_S(x0,s);
        XT_SSIP(x0, pY, sizeof(xtfloat));
    }
}
#endif
