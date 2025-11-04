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
    Vector Sum 
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/
/* Signal Processing Library API. */
#include "NatureDSP_Signal_vector.h"
#include "common.h"
#include "common_fpu.h"

#if (HAVE_VFPU==0 && HAVE_FPU==0)
DISCARD_FUN(void,vec_addf,( float32_t * restrict z,
              const float32_t * restrict x,
              const float32_t * restrict y,
              int N))
#elif (HAVE_VFPU)
/*-------------------------------------------------------------------------
  Vector Sum
  This routine makes pair wise saturated summation of vectors.
  Two versions of routines are available: regular versions (vec_add32x32, 
  vec_add24x24, vec_add16x16, vec_addf) work with arbitrary arguments, 
  faster versions (vec_add32x32_fast, vec_add24x24_fast, vec_add16x16_fast) 
  apply some restrictions.

  Precision: 
  32x32 32-bit inputs, 32-bit output
  24x24 24-bit inputs, 24-bit output
  16x16 16-bit inputs, 16-bit output
  f     single precision floating point

  Input:
  x[N]   input data
  y[N]   input data
  N      length of vectors
  Output:
  z[N]   output data

  Restriction:
  Regular versions (vec_add32x32, vec_add24x24, vec_add16x16, vec_addf):
  x,y,z - should not be overlapped
  Faster versions (vec_add32x32_fast, vec_add24x24_fast, vec_add16x16_fast):
  z,x,y - aligned on 8-byte boundary
  N     - multiple of 4
-------------------------------------------------------------------------*/
void vec_addf ( float32_t * restrict z,
              const float32_t * restrict x,
              const float32_t * restrict y,
              int N)
{
  int n;

  xtfloatx2 vxf, vyf, vzf;
  xtfloat xf, yf, zf;

  const xtfloatx2 * restrict px = (const xtfloatx2 *)x;
  const xtfloatx2 * restrict py = (const xtfloatx2 *)y;
  xtfloatx2* restrict pz = (xtfloatx2 *)z;
  ae_valign x_align, y_align, z_align;
  NASSERT(x);
  NASSERT(y);
  NASSERT(z);
  if (N <= 0) return;

  x_align = AE_LA64_PP(px);
  y_align = AE_LA64_PP(py);
  z_align = AE_ZALIGN64();

  XT_LASX2IP(vxf, x_align, px);
  XT_LASX2IP(vyf, y_align, py);

  for (n = 0; n<N - 1; n += 2)
  {
    vzf = XT_ADD_SX2(vxf, vyf);
    XT_LASX2IP(vxf, x_align, px);
    XT_LASX2IP(vyf, y_align, py);
    XT_SASX2IP(vzf, z_align, pz);
  }
  AE_SA64POS_FP(z_align, pz);  

  if (N & 1)
  {
    xf=XT_LSI( (const xtfloat *)px, -8);
    yf=XT_LSI( (const xtfloat *)py, -8);
    zf = XT_ADD_S(xf, yf);
    XT_SSI(zf, (xtfloat *)pz, 0);
  }
} /* vec_addf() */
#else
void vec_addf ( float32_t * restrict z,
              const float32_t * restrict x,
              const float32_t * restrict y,
              int N)
{
  xtfloat x0,y0;
  int n;
        xtfloat  * restrict pZ = (      xtfloat  *)z;
  const xtfloat  * restrict pX = (const xtfloat  *)x;
  const xtfloat  * restrict pY = (const xtfloat  *)y;
  for (n=0; n<N; n++)
  {
    XT_LSIP(x0, pX, sizeof(xtfloat));
    XT_LSIP(y0, pY, sizeof(xtfloat));
    x0=XT_ADD_S(x0,y0);
    XT_SSIP(x0, pZ, sizeof(xtfloat));
  }
}
#endif
