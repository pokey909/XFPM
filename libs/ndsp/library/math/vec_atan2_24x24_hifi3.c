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
  NatureDSP Signal Processing Library. Math functions
    Full-Quadrant Arc Tangent
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/


/*-------------------------------------------------------------------------
  Full-Quadrant Arc Tangent
  The functions compute the arc tangent of the ratios y[N]/x[N] and store the
  result to output vector z[N]. 
  Floating point functions output is in radians. Fixed point functions
  scale its output by pi.

  NOTE:
  1. Scalar floating point function is compatible with standard ANSI C routines and set 
     errno and exception flags accordingly
  2. Scalar floating point function assigns EDOM to errno whenever y==0 and x==0.

  Accuracy:
  24 bit version: 768 (3.57e-7)
  floating point: 2 ULP

  Special cases:
       y    |   x   |  result   |  extra conditions    
    --------|-------|-----------|---------------------
     +/-0   | -0    | +/-pi     |
     +/-0   | +0    | +/-0      |
     +/-0   |  x    | +/-pi     | x<0
     +/-0   |  x    | +/-0      | x>0
     y      | +/-0  | -pi/2     | y<0
     y      | +/-0  |  pi/2     | y>0
     +/-y   | -inf  | +/-pi     | finite y>0
     +/-y   | +inf  | +/-0      | finite y>0
     +/-inf | x     | +/-pi/2   | finite x
     +/-inf | -inf  | +/-3*pi/4 | 
     +/-inf | +inf  | +/-pi/4   |

  Input:
    y[N]  vector of numerator values, Q31 or floating point
    x[N]  vector of denominator values, Q31 or floating point
    N     length of vectors
  Output:
    z[N]  results, Q31 or floating point

---------------------------------------------------------------------------*/

#include "common.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_math.h"

#include "scl_atan2_24x24_table.h"

void vec_atan2_24x24(f24 * z, const f24 * y, const f24 * x, int N)
{
  const ae_f24x2 * restrict px = (const ae_f24x2 *)x;
  const ae_f24x2 * restrict py = (const ae_f24x2 *)y;
        ae_f24x2 * restrict pz = (      ae_f24x2 *)z;
  const ae_f24   *  pt = (const ae_f24   *)polyatan24x24q23;
  int n, ey;
  ae_valign x_align, y_align, z_align;
  ae_f24x2 tmp, xf, yf, zf, ef, x2f, onef;
  ae_int32x2 x0, y0, t0, zero;
  ae_int32x2 yl, yh, xl, xh;
  ae_f32x2 z0;
  xtbool2 sx, sy, small;
  if (N <= 0) return; 
  x_align = AE_LA64_PP(px);
  y_align = AE_LA64_PP(py);
  z_align = AE_ZALIGN64();
  zero = AE_MOVI(0);
  onef = 0x00400000;
  for (n = 0; n<(N>>1) ; n++)
  {
    AE_LA32X2F24_IP(tmp, x_align, px); x0 = tmp;
    AE_LA32X2F24_IP(tmp, y_align, py); y0 = tmp;
    sx = AE_LT32(x0, zero);
    sy = AE_LT32(y0, zero);
    t0 = AE_NEG32S(x0);
    AE_MOVT32X2(x0, t0, sx);
    t0 = AE_NEG32S(y0);
    AE_MOVT32X2(y0, t0, sy);
    small = AE_LT32(x0, y0);
    t0 = x0; x0 = AE_MIN32(x0, y0); y0 = AE_MAX32(t0, y0);

    ey = AE_NSAZ32_L(y0) - 8;
    yl = AE_SLAA32S(y0, ey);
    xl = AE_SLAA32S(x0, ey);
    y0 = AE_SEL32_LH(y0, y0);
    x0 = AE_SEL32_LH(x0, x0);
    ey = AE_NSAZ32_L(y0) - 8;
    yh = AE_SLAA32S(y0, ey);
    xh = AE_SLAA32S(x0, ey);
    y0 = AE_SEL32_LL(yh, yl);
    x0 = AE_SEL32_LL(xh, xl);

    yf = AE_MOVF24X2_FROMINT32X2(y0);
    xf = AE_MOVF24X2_FROMINT32X2(x0);
   

    /* divide X/Y in Q23 via reciprocal */
    zf = AE_SUB24S(0x00BB5000, yf);

    z0 = AE_MULFP24X2RA(zf, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    ef = AE_SUB24S(onef, ef);
    ef = AE_ADD24S(ef, ef);
    z0 = zf;
    AE_MULAFP24X2RA(z0, zf, ef);
    zf = AE_MOVF24X2_FROMF32X2(z0);

    z0 = AE_MULFP24X2RA(zf, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    ef = AE_SUB24S(onef, ef);
    ef = AE_ADD24S(ef, ef);
    z0 = zf;
    AE_MULAFP24X2RA(z0, zf, ef);
    zf = AE_MOVF24X2_FROMF32X2(z0);

    z0 = AE_MULFP24X2RA(zf, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    ef = AE_SUB24S(onef, ef);
    ef = AE_ADD24S(ef, ef);
    z0 = zf;
    AE_MULAFP24X2RA(z0, zf, ef);
    zf = AE_MOVF24X2_FROMF32X2(z0);

    z0 = AE_MULFP24X2RA(zf, xf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    ef = AE_ADD24S(ef, ef);

    /* compute atan via polynomial */
    z0 = AE_MULFP24X2RA(ef, ef);
    xf = ef;
    x2f = AE_MOVF24X2_FROMF32X2(z0);

    AE_L32F24_IP(yf,pt,+4);
    z0 = AE_MULFP24X2RA(x2f, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    AE_L32F24_IP(yf, pt, +4);
    yf = AE_ADD24S(yf, ef);

    z0 = AE_MULFP24X2RA(x2f, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    AE_L32F24_IP(yf, pt, +4);
    yf = AE_ADD24S(yf, ef);

    z0 = AE_MULFP24X2RA(x2f, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    AE_L32F24_IP(yf, pt, +4);
    yf = AE_ADD24S(yf, ef);

    z0 = AE_MULFP24X2RA(x2f, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    AE_L32F24_IP(yf, pt, +4);
    yf = AE_ADD24S(yf, ef);

    z0 = AE_MULFP24X2RA(x2f, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    AE_L32F24_IP(yf, pt, +4);
    yf = AE_ADD24S(yf, ef);

    z0 = AE_MULFP24X2RA(x2f, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    AE_L32F24_IP(yf, pt, +4);
    yf = AE_ADD24S(yf, ef);

    z0 = AE_MULFP24X2RA(x2f, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    AE_L32F24_IP(yf, pt, -7*4);
    yf = AE_ADD24S(yf, ef);

    z0 = AE_MULFP24X2RA(xf, yf);
    /* move to the right octant */
    y0 = AE_SUB32(0x400000, z0);
    t0 = z0;
    AE_MOVT32X2(t0, y0, small);
    y0 = AE_SUB32(0x800000, t0);
    AE_MOVT32X2(t0, y0, sx);
    y0 = AE_NEG32(t0);
    AE_MOVT32X2(t0, y0, sy);

    tmp = AE_SAT24S(t0); AE_SA32X2F24_IP(tmp, z_align, pz);
  }
  AE_SA64POS_FP(z_align, pz);
  if (N & 1)
  {
    AE_L32F24_IP(tmp, castxcc(ae_f24, px), 0); x0 = tmp;
    AE_L32F24_IP(tmp, castxcc(ae_f24, py), 0); y0 = tmp;
    sx = AE_LT32(x0, zero);
    sy = AE_LT32(y0, zero);
    t0 = AE_NEG32S(x0);
    AE_MOVT32X2(x0, t0, sx);
    t0 = AE_NEG32S(y0);
    AE_MOVT32X2(y0, t0, sy);
    small = AE_LT32(x0, y0);
    t0 = x0; x0 = AE_MIN32(x0, y0); y0 = AE_MAX32(t0, y0);

    ey = AE_NSAZ32_L(y0) - 8;
    y0 = AE_SLAA32S(y0, ey);
    x0 = AE_SLAA32S(x0, ey);

    yf = AE_MOVF24X2_FROMINT32X2(y0);
    xf = AE_MOVF24X2_FROMINT32X2(x0);


    /* divide X/Y in Q23 via reciprocal */
    zf = AE_SUB24S(0x00BB5000, yf);

    z0 = AE_MULFP24X2RA(zf, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    ef = AE_SUB24S(0x00400000, ef);
    ef = AE_ADD24S(ef, ef);
    z0 = zf;
    AE_MULAFP24X2RA(z0, zf, ef);
    zf = AE_MOVF24X2_FROMF32X2(z0);

    z0 = AE_MULFP24X2RA(zf, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    ef = AE_SUB24S(0x00400000, ef);
    ef = AE_ADD24S(ef, ef);
    z0 = zf;
    AE_MULAFP24X2RA(z0, zf, ef);
    zf = AE_MOVF24X2_FROMF32X2(z0);

    z0 = AE_MULFP24X2RA(zf, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    ef = AE_SUB24S(0x00400000, ef);
    ef = AE_ADD24S(ef, ef);
    z0 = zf;
    AE_MULAFP24X2RA(z0, zf, ef);
    zf = AE_MOVF24X2_FROMF32X2(z0);

    z0 = AE_MULFP24X2RA(zf, xf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    ef = AE_ADD24S(ef, ef);

    /* compute atan via polynomial */
    z0 = AE_MULFP24X2RA(ef, ef);
    xf = ef;
    x2f = AE_MOVF24X2_FROMF32X2(z0);

    AE_L32F24_IP(yf, pt, +4);
    // yf = AE_MOVF24X2_FROMF32X2(-11369);
    z0 = AE_MULFP24X2RA(x2f, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    AE_L32F24_IP(yf, pt, +4);
    yf = AE_ADD24S(yf, ef);

    z0 = AE_MULFP24X2RA(x2f, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    AE_L32F24_IP(yf, pt, +4);
    yf = AE_ADD24S(yf, ef);

    z0 = AE_MULFP24X2RA(x2f, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    AE_L32F24_IP(yf, pt, +4);
    yf = AE_ADD24S(yf, ef);

    z0 = AE_MULFP24X2RA(x2f, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    AE_L32F24_IP(yf, pt, +4);
    yf = AE_ADD24S(yf, ef);

    z0 = AE_MULFP24X2RA(x2f, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    AE_L32F24_IP(yf, pt, +4);
    yf = AE_ADD24S(yf, ef);

    z0 = AE_MULFP24X2RA(x2f, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    AE_L32F24_IP(yf, pt, +4);
    yf = AE_ADD24S(yf, ef);

    z0 = AE_MULFP24X2RA(x2f, yf);
    ef = AE_MOVF24X2_FROMF32X2(z0);
    AE_L32F24_IP(yf, pt, -7 * 4);
    yf = AE_ADD24S(yf, ef);

    z0 = AE_MULFP24X2RA(xf, yf);
    /* move to the right octant */
    y0 = AE_SUB32(0x400000, z0);
    t0 = z0;
    AE_MOVT32X2(t0, y0, small);
    y0 = AE_SUB32(0x800000, t0);
    AE_MOVT32X2(t0, y0, sx);
    y0 = AE_NEG32(t0);
    AE_MOVT32X2(t0, y0, sy);

    tmp = AE_SAT24S(t0); AE_S32F24_L_IP(tmp, castxcc(ae_f24, pz), +4);
  }
  return;
} /* vec_atan2_24x24() */
