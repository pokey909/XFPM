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
int64_t vec_dot24x24 (const f24 * restrict x,const f24 * restrict y,int N)
#if 0
{
  int n;
  ae_int32x2  vxw, vyw, vzw;
  ae_f24x2    vxf, vyf;
  ae_f64      vaf;
  ae_int64    vai;

  ae_valign      x_align, y_align;

  const ae_f24x2 * restrict px = (const ae_f24x2 *)x;
  const ae_f24x2 * restrict py = (const ae_f24x2 *)y;

  vaf = AE_ZERO64();
  vzw = 0;

  x_align = AE_LA64_PP(px);
  y_align = AE_LA64_PP(py);

  AE_LA32X2F24_IP(vxf, x_align, px);
  AE_LA32X2F24_IP(vyf, y_align, py);

  for (n=0; n<N-1; n+=2)
  {
    AE_MULAAFD24_HH_LL(vaf, vxf, vyf);
    AE_LA32X2F24_IP(vxf, x_align, px);
    AE_LA32X2F24_IP(vyf, y_align, py);
  }

  if (n<N)
  {
    vxw = AE_SEL32_HL(vxf, vzw);
    vyw = AE_SEL32_HL(vyf, vzw);
    vxf = AE_MOVF24X2_FROMINT32X2(vxw);
    vyf = AE_MOVF24X2_FROMINT32X2(vyw);
    AE_MULAAFD24_HH_LL(vaf, vxf, vyf);
  }
  
  vai = (vaf);
  vai = AE_SRAI64(vai, 16);

  return_int64(vai);
} /* vec_dot24x24() */
#else
{
  int n;
  ae_f24x2    vxf, vyf;
  ae_f64      vaf;
  ae_valign      x_align, y_align;

  const ae_f24x2 * restrict px = (const ae_f24x2 *)x;
  const ae_f24x2 * restrict py = (const ae_f24x2 *)y;

  if (N<=0) return 0;
  vaf = AE_ZERO64();
  x_align = AE_LA64_PP(px);
  y_align = AE_LA64_PP(py);
  for (n=0; n<(N&~1); n+=2)
  {
    AE_LA32X2F24_IP(vxf, x_align, px);
    AE_LA32X2F24_IP(vyf, y_align, py);
    AE_MULAAFD24_HH_LL(vaf, vxf, vyf);
  }
  if (N&1)
  {
    vxf=AE_L32F24_I((const ae_f24*)px,0);
    vyf=AE_L32F24_I((const ae_f24*)py,0);
    vyf=AE_SEL24_HH(vyf,0);
    AE_MULAAFD24_HH_LL(vaf, vxf, vyf);
  }

  return_int64(AE_SRAI64(vaf, 16));
} /* vec_dot24x24() */
#endif
