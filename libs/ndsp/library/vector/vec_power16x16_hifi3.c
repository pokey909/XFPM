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
    Power of a Vector
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/
/* Signal Processing Library API. */
#include "NatureDSP_Signal_vector.h"
#include "NatureDSP_types.h"
#include "common.h"

/*-------------------------------------------------------------------------
  Power of a Vector
  These routines compute power of vector with scaling output result by rsh 
  bits. Fixed point rountines make accumulation in the 64-bit wide 
  accumulator and output may scaled down with saturation by rsh bits. 
  So, if representation of x input is Qx, result will be represented in 
  Q(2x-rsh) format.
  Two versions of routines are available: regular versions (vec_power24x24, 
  vec_power32x32, vec_power16x16, vec_powerf) work with arbitrary arguments, 
  faster versions (vec_power24x24_fast, vec_power32x32_fast, 
  vec_power16x16_fast) apply some restrictions.

  Precision: 
  24x24 24x24-bit data, 64-bit output
  32x32 32x32-bit data, 64-bit output
  16x16 16x16-bit data, 64-bit output
  f     single precision floating point

  Input:
  x[N]  input data, Q31, Q15 or floating point
  rsh   right shift of result
  N     length of vector
  Returns: 
  Sum of squares of a vector, Q(2x-rsh)

  Restrictions:
  for vec_power32x32(): rsh in range 31...62
  for vec_power24x24(): rsh in range 15...46
  for vec_power16x16(): rsh in range 0...31
  For regular versions (vec_power24x24, vec_power32x32, 
  vec_power16x16, vec_powerf):
  none
  For faster versions (vec_power24x24_fast, 
  vec_power32x32_fast, vec_power16x16_fast ):
  x - aligned on 8-byte boundary
  N - multiple of 4
-------------------------------------------------------------------------*/
int64_t vec_power16x16 (const int16_t * restrict x, int rsh, int N)
{
#ifndef AE_MULAAAAQ16
  int n;

  ae_f64      vaf;
  ae_int64    vai, vzi;
  ae_int16x4  vxh;
  ae_f24x2    vxf;
  ae_int32x2  vxw, vzw;
  ae_valign      x_align;
  const ae_int16x4 * restrict px = (const ae_int16x4 *)x;

  NASSERT(x);
  NASSERT(rsh>=0 && rsh<=31);
  if (N <= 0) return 0;

  vzi = AE_ZERO64();
  vaf = (vzi);
  vzw = AE_MOVINT32X2_FROMINT64(vzi);

  x_align = AE_LA64_PP(px);

  AE_LA16X4_IP(vxh, x_align, px);

  for (n=0; n<N-3; n+=4)
  {
    vxw = AE_SEXT32X2D16_32(vxh);
    vxf = AE_MOVF24X2_FROMINT32X2(vxw);
    AE_MULAAFD24_HH_LL(vaf, vxf, vxf);
    vxw = AE_SEXT32X2D16_10(vxh);
    vxf = AE_MOVF24X2_FROMINT32X2(vxw);
    AE_MULAAFD24_HH_LL(vaf, vxf, vxf);
    AE_LA16X4_IP(vxh, x_align, px);
  }

  switch(N&3)
  {
  case 1:
    vxw = AE_SEXT32X2D16_32(vxh);
    vxf = AE_SEL24_HL(vxw, vzw);
    AE_MULAAFD24_HH_LL(vaf, vxf, vxf);
    break;
  case 2:
    vxw = AE_SEXT32X2D16_32(vxh);
    vxf = AE_MOVF24X2_FROMINT32X2(vxw);
    AE_MULAAFD24_HH_LL(vaf, vxf, vxf);
    break;
  case 3:
    vxw = AE_SEXT32X2D16_32(vxh);
    vxf = AE_MOVF24X2_FROMINT32X2(vxw);
    AE_MULAAFD24_HH_LL(vaf, vxf, vxf);
    vxw = AE_SEXT32X2D16_10(vxh);
    vxf = AE_SEL24_HL(vxw, vzw);
    AE_MULAAFD24_HH_LL(vaf, vxf, vxf);
    break;
  default:
    break;
  }

  vai = (vaf);
  vai = AE_SRAA64(vai, rsh+1);

  return_int64(vai);
#else
    int n;
    ae_int64    vai;
    ae_int16x4  vxh;
    ae_valign   x_align;
    const ae_int16x4 * restrict px = (const ae_int16x4 *)x;
    xtbool4 b;

    NASSERT(x);
    NASSERT(rsh >= 0 && rsh <= 31);
    if (N <= 0) return 0;

    vai = 0;

    x_align = AE_LA64_PP(px);

    AE_LA16X4_IP(vxh, x_align, px);

    for (n = 0; n < N - 3; n += 4)
    {
        AE_MULAAAAQ16(vai, vxh, vxh);
        AE_LA16X4_IP(vxh, x_align, px);
    }

    N &= 3;
    if (N)
    {
        b = (1 << (4 - N)) - 1;
        AE_MOVT16X4(vxh, AE_ZERO16(), b);
        AE_MULAAAAQ16(vai, vxh, vxh);
    }

    vai = AE_SRAA64(vai, rsh);
    return_int64(vai);
#endif
} /* vec_power16x16() */
