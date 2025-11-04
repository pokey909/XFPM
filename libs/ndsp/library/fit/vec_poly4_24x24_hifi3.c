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
  NatureDSP Signal Processing Library. Fitting and Interpolation Routines
    Polynomial approximation
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/
#include "common.h"
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fit.h"
/*-------------------------------------------------------------------------
  Polynomial approximation
  Functions calculate polynomial approximation for all values from given 
  vector. Fixed point functions take polynomial coefficients in Q31 precision. 
  NOTE:
  approximation is calculated like Taylor series that is why overflow may 
  potentially occur if cumulative sum of coefficients given from the last 
  to the first coefficient is bigger that 1. To avoid this negative effect,
  all the coefficients may be scaled down and result will be shifted left 
  after all intermediate computations.

  Precision: 
  24x24  24-bit inputs, 24-bit coefficients, 24-bit output. 
  32x32  32-bit inputs, 32-bit coefficients, 32-bit output.
  f      floating point

  Input:
  x[N]    input data, Q31 or floating point
  N       length of vector
  lsh     additional left shift for result
  c[M+1]  coefficients (M=4 coefficients for vec_poly4_xxx 
          and M=8 coefficients for vec_poly8_xxx), Q31 or floating point
  Output:			
  z[N]    result, Q31 or floating point

  Restriction:
  x,c,z should not overlap
  lsh   should be in range 0...31

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  x aligned on 8-byte boundary
  N   - multiple of 2
-------------------------------------------------------------------------*/
void vec_poly4_24x24 (f24 * restrict z, const f24 * restrict x, const f24 * restrict c, int lsh, int N )
{
    ae_int32x2 ALIGN(16) temp[3];
    int         n;
    const ae_f24x2 * px = (const ae_f24x2  *)x;
    const ae_int32 * pc = (const ae_int32  *)c;
          ae_f32x2 * pz = (       ae_f32x2 *)z;

    ae_valign      x_align, z_align, c_align;
    ae_int32x2    c0f, c1f, c2f, c3f, c4f, s, vxf;
    ae_f32x2 t;
    ae_f24x2 tt;

    NASSERT(x);
    NASSERT(c);
    NASSERT(z);
    NASSERT(lsh>=0 && lsh<=31);

    if (N <= 0) return;
    c_align = AE_LA64_PP(pc);
    AE_LA32X2F24_IP(tt, c_align, castxcc(const ae_f24x2,pc)); s=AE_SLAI32(tt, 8); AE_S32X2_I(s,temp,0);
    AE_LA32X2F24_IP(tt, c_align, castxcc(const ae_f24x2,pc)); s=AE_SLAI32(tt, 8); AE_S32X2_I(s,temp,1*sizeof(ae_int32x2));
    AE_LA32X2F24_IP(tt, c_align, castxcc(const ae_f24x2,pc)); s=AE_SLAI32(tt, 8); AE_S32X2_I(s,temp,2*sizeof(ae_int32x2));
    pc=(const ae_int32  *)temp;
    __Pragma("no_reorder");
  
    x_align = AE_LA64_PP(px);
    z_align = AE_ZALIGN64();
    vxf=0; c4f=0;
    c3f=AE_L32_I (pc, 3*sizeof(*pc)); t=c3f; AE_MULAFP32X2RAS(t, vxf, c4f); s=(t);
    c2f=AE_L32_I (pc, 2*sizeof(*pc)); t=c2f; AE_MULAFP32X2RAS(t, vxf, s  ); s=(t);
    c4f=AE_L32_X (pc, 4*sizeof(*pc));
    WUR_AE_SAR(lsh);
    for (n=0; n<(N&~1); n+=2)
    { 
        AE_LA32X2F24_IP(tt, x_align, px); vxf=tt; vxf=AE_SLAI32(vxf, 8);
                                          t=c3f; AE_MULAFP32X2RAS(t, vxf, c4f); s=(t);
                                          t=c2f; AE_MULAFP32X2RAS(t, vxf, s  ); s=(t);
        c1f=AE_L32_I (pc, 1*sizeof(*pc)); t=c1f; AE_MULAFP32X2RAS(t, vxf, s  ); s=(t);
        AE_L32_XP(c0f,pc, 0*sizeof(*pc)); t=c0f; AE_MULAFP32X2RAS(t, vxf, s  ); 
        t=AE_SLAS32S(t);
        AE_SA32X2_IP(t, z_align, pz);//put answer
    }
    AE_SA64POS_FP(z_align, pz);
    if (N&1)
    { 
        AE_LA32X2F24_IP(tt, x_align, px); vxf=tt; vxf=AE_SLAI32(vxf, 8);
                                          t=c3f; AE_MULAFP32X2RAS(t, vxf, c4f); s=(t);
                                          t=c2f; AE_MULAFP32X2RAS(t, vxf, s  ); s=(t);
        c1f=AE_L32_I (pc, 1*sizeof(*pc)); t=c1f; AE_MULAFP32X2RAS(t, vxf, s  ); s=(t);
        AE_L32_XP(c0f,pc, 0*sizeof(*pc)); t=c0f; AE_MULAFP32X2RAS(t, vxf, s  ); 
        s=AE_SLAS32S(t);
        s=AE_SEL32_HH(s,s);
        AE_S32_L_I(s, (ae_int32 *)pz, 0);//put answer
    }
} /* vec_poly4_24x24() */
