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
    Blockwise Adaptive LMS Algorithm for Real Data, 24x24-bit
    C code optimized for HiFi3
	Integrit, 2006-2017
*/

/*-------------------------------------------------------------------------
  Blockwise Adaptive LMS Algorithm for Real Data
  Blockwise LMS algorithm performs filtering of reference samples x[N+M-1],
  computation of error e[N] over a block of input samples r[N] and makes
  blockwise update of IR to minimize the error output.
  Algorithm includes FIR filtering, calculation of correlation between the 
  error output e[N] and reference signal x[N+M-1] and IR taps update based
  on that correlation.
NOTES: 
  1. The algorithm must be provided with the normalization factor, which is
     the power of the reference signal times N - the number of samples in a
     data block. This can be calculated using the vec_power24x24() or 
     vec_power16x16() function. In order to avoid the saturation of the 
     normalization factor, it may be biased, i.e. shifted to the right.
     If it's the case, then the adaptation coefficient must be also shifted
     to the right by the same number of bit positions.
  2. This algorithm consumes less CPU cycles per block than single 
     sample algorithm at similar convergence rate.
  3. Right selection of N depends on the change rate of impulse response:
     on static or slow varying channels convergence rate depends on
     selected mu and M, but not on N.
  4. 16x16 routine may converge slower on small errors due to roundoff 
     errors. In that cases, 16x32 routine will give better results although
     convergence rate on bigger errors is the same.

  Precision: 
  16x16    16-bit coefficients, 16-bit data, 16-bit output
  24x24    24-bit coefficients, 24-bit data, 24-bit output
  16x32    32-bit coefficients, 16-bit data, 16-bit output
  32x32    32-bit coefficients, 32-bit data, 32-bit output
  f        floating point
  Input:
  h[M]     impulse response, Q15, Q31 or floating point
  r[N]	   input data vector (near end). First in time value is in 
           r[0], Q15, Q31 or floating point
  x[N+M-1] reference data vector (far end). First in time value is in x[0],  
           Q15, Q31 or floating point
  norm     normalization factor: power of signal multiplied by N, Q15, Q31  
           or floating point
           Fixed-point format for the 24x24-bit variant: Q(2*x-31-bias)
           Fixed-point format for the 32x16-bit variant: Q(2*x+1-bias)
  mu       adaptation coefficient (LMS step), Q(31-bias) or Q(15-bias)
  N        length of data block
  M        length of h
  Output:
  e[N]     estimated error, Q15, Q31 or floating point
  h[M]     updated impulse response, Q15, Q31 or floating point

  Restriction:
  x,r,h,e  should not overlap
  x,r,h,e  aligned on a 8-bytes boundary
  N,M      multiples of 8 and >0
-------------------------------------------------------------------------*/
#include "NatureDSP_types.h"
#include "NatureDSP_Signal_fir.h"
#include "common.h"

void fir_blms24x24( f24 * restrict e,
                    f24 * restrict h,
              const f24 * restrict r,
              const f24 * restrict x,
              f24 norm, f24 mu,
              int N, int M )
{

        ae_f24x2   * restrict pe = (      ae_f24x2 *)e;
        ae_f24x2   * restrict ph = (      ae_f24x2 *)h;
  const ae_f24x2   * restrict px = (const ae_f24x2 *)x;
  const ae_f24x2   * restrict pr = (const ae_f24x2 *)r;
  const ae_f24x2   * restrict px1;  
  const ae_f24x2   * restrict ph1;  
  int         n, m, s_exp, x_nsa, mu_exp;
  ae_f64      A0, A1, A2, A3 ,A4, A5, A6, A7;
  ae_int64    B0, B1, B2, B3, B4, B5, B6, B7, z1;
  ae_int24x2  d0, d1, d2, d3, d4, d5, d6, d7;
  ae_int32x2  zl, vxw, y1, vxa, y0, vrw, y2, y3;
  ae_f32x2    vxf, vaf, vbf, vcf, vdf;
  ae_f24x2    x0, x1, x2, x3, x4, x5, x6, x7, val, vbl;

  zl = AE_MOVDA32X2(0,1);//1
  z1 = AE_MOVINT64_FROMINT32X2(zl);

  NASSERT(e);
  NASSERT(h);
  NASSERT(r);
  NASSERT(x);
  NASSERT_ALIGN(e,8);
  NASSERT_ALIGN(h,8);
  NASSERT_ALIGN(r,8);
  NASSERT_ALIGN(x,8);
  NASSERT(N>0 && M>0);
  NASSERT(M%8==0 && N%8==0);
  //
  // Pass the reference signal x[] through the adaptive filter to obtain the
  // predicted signal and calculate the error, i.e. the distance to the 
  // actual input signal r[].
  //
  for ( n=0; n<N; n+=8 )
  {
    ph = (      ae_f24x2 *)h + (M/2)-1;
    AE_L32X2F24_IP(x0, px, +8);
    AE_L32X2F24_IP(x1, px, +8);
    AE_L32X2F24_IP(x2, px, +8);
    AE_L32X2F24_IP(x3, px, +8);
    px1 = (ae_f24x2 *) px;

    AE_L32X2F24_IP(x4, px1, +8);
    AE_L32X2F24_IP(x5, px1, +8);

    AE_L32X2F24_RIP(val, ph);
    AE_L32X2F24_RIP(vbl, ph);

    AE_MULFD24X2_FIR_H(A0,A1,x0,x1,val);
    AE_MULFD24X2_FIR_H(A2,A3,x1,x2,val);
    AE_MULFD24X2_FIR_H(A4,A5,x2,x3,val);
    AE_MULFD24X2_FIR_H(A6,A7,x3,x4,val);

    AE_MULAFD24X2_FIR_H(A0,A1,x1,x2,vbl);
    AE_MULAFD24X2_FIR_H(A2,A3,x2,x3,vbl);
    AE_MULAFD24X2_FIR_H(A4,A5,x3,x4,vbl);
    AE_MULAFD24X2_FIR_H(A6,A7,x4,x5,vbl);
    x0 = x2;
    x1 = x3;
    x2 = x4;
    x3 = x5;
    for (m=0; m<((M>>2)-1); m++ )
    {
      AE_L32X2F24_IP(x4, px1, +8);
      AE_L32X2F24_IP(x5, px1, +8);

      AE_L32X2F24_RIP(val, ph);
      AE_L32X2F24_RIP(vbl, ph);
      // Q16.47 <- Q16.47 + ( Q(23+8)*(23+8) - 16 + 1 )
      AE_MULAFD24X2_FIR_H(A0,A1,x0,x1,val);
      AE_MULAFD24X2_FIR_H(A2,A3,x1,x2,val);
      AE_MULAFD24X2_FIR_H(A4,A5,x2,x3,val);
      AE_MULAFD24X2_FIR_H(A6,A7,x3,x4,val);

      AE_MULAFD24X2_FIR_H(A0,A1,x1,x2,vbl);
      AE_MULAFD24X2_FIR_H(A2,A3,x2,x3,vbl);
      AE_MULAFD24X2_FIR_H(A4,A5,x3,x4,vbl);
      AE_MULAFD24X2_FIR_H(A6,A7,x4,x5,vbl);
      x0 = x2;
      x1 = x3;
      x2 = x4;
      x3 = x5;
    }
    // Q(23+8) <- [ Q16.47 - 24 w/ rounding and saturation ] + 8
    x0 = AE_ROUND24X2F48SASYM(A0,A1);
    x1 = AE_ROUND24X2F48SASYM(A2,A3);
    x2 = AE_ROUND24X2F48SASYM(A4,A5);
    x3 = AE_ROUND24X2F48SASYM(A6,A7);
    
    d0 = (x0);
    d1 = (x1);
    d2 = (x2);
    d3 = (x3);
        
    AE_L32X2F24_IP(x0, pr, +8);
    AE_L32X2F24_IP(x1, pr, +8);
    AE_L32X2F24_IP(x2, pr, +8);
    AE_L32X2F24_IP(x3, pr, +8);
    d4 = (x0);
    d5 = (x1);
    d6 = (x2);
    d7 = (x3);
    // Q(23+8) <- Q(23+8) - Q(23+8) w/ saturation
    d0 = AE_SUBSP24S(d4,d0);
    d1 = AE_SUBSP24S(d5,d1);
    d2 = AE_SUBSP24S(d6,d2);
    d3 = AE_SUBSP24S(d7,d3);
    x0 =  (d0);
    x1 =  (d1);
    x2 =  (d2);
    x3 =  (d3);
    AE_S32X2F24_IP(x0, pe, +8);
    AE_S32X2F24_IP(x1, pe, +8);
    AE_S32X2F24_IP(x2, pe, +8);
    AE_S32X2F24_IP(x3, pe, +8);
  }
  
  //
  // Compute the reciprocal for the normalization factor.
  //
  norm&=0xffffff00;
  mu  &=0xffffff00;
  s_exp = AE_NSAZ32_L((int32_t) norm);
  vxw = (int32_t) norm;
  vxw = AE_SLAA32(vxw, s_exp);
  /* reciprocal of Q31 in Q30: 6 LSB accuracy */
  {
      ae_int32x2 y,e,x,_1Q30;
      ae_f32x2 t;
      x=vxw;
       _1Q30=AE_MOVDA32((int)0x40000000);
      y=AE_SUB32 (_1Q30,x);/* no saturation here!!! */
      y=AE_ADD32 (_1Q30,y);
      y=AE_INT32X2_ADD32S(_1Q30,y);

      t=_1Q30; AE_MULSFP32X2RAS(t,x,y); e=t;
      t=AE_MULFP32X2RAS(e,y); e=t;
      y=AE_ADD32 (y,e);
      y=AE_ADD32 (y,e);

      t=_1Q30; AE_MULSFP32X2RAS(t,x,y); e=t;
      t=AE_MULFP32X2RAS(e,y); e=t;
      y=AE_ADD32 (y,e);
      y=AE_ADD32 (y,e);

      t=_1Q30; AE_MULSFP32X2RAS(t,x,y); e=t;
      t=AE_MULFP32X2RAS(e,y); e=t;
      y=AE_ADD32 (y,e);
      y=AE_ADD32 (y,e);
      vrw=y;
  }

  vxa = AE_SEL32_HH(vrw, vrw);
  x_nsa = AE_NSAZ32_L(vxa);
  vxa = AE_SLAA32(vxa, x_nsa);
  vxa = AE_SRAI32(vxa, 8);
  /*------------------------------------*/
  mu_exp = AE_NSAZ32_L((int32_t) mu);
  vxw = (int32_t) mu;
  vxw = AE_SRAA32(vxw,8);
  // Q(23+8+mu_exp) <- Q(23+8) + mu_exp
  vxw = AE_SLAA32S(vxw, mu_exp);
  // Q(22+8-s_exp+mu_exp) <- Q(22+8-s_exp)*Q(23+8+mu_exp) - 16 + 1 - 24 + 8
  B0 = AE_MUL32_LL(vxw,vxa);
  B0 = AE_SRAI64(B0,23);
  vxw = AE_MOVINT32X2_FROMINT64(B0);
  vxw = AE_SEL32_LL(vxw,vxw);
  
  s_exp -= mu_exp; // -31..30
  //
  // Calculate the cross-correlation between the error signal and the 
  // reference signal. Scale the result and update the estimation of the
  // impulse response.
  //
  z1 = AE_SLAA64(z1,(14-s_exp)); 
  ph = (      ae_f24x2 *)h;
  ph = ph + (M/2)-1;
  px = (const ae_f24x2 *)x;
  for ( m=0; m<M; m+=8 )
  {
    pe = (      ae_f24x2 *)e;
    AE_L32X2F24_IP(x0, px, +8);
    AE_L32X2F24_IP(x1, px, +8);
    AE_L32X2F24_IP(x2, px, +8);
    AE_L32X2F24_IP(x3, px, +8);
    px1 = (ae_f24x2 *) px;

    AE_L32X2F24_IP(x4, px1, +8);
    AE_L32X2F24_IP(x5, px1, +8);
    AE_L32X2F24_IP(val, pe, +8);
    AE_L32X2F24_IP(vbl, pe, +8);
    // Q16.47 <- Q16.47 + ( Q(23+8)*(23+8) - 16 + 1 )
    AE_MULFD24X2_FIR_H(A0,A1,x0,x1,val);
    AE_MULFD24X2_FIR_H(A2,A3,x1,x2,val);
    AE_MULFD24X2_FIR_H(A4,A5,x2,x3,val);
    AE_MULFD24X2_FIR_H(A6,A7,x3,x4,val);

    AE_MULAFD24X2_FIR_H(A0,A1,x1,x2,vbl);
    AE_MULAFD24X2_FIR_H(A2,A3,x2,x3,vbl);
    AE_MULAFD24X2_FIR_H(A4,A5,x3,x4,vbl);
    AE_MULAFD24X2_FIR_H(A6,A7,x4,x5,vbl);

    x0 = x2;
    x1 = x3;
    x2 = x4;
    x3 = x5;
    for (n=0; n<((N>>2)-1); n++ )
    {
      AE_L32X2F24_IP(x4, px1, +8);
      AE_L32X2F24_IP(x5, px1, +8);
      AE_L32X2F24_IP(val, pe, +8);
      AE_L32X2F24_IP(vbl, pe, +8);
      // Q16.47 <- Q16.47 + ( Q(23+8)*(23+8) - 16 + 1 )
      AE_MULAFD24X2_FIR_H(A0,A1,x0,x1,val);
      AE_MULAFD24X2_FIR_H(A2,A3,x1,x2,val);
      AE_MULAFD24X2_FIR_H(A4,A5,x2,x3,val);
      AE_MULAFD24X2_FIR_H(A6,A7,x3,x4,val);

      AE_MULAFD24X2_FIR_H(A0,A1,x1,x2,vbl);
      AE_MULAFD24X2_FIR_H(A2,A3,x2,x3,vbl);
      AE_MULAFD24X2_FIR_H(A4,A5,x3,x4,vbl);
      AE_MULAFD24X2_FIR_H(A6,A7,x4,x5,vbl);

      x0 = x2;
      x1 = x3;
      x2 = x4;
      x3 = x5;
    }
    // Q(23+8) <- [ Q16.47 - 24 w/ rounding and saturation ] + 8
    x0 = AE_ROUND24X2F48SASYM(A0,A1);
    x1 = AE_ROUND24X2F48SASYM(A2,A3);
    x2 = AE_ROUND24X2F48SASYM(A4,A5);
    x3 = AE_ROUND24X2F48SASYM(A6,A7);
    
    d0 = (x0);
    d1 = (x1);
    d2 = (x2);
    d3 = (x3);

    val = AE_MOVF24X2_FROMINT32X2(vxw);
    vxf = (val);
    
    vaf = (x0);
    vbf = (x1);
    vcf = (x2);
    vdf = (x3);
    // Q(17.46-s_exp) <- Q(23+8)*Q(22+8-s_exp) - 16 + 1
    A0 = AE_MULF32S_HH(vxf,vaf);
    A1 = AE_MULF32S_LL(vxf,vaf);
    A2 = AE_MULF32S_HH(vxf,vbf);
    A3 = AE_MULF32S_LL(vxf,vbf);
    A4 = AE_MULF32S_HH(vxf,vcf);
    A5 = AE_MULF32S_LL(vxf,vcf);
    A6 = AE_MULF32S_HH(vxf,vdf);
    A7 = AE_MULF32S_LL(vxf,vdf);
    B0 = (A0);
    B1 = (A1);
    B2 = (A2);
    B3 = (A3);
    B4 = (A4);
    B5 = (A5);
    B6 = (A6);
    B7 = (A7);

    ASSERT( s_exp >= -31 && s_exp < 31 );
    // Q31 <- Q(17.46-s_exp) + s_exp - 15 w/ rounding and saturation
    if ( s_exp < 15 )
    {
      B0 = AE_ADD64S(B0,z1);
      B1 = AE_ADD64S(B1,z1);
      B2 = AE_ADD64S(B2,z1);
      B3 = AE_ADD64S(B3,z1);
      B4 = AE_ADD64S(B4,z1);
      B5 = AE_ADD64S(B5,z1);
      B6 = AE_ADD64S(B6,z1);
      B7 = AE_ADD64S(B7,z1);
      B0 = AE_SRAA64(B0,(15 - s_exp));
      B1 = AE_SRAA64(B1,(15 - s_exp));
      B2 = AE_SRAA64(B2,(15 - s_exp));
      B3 = AE_SRAA64(B3,(15 - s_exp));
      B4 = AE_SRAA64(B4,(15 - s_exp));
      B5 = AE_SRAA64(B5,(15 - s_exp));
      B6 = AE_SRAA64(B6,(15 - s_exp));
      B7 = AE_SRAA64(B7,(15 - s_exp));

      y0 = AE_TRUNCA32X2F64S(B0,B1,24);
      y1 = AE_TRUNCA32X2F64S(B2,B3,24);
      y2 = AE_TRUNCA32X2F64S(B4,B5,24);
      y3 = AE_TRUNCA32X2F64S(B6,B7,24);

      x0 = AE_MOVF24X2_FROMINT32X2(y0);
      x1 = AE_MOVF24X2_FROMINT32X2(y1);
      x2 = AE_MOVF24X2_FROMINT32X2(y2);
      x3 = AE_MOVF24X2_FROMINT32X2(y3);
    }
    else
    {
      
      y0 = AE_TRUNCA32X2F64S(B0,B1,32);
      y1 = AE_TRUNCA32X2F64S(B2,B3,32);
      y2 = AE_TRUNCA32X2F64S(B4,B5,32);
      y3 = AE_TRUNCA32X2F64S(B6,B7,32);

      y0 = AE_SLAA32S(y0,s_exp - 15 - 8);
      y1 = AE_SLAA32S(y1,s_exp - 15 - 8);
      y2 = AE_SLAA32S(y2,s_exp - 15 - 8);
      y3 = AE_SLAA32S(y3,s_exp - 15 - 8);
      x0 = AE_MOVF24X2_FROMINT32X2(y0);
      x1 = AE_MOVF24X2_FROMINT32X2(y1);
      x2 = AE_MOVF24X2_FROMINT32X2(y2);
      x3 = AE_MOVF24X2_FROMINT32X2(y3);
      
    }
    
    ph1 = (ae_f24x2 *) ph;
    AE_L32X2F24_RIP(x4, ph1);
    AE_L32X2F24_RIP(x5, ph1);
    AE_L32X2F24_RIP(x6, ph1);
    AE_L32X2F24_RIP(x7, ph1);
    d0 = (x0);
    d1 = (x1);
    d2 = (x2);
    d3 = (x3);
    d4 = (x4);
    d5 = (x5);
    d6 = (x6);
    d7 = (x7);
    // Q31 <- Q31 + Q31 w/ saturation
    d0 = AE_ADDSP24S(d0,d4);
    d1 = AE_ADDSP24S(d1,d5);
    d2 = AE_ADDSP24S(d2,d6);
    d3 = AE_ADDSP24S(d3,d7);

    x0 =  (d0);
    x1 =  (d1);
    x2 =  (d2);
    x3 =  (d3);

    AE_S32X2F24_RIP(x0,ph); 
    AE_S32X2F24_RIP(x1,ph);
    AE_S32X2F24_RIP(x2,ph);
    AE_S32X2F24_RIP(x3,ph);
  }
} /* fir_blms24x24() */
