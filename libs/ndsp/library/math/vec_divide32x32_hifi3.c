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
    Division
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/
/* Signal Processing Library API. */
#include "NatureDSP_Signal_math.h"
#include "common.h"

/*-------------------------------------------------------------------------
  Division
  These routines perform pair-wise division of vectors written in Q31 or Q15 
  format. They return the fractional and exponential portion of the division 
  result. Since the division may generate result greater than 1, it returns 
  fractional portion frac in Q(31-exp) or Q(15-exp) format and exponent 
  exp so true division result in the Q0.31 may be found by shifting 
  fractional part left by exponent value.
  Additional routine makes integer division of 64-bit number to 32-bit 
  denominator forming 32-bit result. If result is overflown, 0x7fffffff 
  or 0x80000000 is returned depending on the signs of inputs.
  For division to 0, the result is not defined.

  Two versions of routines are available: regular versions (vec_divide64x32i,
  vec_divide32x32, vec_divide24x24, vec_divide16x16) work with arbitrary
  arguments, faster versions (vec_divide32x32_fast, vec_divide24x24_fast, 
  vec_divide16x16_fast) apply some restrictions.

  Accuracy is measured as accuracy of fractional part (mantissa):
  vec_divide64x32i, scl_divide64x32                      :  1 LSB   
  vec_divide32x32, vec_divide32x32_fast                  :  2 LSB (1.8e-9) 
  scl_divide32x32, vec_divide24x24, scl_divide24x24      :  2 LSB (4.8e-7) 
  vec_divide16x16, scl_divide16x16, vec_divide16x16_fast :  2 LSB (1.2e-4)

  Precision: 
  64x32i integer division, 64-bit nominator, 32-bit denominator, 32-bit output. 
  32x32  fractional division, 32-bit inputs, 32-bit output. 
  24x24  fractional division, 24-bit inputs, 24-bit output. 
  16x16  fractional division, 16-bit inputs, 16-bit output. 

  Input:
  x[N]    nominator, 64-bit integer, Q31 or Q15
  y[N]    denominator, 32-bit integer, Q31 or Q15
  N       length of vectors
  Output:
  frac[N] fractional parts of result, Q(31-exp) or Q(15-exp)
  exp[N]  exponents of result 

  Restriction:
  For regular versions (vec_divide64x32i, vec_divide32x32,
  vec_divide24x24, vec_divide16x16) :
  x,y,frac,exp should not overlap

  For faster versions (vec_divide32x32_fast, vec_divide24x24_fast, 
  vec_divide16x16_fast) :
  x,y,frac,exp  should not overlap
  x,y,frac      to be aligned by 8-byte boundary, N - multiple of 4.

  Scalar versions:
  ----------------
  scl_divide64x32(): integer remainder
  Return packed value: 
  scl_divide24x24(),scl_divide32x32():
  bits 23:0 fractional part
  bits 31:24 exponent
  scl_divide16x16():
  bits 15:0 fractional part
  bits 31:16 exponent
-------------------------------------------------------------------------*/
void vec_divide32x32 
(
  int32_t * restrict        frac,
  int16_t *                 exp,
  const int32_t * restrict  x,
  const int32_t * restrict  y,
  int                       M
)
{
    ae_int32x2 X,Y,Z,E;
    int n,N=M;
    const ae_int32x2  * restrict px;
    const ae_int32x2  * restrict py;
          ae_int32x2 * restrict pfWr;
          ae_valign   wr_align;
          ae_valign   x_align ,y_align ;
    if(N<=0) return;
    /* take exponent and normalize inputs. Y is saved to the scratch */
    px=(const ae_int32x2*)x;
    py=(const ae_int32x2*)y;
    pfWr = (      ae_int32x2*)frac;
    wr_align=AE_ZALIGN64();    
    x_align =AE_LA64_PP(px);
    y_align =AE_LA64_PP(py);

    ae_int32x2 Xl,Xh,Yl,Yh;
    ae_f32x2 t;
    xtbool2 sy;
    ae_int32x2 _0x40000000=AE_MOVDA32(0x40000000);
    int expx,expy;

    for (n=0; n<(N&(~1))/2; n++)
    {        
        /* normalization */
        AE_LA32X2_IP(X,x_align,px);
        AE_LA32X2_IP(Y,y_align,py);
        expx = AE_NSAZ32_L(X);
        expy = AE_NSAZ32_L(Y);
        Xl=AE_SLAA32S(X,expx);
        Yl=AE_SLAA32S(Y,expy);
        exp[1]=(int16_t)(expy-expx+1);
        X=AE_SEL32_LH(X,X);
        Y=AE_SEL32_LH(Y,Y);
        expx = AE_NSAZ32_L(X);
        expy = AE_NSAZ32_L(Y);
        Xh=AE_SLAA32S(X,expx);
        Yh=AE_SLAA32S(Y,expy);
        exp[0]=(int16_t)(expy-expx+1);
        X=AE_SEL32_LL(Xh,Xl);
        Y=AE_SEL32_LL(Yh,Yl);
        exp+=2;

        sy=AE_LT32(Y,AE_ZERO32());
        Y=AE_INT32X2_ABS32S(Y);
        /* first approximation */
        Z=AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000),Y); 
            /* 4 iterations to achieve 1 LSB accuracy in mantissa */
        t=_0x40000000; AE_MULSFP32X2RAS(t,Y,Z); E=t;
        E=AE_ADD32(E,E);
        t=Z; AE_MULAFP32X2RAS(t,Z,E); Z=t;
        t=_0x40000000; AE_MULSFP32X2RAS(t,Y,Z); E=t;
        E=AE_ADD32(E,E);
        t=Z; AE_MULAFP32X2RAS(t,Z,E); Z=t;
        t=_0x40000000; AE_MULSFP32X2RAS(t,Y,Z); E=t;
        E=AE_ADD32(E,E);
        t=Z; AE_MULAFP32X2RAS(t,Z,E); Z=t;
        t=_0x40000000; AE_MULSFP32X2RAS(t,Y,Z); E=t;
        E=AE_ADD32(E,E);
        t=Z; AE_MULAFP32X2RAS(t,Z,E); Z=t;
        /* restore original sign */
        Y=AE_INT32X2_NEG32S(Z);
        AE_MOVT32X2(Z,Y,sy);
        /* multiply by X */
        Z=AE_MULFP32X2RAS(X,Z);
        AE_SA32X2_IP(Z,wr_align,pfWr);
    }
    if( N&1 )
    {        
        /* normalization */
        AE_LA32X2_IP(X,x_align,px);
        AE_LA32X2_IP(Y,y_align,py);
        expx = AE_NSAZ32_L(X);
        expy = AE_NSAZ32_L(Y);
        Xl=AE_SLAA32S(X,expx);
        Yl=AE_SLAA32S(Y,expy);
        X=AE_SEL32_LH(X,X);
        Y=AE_SEL32_LH(Y,Y);
        expx = AE_NSAZ32_L(X);
        expy = AE_NSAZ32_L(Y);
        Xh=AE_SLAA32S(X,expx);
        Yh=AE_SLAA32S(Y,expy);
        exp[0]=(int16_t)(expy-expx+1);
        X=AE_SEL32_LL(Xh,Xl);
        Y=AE_SEL32_LL(Yh,Yl);
        exp+=2;

        sy=AE_LT32(Y,AE_ZERO32());
        Y=AE_INT32X2_ABS32S(Y);
        /* first approximation */
        Z=AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000),Y); 
            /* 4 iterations to achieve 1 LSB accuracy in mantissa */
        t=_0x40000000; AE_MULSFP32X2RAS(t,Y,Z); E=t;
        E=AE_ADD32(E,E);
        t=Z; AE_MULAFP32X2RAS(t,Z,E); Z=t;
        t=_0x40000000; AE_MULSFP32X2RAS(t,Y,Z); E=t;
        E=AE_ADD32(E,E);
        t=Z; AE_MULAFP32X2RAS(t,Z,E); Z=t;
        t=_0x40000000; AE_MULSFP32X2RAS(t,Y,Z); E=t;
        E=AE_ADD32(E,E);
        t=Z; AE_MULAFP32X2RAS(t,Z,E); Z=t;
        t=_0x40000000; AE_MULSFP32X2RAS(t,Y,Z); E=t;
        E=AE_ADD32(E,E);
        t=Z; AE_MULAFP32X2RAS(t,Z,E); Z=t;
        /* restore original sign */
        Y=AE_INT32X2_NEG32S(Z);
        AE_MOVT32X2(Z,Y,sy);
        /* multiply by X */
        Z=AE_MULFP32X2RAS(X,Z);
        X=AE_SEL32_LH(Z,Z);
        AE_S32_L_I(X,(ae_int32 *)pfWr,0);
    }
    AE_SA64POS_FP(wr_align,pfWr);
} /* vec_divide32x32() */
