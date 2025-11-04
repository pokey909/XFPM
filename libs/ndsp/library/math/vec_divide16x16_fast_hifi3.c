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
void vec_divide16x16_fast 
(
  int16_t *       restrict  frac,
  int16_t *       restrict  exp,
  const int16_t * restrict  x,
  const int16_t * restrict  y,
  int M)
{
    ae_int16x4 X,Y,Z,E;
    int n,N=M;
    xtbool4 mask5;
    const ae_int16x4 * restrict px;
    const ae_int16x4 * restrict py;
          ae_int16x4 * restrict pf;
    static const int16_t ALIGN(8) _0101[]={0,1,0,1};
    mask5=AE_EQ16(AE_L16X4_I((const ae_int16x4*)_0101,0),AE_ZERO16());
    NASSERT_ALIGN(x,8);
    NASSERT_ALIGN(y,8);
    NASSERT_ALIGN(frac,8);
    NASSERT(N%4==0);

    if(N<=0) return;
    /* take exponent and normalize inputs. Y is saved to the scratch */
    px=(const ae_int16x4 *)x;
    py=(const ae_int16x4 *)y;
    pf=(      ae_int16x4*)frac;
    for (n=0; n<(N>>2); n++)
    {
        ae_int16x4 X0,X1,X2,X3,Y0,Y1,Y2,Y3;
        xtbool4 sy;
        ae_int16x4 _0x4000=AE_MOVDA16(16384);
        int expx,expy;
        AE_L16X4_IP(X,px,sizeof(X));
        AE_L16X4_IP(Y,py,sizeof(Y));
        expx = AE_NSAZ16_0(X);
        expy = AE_NSAZ16_0(Y);
        X0=AE_SLAA16S(X,expx);
        Y0=AE_SLAA16S(Y,expy);
        exp[3]=(int16_t)(expy-expx+1);
        X=AE_SEL16_4321(X,X);
        Y=AE_SEL16_4321(Y,Y);

        expx = AE_NSAZ16_0(X);
        expy = AE_NSAZ16_0(Y);
        X1=AE_SLAA16S(X,expx);
        Y1=AE_SLAA16S(Y,expy);
        exp[2]=(int16_t)(expy-expx+1);
        X=AE_SEL16_4321(X,X);
        Y=AE_SEL16_4321(Y,Y);

        expx = AE_NSAZ16_0(X);
        expy = AE_NSAZ16_0(Y);
        X2=AE_SLAA16S(X,expx);
        Y2=AE_SLAA16S(Y,expy);
        exp[1]=(int16_t)(expy-expx+1);
        X=AE_SEL16_4321(X,X);
        Y=AE_SEL16_4321(Y,Y);

        expx = AE_NSAZ16_0(X);
        expy = AE_NSAZ16_0(Y);
        X3=AE_SLAA16S(X,expx);
        Y3=AE_SLAA16S(Y,expy);
        exp[0]=(int16_t)(expy-expx+1);
        exp+=4;

        X =AE_SEL16_6420(X2,X0);
        X1=AE_SEL16_6420(X3,X1);
        X1=AE_SEL16_6543(X1,X1);
        AE_MOVT16X4(X ,X1,mask5);
        Y =AE_SEL16_6420(Y2,Y0);
        Y1=AE_SEL16_6420(Y3,Y1);
        Y1=AE_SEL16_6543(Y1,Y1);
        AE_MOVT16X4(Y ,Y1,mask5);

        sy=AE_LT16(Y,AE_ZERO16());
        Y=AE_ABS16S_vector(Y);
        /* first approximation */
        Z=AE_SUB16(AE_MOVDA16((int16_t)47852),Y); 
            /* 3 iterations to achieve 1 LSB accuracy in mantissa */
        E=AE_SUB16(_0x4000,AE_MULFP16X4S_vector(Y,Z)); 
        E=AE_ADD16(E,E);
        Z=AE_ADD16S_vector(Z,AE_MULFP16X4S_vector(Z,E));
        E=AE_SUB16(_0x4000,AE_MULFP16X4S_vector(Y,Z)); 
        E=AE_ADD16(E,E);
        Z=AE_ADD16S_vector(Z,AE_MULFP16X4S_vector(Z,E));
        E=AE_SUB16(_0x4000,AE_MULFP16X4S_vector(Y,Z)); 
        E=AE_ADD16(E,E);
        Z=AE_ADD16S_vector(Z,AE_MULFP16X4S_vector(Z,E));
        /* restore original sign */
        Y=AE_NEG16S_vector(Z);
        AE_MOVT16X4(Z,Y,sy);
        /* multiply by X */
        Z=AE_MULFP16X4RAS(X,Z);
        AE_S16X4_IP(Z,pf,sizeof(Z));
    }
} /* vec_divide16x16_fast() */
