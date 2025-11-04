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
    Reciprocal
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/
/* Signal Processing Library API. */
#include "NatureDSP_Signal_math.h"
#include "common.h"

/*-------------------------------------------------------------------------
  Reciprocal on Q31/Q15 Numbers
  These routines return the fractional and exponential portion of the 
  reciprocal of a vector x of Q31 or Q15 numbers. Since the reciprocal is 
  always greater than 1, it returns fractional portion frac in Q(31-exp) 
  or Q(15-exp) format and exponent exp so true reciprocal value in the 
  Q0.31/Q0.15 may be found by shifting fractional part left by exponent 
  value.

  Mantissa accuracy is 1 LSB, so relative accuracy is:
  vec_recip16x16, scl_recip16x16                   6.2e-5 
  vec_recip24x24, scl_recip32x32, scl_recip24x24   2.4e-7 
  vec_recip32x32                                   9.2e-10

  Precision: 
  32x32  32-bit input, 32-bit output. 
  24x24  24-bit input, 24-bit output. 
  16x16  16-bit input, 16-bit output. 

  Input:
  x[N]    input data, Q31 or Q15
  N       length of vectors

  Output:
  frac[N] fractional part of result, Q(31-exp) or Q(15-exp)
  exp[N]  exponent of result 

  Restriction:
  x,frac,exp should not overlap

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  frac,x - aligned on 8-byte boundary
  N      - multiple of 4 and >4

  Scalar versions:
  ----------------
  Return packed value: 
  scl_recip24x24(),scl_recip32x32():
  bits 23:0 fractional part
  bits 31:24 exponent
  scl_recip16x16():
  bits 15:0 fractional part
  bits 31:16 exponent
-------------------------------------------------------------------------*/
uint32_t scl_recip32x32 (int32_t x)
{
    uint32_t res;
    /* compute exponent and normalize inputs */
    xtbool2 sx,isZero;
    ae_f32x2 t;
    ae_int32x2 X,Y,E,_0x40000000;
    unsigned nsa;
    X=AE_MOVDA32(x);
    isZero=AE_EQ32(X,AE_ZERO32());
    nsa=AE_NSAZ32_L(X);
    XT_MOVT(nsa,0x20,xtbool2_extract_0(isZero));  /* return biggest exponent for zero input */
    X=AE_SLAA32S(X,nsa);
    res=((uint32_t)(nsa+1))<<24;

    _0x40000000=AE_MOVDA32(0x40000000);
    sx=AE_LT32(X,AE_ZERO32());
    X=AE_INT32X2_ABS32S(X);
    AE_MOVT32X2(X,AE_MOVDA32(0x7fffffff),isZero);
    Y=AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000),X); /* first approximation */
    /* 3 iterations to achieve 1 LSB accuracy in mantissa */
    t=_0x40000000; AE_MULSFP32X2RAS(t,X,Y); E=t;
    E=AE_ADD32(E,E);
    t=Y; AE_MULAFP32X2RAS(t,Y,E); Y=t;
    t=_0x40000000; AE_MULSFP32X2RAS(t,X,Y); E=t;
    E=AE_ADD32(E,E);
    t=Y; AE_MULAFP32X2RAS(t,Y,E); Y=t;
    t=_0x40000000; AE_MULSFP32X2RAS(t,X,Y); E=t;
    E=AE_ADD32(E,E);
    t=Y; AE_MULAFP32X2RAS(t,Y,E); Y=t;

    X=AE_INT32X2_NEG32S(Y); /* restore original sign */
    AE_MOVT32X2(Y,X,sx);

    res|=((uint32_t)AE_MOVAD32_L(Y))>>8;
    return res;
} /* scl_recip32x32() */

