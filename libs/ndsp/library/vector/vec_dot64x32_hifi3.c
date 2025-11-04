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

/* Library API */
#include "NatureDSP_Signal_vector.h"
#include "common.h"


/*===========================================================================
  Vector matematics:
  vec_dot              Vector Dot Product
===========================================================================*/

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
int64_t vec_dot64x32 (const int64_t * restrict x,const int32_t * restrict y,int N)
{
    int n;
    const ae_int16x4 * restrict px;
    const ae_int32x2 * restrict py;
    ae_valign aly;
    ae_int64 ACChi, ACClo, ACC;
    ae_int32x2 yin;
    ae_int32x2 xh0, xl0, xh1, xl1;
    ae_int16x4 tx0, tx1, zero;
    ae_int64   A0,A1,A2,A3;
    xtbool cf, sat;

    if (N <= 0) return 0;

    px = (const ae_int16x4 *)x;
    py = (const ae_int32x2 *)y;
    aly = AE_LA64_PP(py);
    A0=A1=A2=A3 = AE_ZERO64();
    zero = AE_ZERO16();

    for (n=0; n<(N>>1); n++)
    {
        AE_L16X4_IP(tx0, px, sizeof(int64_t));
        AE_L16X4_IP(tx1, px, sizeof(int64_t));
        AE_LA32X2_IP(yin, aly, py);

        xl0 = AE_MOVINT32X2_FROMINT16X4(AE_SEL16_7362(zero, tx0));
        xh0 = AE_MOVINT32X2_FROMINT16X4(AE_SEL16_5146(zero, tx0));
        xl1 = AE_MOVINT32X2_FROMINT16X4(AE_SEL16_7362(zero, tx1));
        xh1 = AE_MOVINT32X2_FROMINT16X4(AE_SEL16_5146(zero, tx1));

        AE_MULA32X16_H0(A0,yin,tx0);/* Q46 <- Q15*Q31 */
        AE_MULA32_HH(A1,yin,xh0);   /* Q62 <- Q31*Q31 */
        AE_MULA32_HL(A2,yin,xl0);   /* Q78 <- Q47*Q31 */
        AE_MULA32_HH(A3,yin,xl0);   /* Q94 <- Q63*Q31 */

        AE_MULA32X16_L0(A0,yin,tx1);/* Q46 <- Q15*Q31 */
        AE_MULA32_LH(A1,yin,xh1);   /* Q62 <- Q31*Q31 */
        AE_MULA32_LL(A2,yin,xl1);   /* Q78 <- Q47*Q31 */
        AE_MULA32_LH(A3,yin,xl1);   /* Q94 <- Q63*Q31 */
    }
    if (N&1)
    {
        tx1 = AE_L16X4_I(px, 0);
        yin = AE_L32_I((const ae_int32 *)py, 0);
        xl1 = AE_MOVINT32X2_FROMINT16X4(AE_SEL16_7362(zero, tx1));
        xh1 = AE_MOVINT32X2_FROMINT16X4(AE_SEL16_5146(zero, tx1));

        AE_MULA32X16_L0(A0,yin,tx1);/* Q46 <- Q15*Q31 */
        AE_MULA32_LH(A1,yin,xh1);   /* Q62 <- Q31*Q31 */
        AE_MULA32_LL(A2,yin,xl1);   /* Q78 <- Q47*Q31 */
        AE_MULA32_LH(A3,yin,xl1);   /* Q94 <- Q63*Q31 */
    }

    ACChi  = A3>>63;  ACClo  = AE_AND64( A3    , MAX_INT64);
    ACChi += A2>>47;  ACClo += AE_AND64( A2<<16, MAX_INT64);
    cf = AE_INT64_LT(ACClo, 0);
    AE_MOVT64(ACChi, AE_ADD64(ACChi, 1), cf);
    AE_MOVT64(ACClo, AE_AND64(ACClo, MAX_INT64), cf);
    ACChi += A1>>31;  ACClo += AE_AND64( A1<<32, MAX_INT64);
    cf = AE_INT64_LT(ACClo, 0);
    AE_MOVT64(ACChi, AE_ADD64(ACChi, 1), cf);
    AE_MOVT64(ACClo, AE_AND64(ACClo, MAX_INT64), cf);
    ACChi += A0>>15;  ACClo += AE_AND64( A0<<48, MAX_INT64);
    cf = AE_INT64_LT(ACClo, 0);
    AE_MOVT64(ACChi, AE_ADD64(ACChi, 1), cf);
    AE_MOVT64(ACClo, AE_AND64(ACClo, MAX_INT64), cf);

    ACClo = AE_SRLI64(ACClo, 31);
    ACC = (ACChi<<32) | ACClo;
    sat = AE_INT64_GE( ACChi, 1LL<<31 );
    AE_MOVT64(ACC, MAX_INT64, sat);
    sat = AE_INT64_LT( ACChi,-2147483648LL );
    AE_MOVT64(ACC, MIN_INT64, sat);
    return ACC;
}
