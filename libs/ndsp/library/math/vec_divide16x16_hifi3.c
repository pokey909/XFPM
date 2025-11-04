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
void vec_divide16x16 
(
  int16_t *       restrict  frac,
  int16_t *       restrict  exp,
  const int16_t * restrict  x,
  const int16_t * restrict  y,
  int M)
{
#define SCR_SZ (MAX_ALLOCA_SZ/sizeof(int16_t))
    int16_t ALIGN(16) scr[SCR_SZ];   /* local scratch */
    static const ALIGN(8) int16_t cnst0123[]={0,1,2,3};
    ae_int16x4 _0123;
    ae_int16x4 X,Y,Z,E;
    int n,N;
    const ae_int16   * restrict px;
    const ae_int16   * restrict py;
          ae_int16   * restrict pf;
          ae_int16   * restrict ps;
    const ae_int16x4 * restrict pfRd;
          ae_int16x4 * restrict pfWr;
    const ae_int16x4 * restrict psRd;
          ae_valign   wr_align,rd_align;

    while(M>0)
    {
        N=XT_MIN(SCR_SZ,M); /* size of portion */
        /* take exponent and normalize inputs. Y is saved to the scratch */
        px=(const ae_int16 *)x;
        py=(const ae_int16 *)y;
        pf=(      ae_int16 *)frac;
        ps=(      ae_int16 *)scr;
        for (n=0; n<N; n++)
        {
            int expx,expy;
            AE_L16_IP(X,px,sizeof(int16_t));
            AE_L16_IP(Y,py,sizeof(int16_t));
            expx = AE_NSAZ16_0(X);
            expy = AE_NSAZ16_0(Y);
            X=AE_SLAA16S(X,expx);
            Y=AE_SLAA16S(Y,expy);
            AE_S16_0_IP(X,pf,sizeof(int16_t));
            AE_S16_0_IP(Y,ps,sizeof(int16_t));
            exp[n]=(int16_t)(expy-expx+1);
        }
        __Pragma("no_reorder");
        pfRd = (const ae_int16x4*)frac;
        pfWr = (      ae_int16x4*)frac;
        psRd = (const ae_int16x4*)scr;
        wr_align=AE_ZALIGN64();
        rd_align=AE_LA64_PP(pfRd);
        _0123 = AE_L16X4_I((const ae_int16x4*)cnst0123,0);
        for(n=N; n>0; n-=4)
        {
            xtbool4 sy,wr_mask;
            ae_int16x4 _0x4000=AE_MOVDA16(16384);
            AE_L16X4_IP(Y,psRd,sizeof(Y));
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
            AE_LA16X4_IP(X,rd_align,pfRd);
            Z=AE_MULFP16X4RAS(X,Z);
            wr_mask=AE_LT16(_0123,AE_MOVDA16(n)); /* compute mask for last incomplete iteration */
            AE_MOVT16X4(X,Z,wr_mask);
            AE_SA16X4_IP(X,wr_align,pfWr);
        }
        AE_SA64POS_FP(wr_align,pfWr);
        /* process next portion */
        M-=N;
        x+=N;
        y+=N;
        frac+=N;
        exp+=N;
    }
} /* vec_divide16x16() */
