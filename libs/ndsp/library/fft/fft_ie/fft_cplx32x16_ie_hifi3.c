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
	NatureDSP Signal Processing Library. FFT part
    FFT on Complex Data with Optimized Memory Usage
    C code optimized for HiFi3
	Integrit, 2006-2017
*/

#include "NatureDSP_Signal_fft.h"
#include "common.h"

extern int fft_stage_last_ie( int32_t *x, 
                               int32_t *y, 
                               int N); 
/*-------------------------------------------------------------------------
  FFT on Complex Data with Optimized Memory Usage
  These functions make FFT on complex data with optimized memory usage.
  Scaling: 
      +-------------------+----------------------------------------+
      |      Function     |           Scaling options              |
      +-------------------+----------------------------------------+
      |  fft_cplx16x16_ie |  2 - 16-bit dynamic scaling            | 
      |  fft_cplx24x24_ie |  3 - fixed scaling before each stage   | 
      |  fft_cplx32x16_ie |  3 - fixed scaling before each stage   | 
      |  fft_cplx32x32_ie |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      +-------------------+----------------------------------------+
  NOTES:
  1. Bit-reversing reordering is done here.
  2. FFT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED after 
     the call.
  3. FFT of size N may be supplied with constant data
     (twiddle factors) of a larger-sized FFT = N*twdstep.

  Precision: 
  16x16_ie      16-bit input/outputs, 16-bit twiddles
  24x24_ie      24-bit input/outputs, 24-bit twiddles
  32x16_ie      32-bit input/outputs, 16-bit twiddles
  32x32_ie      32-bit input/outputs, 32-bit twiddles
  f_ie          floating point
 
  Input:
  x[N]                  complex input signal. Real and imaginary data 
                        are interleaved and real data goes first
  twd[N*twdstep*3/4]    twiddle factor table of a complex-valued FFT of 
                        size N*twdstep
  N                     FFT size
  twdstep               twiddle step 
  scalingOpt            scaling option (see table above), not applicable
                        to the floating point function 
  Output:
  y[N]                  output spectrum. Real and imaginary data are 
                        interleaved and real data goes first

  Returned value: total number of right shifts occurred during scaling 
                  procedure. Floating point function always return 0.

  Restrictions:
  x,y   should not overlap
  x,y   aligned on 8-bytes boundary
-------------------------------------------------------------------------*/
int fft_cplx32x16_ie( complex_fract32* y,complex_fract32* x, const complex_fract16* twd, int twdstep, int N, int scalingOpt)
#if 0
{
    int shift = 0;
    int M = twdstep*N;
    int s;
    int tw_step = twdstep;
    int stride = N>>2;     

    const ae_int32x2 * restrict px;
          ae_int32x2 * restrict py;
    const ae_int32   * restrict p16tw1; 
 
    ae_int32x2  vA0, vA1, vA2, vA3, 
                vB0, vB1, vB2, vB3, 
                vC0, vC1, vC2, vC3;
    int i; 
    int32_t acc_inc=0;

    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT(scalingOpt==3);
    NASSERT(N==128||N==256||N==512||N==1024);
    WUR_AE_CBEGIN0((unsigned)x); 
    WUR_AE_CEND0((unsigned)(uintptr_t)(x+(N-1)));     
    s = 3; 

    while( stride > 1 )
    {
        ae_int32x2 t1, t2, t3; 
        ae_f16x4 t1_f16x4; 
        ae_f16x4 t2_f16x4; 
        ae_f16x4 t3_f16x4; 
        int acc = 0; 
        WUR_AE_SAR(s); 
        p16tw1 = (ae_int32*)  twd; 
        shift += s;
        s = 2; 
        px = (const ae_int32x2*)x; 
        py = (ae_int32x2*)x; 
  
        for (i=0; i<(N>>2); i++)
        {
            {
                int offset_inc = 0; 
                acc+=acc_inc;
                XT_MOVEQZ(offset_inc, tw_step*4, acc); 
                t2 =AE_L32_X ( p16tw1, 1*M);
                t3 =AE_L32_X ( p16tw1, 2*M);
                AE_L32_XP( t1 ,p16tw1, offset_inc);
            }
//            vA3 = AE_L32X2_X(px, 3*8*stride); 
//            vA2 = AE_L32X2_X(px, 2*8*stride); 
//            vA1 = AE_L32X2_X(px, 1*8*stride); 
//            AE_L32X2_XC(vA0, px, 4*8*stride); 
            px=(const ae_int32x2*)py;
            vA3 = AE_L32X2_X(px, 3*8*stride); 
            vA2 = AE_L32X2_X(px, 2*8*stride); 
            vA1 = AE_L32X2_X(px, 1*8*stride); 
            vA0=AE_L32X2_I(px, 0); 

            vA3 = AE_SRAS32(vA3); 
            vA2 = AE_SRAS32(vA2); 
            vA1 = AE_SRAS32(vA1); 
            vA0 = AE_SRAS32(vA0); 

            vB0 = AE_ADD32S(vA0, vA2);
            vB2 = AE_SUB32S(vA0, vA2);
            vB1 = AE_ADD32S(vA1, vA3);
            vB3 = AE_SUB32S(vA1, vA3);
            
            vB3 = AE_SEL32_LH(vB3, vB3);

            vC0 = AE_ADD32S(vB0, vB1);
            vC2 = AE_SUB32S(vB0, vB1);
            vC3 = AE_SUBADD32S(vB2, vB3);
            vC1 = AE_ADDSUB32S(vB2, vB3);

            t1_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
            t2_f16x4 = AE_MOVF16X4_FROMINT32X2(t2);
            t3_f16x4 = AE_MOVF16X4_FROMINT32X2(t3);

            vC1 = AE_MOVINT32X2_FROMF32X2(AE_MULFC32X16RAS_L(AE_MOVF32X2_FROMINT32X2(vC1), t1_f16x4)); 
            vC2 = AE_MOVINT32X2_FROMF32X2(AE_MULFC32X16RAS_L(AE_MOVF32X2_FROMINT32X2(vC2), t2_f16x4)); 
            vC3 = AE_MOVINT32X2_FROMF32X2(AE_MULFC32X16RAS_L(AE_MOVF32X2_FROMINT32X2(vC3), t3_f16x4)); 

            AE_S32X2_X (vC3, py, 3*8*stride);
            AE_S32X2_X (vC2, py, 1*8*stride);
            AE_S32X2_X (vC1, py, 2*8*stride);
            AE_S32X2_XC(vC0, py, 4*8*stride);

        }
        stride>>=2;  
        tw_step<<=2;
        acc_inc  = (acc_inc==0)? 0x40000000: acc_inc>>2; 
    }
    shift += fft_stage_last_ie((int32_t*)x, (int32_t*)y, N); 
    return shift;
} /* fft_cplx32x16_ie() */
#else
{
    int shift = 0;
    int M = twdstep*N;
    int tw_step = twdstep;
    int stride = N>>2;     

    const ae_int32x2 * restrict px;
          ae_int32x2 * restrict py;
    const ae_int32x2 * restrict px1;
          ae_int32x2 * restrict py1;
    const ae_int32   * restrict p16tw1; 
 
    ae_int32x2  vA0, vA1, vA2, vA3, 
                vB0, vB1, vB2, vB3, 
                vC0, vC1, vC2, vC3;
    int i; 
    uint32_t acc_inc=0;

    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT(scalingOpt==3);
    NASSERT(N==128||N==256||N==512||N==1024);
    WUR_AE_CBEGIN0((unsigned)x); 
    WUR_AE_CEND0((unsigned)(uintptr_t)(x+(N-1)));     
    // on the first stage we are reading all the twiddles, on all subsequent 
    // stages we will read just each of 2
    {
        ae_int32x2 t1, t2, t3; 
        ae_f16x4 t1_f16x4; 
        ae_f16x4 t2_f16x4; 
        ae_f16x4 t3_f16x4; 
        WUR_AE_SAR(3); 
        p16tw1 = (ae_int32*)  twd; 
        shift += 3;
        px = (const ae_int32x2*)x; 
        py = (ae_int32x2*)x; 
  
        for (i=0; i<(N>>2); i++)
        {
            t2 =AE_L32_X ( p16tw1, 1*M);
            t3 =AE_L32_X ( p16tw1, 2*M);
            AE_L32_XP( t1 ,p16tw1, tw_step*4);
//            vA3 = AE_L32X2_X(px, 3*8*stride); 
//            vA2 = AE_L32X2_X(px, 2*8*stride); 
//            vA1 = AE_L32X2_X(px, 1*8*stride); 
//            AE_L32X2_XC(vA0, px, 4*8*stride); 
            px=(const ae_int32x2*)py;
            vA3 = AE_L32X2_X(px, 3*8*stride); 
            vA2 = AE_L32X2_X(px, 2*8*stride); 
            vA1 = AE_L32X2_X(px, 1*8*stride); 
            vA0=AE_L32X2_I(px, 0); 

            vA3 = AE_SRAS32(vA3); 
            vA2 = AE_SRAS32(vA2); 
            vA1 = AE_SRAS32(vA1); 
            vA0 = AE_SRAS32(vA0); 

            vB0 = AE_ADD32S(vA0, vA2);
            vB2 = AE_SUB32S(vA0, vA2);
            vB1 = AE_ADD32S(vA1, vA3);
            vB3 = AE_SUB32S(vA1, vA3);
            
            vB3 = AE_SEL32_LH(vB3, vB3);

            vC0 = AE_ADD32S(vB0, vB1);
            vC2 = AE_SUB32S(vB0, vB1);
            vC3 = AE_SUBADD32S(vB2, vB3);
            vC1 = AE_ADDSUB32S(vB2, vB3);

            t1_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
            t2_f16x4 = AE_MOVF16X4_FROMINT32X2(t2);
            t3_f16x4 = AE_MOVF16X4_FROMINT32X2(t3);

            vC1 = AE_MOVINT32X2_FROMF32X2(AE_MULFC32X16RAS_L(AE_MOVF32X2_FROMINT32X2(vC1), t1_f16x4)); 
            vC2 = AE_MOVINT32X2_FROMF32X2(AE_MULFC32X16RAS_L(AE_MOVF32X2_FROMINT32X2(vC2), t2_f16x4)); 
            vC3 = AE_MOVINT32X2_FROMF32X2(AE_MULFC32X16RAS_L(AE_MOVF32X2_FROMINT32X2(vC3), t3_f16x4)); 

            AE_S32X2_X (vC3, py, 3*8*stride);
            AE_S32X2_X (vC2, py, 1*8*stride);
            AE_S32X2_X (vC1, py, 2*8*stride);
            AE_S32X2_XC(vC0, py, 4*8*stride);

        }
        stride>>=2;  
        tw_step<<=2;
    }
    // intermediate stages
    acc_inc  = 0x80000000; 
    while( stride > 1 )
    {
        ae_int32x2 t1, t2, t3; 
        ae_f16x4 t1_f16x4; 
        ae_f16x4 t2_f16x4; 
        ae_f16x4 t3_f16x4; 
        uint32_t acc = 0; 
        __Pragma("no_reorder")
        WUR_AE_SAR(2); 
        p16tw1 = (ae_int32*)  twd; 
        shift += 2;
        py = (ae_int32x2*)x; 
        py1=py+4*stride;
        for (i=0; i<(N>>3); i++)
        {
            {
                int offset_inc = 0; 
                acc+=acc_inc;
                XT_MOVEQZ(offset_inc, tw_step*4, acc); 
                t2 =AE_L32_X ( p16tw1, 1*M);
                t3 =AE_L32_X ( p16tw1, 2*M);
                AE_L32_XP( t1 ,p16tw1, offset_inc);
            }
            px=(const ae_int32x2*)py;
            vA3 = AE_L32X2_X(px, 3*8*stride); 
            vA2 = AE_L32X2_X(px, 2*8*stride); 
            vA1 = AE_L32X2_X(px, 1*8*stride); 
            vA0=AE_L32X2_I(px, 0); 

            vA3 = AE_SRAS32(vA3); 
            vA2 = AE_SRAS32(vA2); 
            vA1 = AE_SRAS32(vA1); 
            vA0 = AE_SRAS32(vA0); 

            vB0 = AE_ADD32S(vA0, vA2);
            vB2 = AE_SUB32S(vA0, vA2);
            vB1 = AE_ADD32S(vA1, vA3);
            vB3 = AE_SUB32S(vA1, vA3);
            
            vB3 = AE_SEL32_LH(vB3, vB3);

            vC0 = AE_ADD32S(vB0, vB1);
            vC2 = AE_SUB32S(vB0, vB1);
            vC3 = AE_SUBADD32S(vB2, vB3);
            vC1 = AE_ADDSUB32S(vB2, vB3);

            t1_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
            t2_f16x4 = AE_MOVF16X4_FROMINT32X2(t2);
            t3_f16x4 = AE_MOVF16X4_FROMINT32X2(t3);

            vC1 = AE_MOVINT32X2_FROMF32X2(AE_MULFC32X16RAS_L(AE_MOVF32X2_FROMINT32X2(vC1), t1_f16x4)); 
            vC2 = AE_MOVINT32X2_FROMF32X2(AE_MULFC32X16RAS_L(AE_MOVF32X2_FROMINT32X2(vC2), t2_f16x4)); 
            vC3 = AE_MOVINT32X2_FROMF32X2(AE_MULFC32X16RAS_L(AE_MOVF32X2_FROMINT32X2(vC3), t3_f16x4)); 

            AE_S32X2_X (vC3, py, 3*8*stride);
            AE_S32X2_X (vC2, py, 1*8*stride);
            AE_S32X2_X (vC1, py, 2*8*stride);
            AE_S32X2_XC(vC0, py, 2*4*8*stride);

            px1=(const ae_int32x2*)py1;
            vA3 = AE_L32X2_X(px1, 3*8*stride); 
            vA2 = AE_L32X2_X(px1, 2*8*stride); 
            vA1 = AE_L32X2_X(px1, 1*8*stride); 
            vA0=AE_L32X2_I(px1, 0); 

            vA3 = AE_SRAS32(vA3); 
            vA2 = AE_SRAS32(vA2); 
            vA1 = AE_SRAS32(vA1); 
            vA0 = AE_SRAS32(vA0); 

            vB0 = AE_ADD32S(vA0, vA2);
            vB2 = AE_SUB32S(vA0, vA2);
            vB1 = AE_ADD32S(vA1, vA3);
            vB3 = AE_SUB32S(vA1, vA3);
            
            vB3 = AE_SEL32_LH(vB3, vB3);

            vC0 = AE_ADD32S(vB0, vB1);
            vC2 = AE_SUB32S(vB0, vB1);
            vC3 = AE_SUBADD32S(vB2, vB3);
            vC1 = AE_ADDSUB32S(vB2, vB3);

            t1_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
            t2_f16x4 = AE_MOVF16X4_FROMINT32X2(t2);
            t3_f16x4 = AE_MOVF16X4_FROMINT32X2(t3);

            vC1 = AE_MOVINT32X2_FROMF32X2(AE_MULFC32X16RAS_L(AE_MOVF32X2_FROMINT32X2(vC1), t1_f16x4)); 
            vC2 = AE_MOVINT32X2_FROMF32X2(AE_MULFC32X16RAS_L(AE_MOVF32X2_FROMINT32X2(vC2), t2_f16x4)); 
            vC3 = AE_MOVINT32X2_FROMF32X2(AE_MULFC32X16RAS_L(AE_MOVF32X2_FROMINT32X2(vC3), t3_f16x4)); 

            AE_S32X2_X (vC3, py1, 3*8*stride);
            AE_S32X2_X (vC2, py1, 1*8*stride);
            AE_S32X2_X (vC1, py1, 2*8*stride);
            AE_S32X2_XC(vC0, py1, 2*4*8*stride);
        }
        stride>>=2;  
        tw_step<<=2;
        acc_inc>>=2; 
    }
    shift += fft_stage_last_ie((int32_t*)x, (int32_t*)y, N); 
    return shift;
} /* fft_cplx32x16_ie() */
#endif
