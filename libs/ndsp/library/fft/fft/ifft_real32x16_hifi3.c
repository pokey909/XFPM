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
    C code optimized for HiFi3
    Real/halfcomplex FFT 32x16 with scaling option 3
	Integrit, 2006-2017
*/

#include "NatureDSP_Signal_fft.h"
#include "fft_real_twiddles.h"
#include "fft_cplx_twiddles.h"

/*
	in-place inverse split part of FFT:
	x[N+2]  input (N+2 samples)/output(N samples)
	N       size of FFT
*/
static void isplitPart_x2(int32_t *X,int N)
{
  int i, step;

  const int step_back = -8;

  ae_int32x2 * restrict p_x0, * restrict p_x1;
  const ae_p16x2s * restrict p_twd;

  ae_int32x2  vA0, vA1, vB0, vB1, vC0, vC1;
  ae_int16x4  vT;
  ae_f32x2    vF0, vF1;
  ae_f16x4    vF2;

  NASSERT_ALIGN8(X);

  // setup table step
  {
    ae_q56s tmp = AE_CVTQ48A32S(N);
    step = 1<<(AE_NSAQ56S(tmp)-(38-MAX_RFFT_PWR));
  }
  step <<= 2;

  p_twd = (const ae_p16x2s *)twiddleSplit;
  p_x0 = (ae_int32x2 *)X;
  p_x1 = (ae_int32x2 *)(X+N);


  // first point
  vA0 = AE_L32X2_I(p_x0, 0);
  vA1 = AE_L32X2_I(p_x1, 0);
  // Additional scaling
  vA0 = AE_SRAI32(vA0, 1);
  vA1 = AE_SRAI32(vA1, 1);

  vB0 = AE_ADD32S(vA0, vA1);
  vB1 = AE_SUB32S(vA0, vA1);
  vB0 = AE_SEL32_HH(vB0, vB1);

  vB1 = AE_MOVI(0);
  vB0 = AE_SRAI32(vB0, 1);

  AE_S32X2_IP(vB0, p_x0, 8);
  AE_S32X2_XP(vB1, p_x1, step_back);

  vA0 = AE_L32X2_I(p_x0, 0);
  vA1 = AE_L32X2_I(p_x1, 0);
  // Additional scaling
  vA0 = AE_SRAI32(vA0, 1);
  vA1 = AE_SRAI32(vA1, 1);

  for(i = 1; i < (N>>2); i++) 
  {
    // load twiddle
    AE_L16X2M_XU(vB1, p_twd, step);
    vB1 = AE_SRAI32(vB1, 8);
    vT = AE_CVT16X4(vB1, vB1);

    // ADD/SUBB
    vB0 = AE_ADD32S(vA0, vA1);
    vB1 = AE_SUB32S(vA0, vA1);
    vB0 = AE_SRAI32(vB0, 1);
    vB1 = AE_SRAI32(vB1, 1);

    vA0 = AE_SEL32_LH(vB0, vB1);
    vB1 = AE_SEL32_HL(vB0, vB1);

    // do rotation
    vF1 = (vA0);
    vF2 = (vT);
    vF0 = AE_MULFC32X16RAS_H(vF1, vF2);
    vB0 = (vF0);
    vA1 = AE_NEG32S(vB0);
    vB0 = AE_SEL32_LH(vA1, vA1);

    // load next data
    vA0 = AE_L32X2_I(p_x0, 8);
    vA1 = AE_L32X2_I(p_x1, -8);
    // Additional scaling
    vA0 = AE_SRAI32(vA0, 1);
    vA1 = AE_SRAI32(vA1, 1);

    // SUM/DIFF
    vC0 = AE_ADD32S(vB1, vB0);
    vC1 = AE_SUB32S(vB1, vB0);
    vB1 = AE_NEG32S(vC1);
    vC1 = AE_SEL32_HL(vC1, vB1);

    AE_S32X2_IP(vC0, p_x0, 8);
    AE_S32X2_XP(vC1, p_x1, step_back);
  }

  // middle sample
  vB0 = AE_NEG32S(vA0);
  vC0 = AE_SEL32_HL(vA0, vB0);
  AE_S32X2_I(vC0, p_x0, 0);
}

/*-------------------------------------------------------------------------
  Inverse FFT on Real Data
  These functions make inverse FFT on half spectral data forming real
  data samples.
      Scaling: 
      +-------------------+----------------------------------------+
      |      Function     |           Scaling options              |
      +-------------------+----------------------------------------+
      |  ifft_real16x16   |  2 - 16-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  ifft_real32x32   |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  ifft_real32x16   |  3 - fixed scaling before each stage   | 
      |  ifft_real24x24   |  0 - no scaling                        | 
      |                   |  1 - 24-bit scaling                    |
      |                   |  2 - 32-bit scaling on the first stage |
      |                   |  and 24-bit scaling later              |
      |                   |  3 - fixed scaling before each stage   |
      +-------------------+----------------------------------------+

  NOTES:
  1. Bit-reversing reordering is done here. 
  2. IFFT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED after
     the call.
  3. Inverse FFT function for real signal transforms the input spectrum  
     and then calls ifft_cplx() with FFT size set to N/2.
  4. 32x32 FFTs support mixed radix transforms
  5. N - FFT size

  Precision:
  32x32  32-bit input/outputs, 32-bit twiddles
  24x24  24-bit input/outputs, 24-bit twiddles
  32x16  32-bit input/outputs, 16-bit twiddles
  16x16  16-bit input/outputs, 16-bit twiddles

  Input:
  x[(N/2+1)*2]	input spectrum. Real and imaginary data are interleaved  
                and real data goes first
  scalingOpt	scaling option (see table above)

  Output:			
  y[N]	        real output signal

  Returned value: total number of right shifts occurred during scaling 
                  procedure

  Restrictions:
  x,y           should not overlap
  x,y           aligned on a 8-bytes boundary
-------------------------------------------------------------------------*/
int ifft_real32x16( 
              f24* y,
              int32_t* x,
              fft_handle_t h,
              int scalingOpt)
{
    int scale=2; // Scaling of isplitPart_x2
    int N;
    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT(scalingOpt==3);
    N=(((const tFftDescr*)h)->N)<<1;

    isplitPart_x2(x, N);
    scale +=  ifft_cplx32x16(y,x,h,3);
    return scale;
}
