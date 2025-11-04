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
    FFT on Real Data with Optimized Memory Usage
    C code optimized for HiFi3
	Integrit, 2006-2017
*/

#include "NatureDSP_Signal_fft.h"
#include "common.h"
#include "fft_real_twiddles_24x24.h"
#include "fft_cplx_twiddles.h"


/*
    in-place split part of FFT:
    x[N+2]    input/output
    N        size of FFT
*/
static void splitPart_x2_24x24(f24 * x, int N, const f24 *tw, int step)
{
  int i;
  const int step_back = -8;

  ae_f24x2 * restrict p_x0, * restrict p_x1;
  const ae_f24x2 * restrict p_twd = (const ae_f24x2 *)tw;

  ae_int24x2  vA0, vA1, vB0, vB1, vC0, vC1, vR;
  ae_f24x2    vT;
  ae_f32x2    vF2;
  ae_f24x2    vF0, vF1;

  NASSERT_ALIGN8(x);

  step*=sizeof(f24)*2;
  
  p_x0 = (ae_f24x2 *)x;
  p_x1 = (ae_f24x2 *)(x+N);

  // load data and prepare pointers for pre-increment
  // first and last samples
  vF0 = AE_L32X2F24_I(p_x0, 0);
  vA0 = (vF0);
  vA1 = AE_SELP24_LH(vA0, vA0);
  vB0 = AE_ADDSP24S(vA0, vA1);
  vB1 = AE_SUBSP24S(vA0, vA1);
  vA1 = 0;
  vB0 = AE_SELP24_HH(vB0, vA1);
  vB1 = AE_SELP24_HH(vB1, vA1);
  vF0 = (vB0);
  vF1 = (vB1);
  AE_S32X2F24_IP(vF0, p_x0, 8);
  AE_S32X2F24_XP(vF1, p_x1, step_back);

  vF0 = AE_L32X2F24_I(p_x0, 0);
  vF1 = AE_L32X2F24_I(p_x1, 0);
  vA0 = (vF0);
  vA1 = (vF1);

  vR = 1;

  AE_L32X2F24_XP(vT, p_twd, step);

  for (i = 1; i < N>>2; i++) 
  {
    ae_int32x2 tmp; 
    // load twiddle
    AE_L32X2F24_XP(vT, p_twd, step);

    // ADD/SUBB
    vB0 = AE_ADDSP24S(vA0, vA1);
    vB1 = AE_SUBSP24S(vA0, vA1);

    vA0 = AE_SELP24_HL(vB1, vB0);
    vB1 = AE_SELP24_HL(vB0, vB1);

    // do rotation
    vF1 =  (vA0);

    // conj (tw)
    tmp = AE_MOVINT32X2_FROMF24X2(vT); 
    tmp = AE_ADDSUB32S(0, tmp); 
    vT =  AE_MOVF24X2_FROMINT32X2(tmp);
    vT = AE_SELP24_LH(vT, vT); 

    vF2 = AE_MULFC24RA(vF1, vT);
    vB0 = AE_MOVINT24X2_FROMF32X2(vF2);

    // load next data
    vF0 = AE_L32X2F24_I(p_x0, 8);
    vA0 = (vF0);
    vF1 = AE_L32X2F24_X(p_x1, step_back);
    vA1 = (vF1);

    vB0 = AE_ADDSP24S(vB0, vR);
    vB1 = AE_ADDSP24S(vB1, vR);
    vB0 = AE_SRAIP24(vB0, 1);
    vB1 = AE_SRAIP24(vB1, 1);

    // SUM/DIFF
    vC0 = AE_SUBSP24S(vB1, vB0);
    vC1 = AE_ADDSP24S(vB1, vB0);
    vB1 = AE_NEGSP24S(vC1);
    vC1 = AE_SELP24_HL(vC1, vB1);

    vF0 =  (vC0);
    AE_S32X2F24_IP(vF0, p_x0, 8);
    vF1 =  (vC1);
    AE_S32X2F24_XP(vF1, p_x1, step_back);
  }

  // middle sample
  vB0 = AE_NEGSP24S(vA0);
  vC0 = AE_SELP24_HL(vA0, vB0);
  vF0 =  (vC0);
  AE_S32X2F24_I(vF0, p_x0, 0);
} /* splitPart_x2_24x24() */

/*-------------------------------------------------------------------------
  FFT on Real Data with Optimized Memory Usage
  These functions make FFT on real data forming half of spectrum with
  optimized memory usage.
  Scaling: 
      +-----------------------+--------------------------------------+
      |      Function         |           Scaling options            |
      +-----------------------+--------------------------------------+
      |  fft_real16x16_ie     |  2 - 16-bit dynamic scaling          |
      |  fft_real32x16_ie     |  3 - fixed scaling before each stage |
      |  fft_real24x24_ie     |  3 - fixed scaling before each stage |
      |  fft_real24x24_ie_24p |  3 - fixed scaling before each stage |
      |  fft_real32x16_ie_24p |  1 - 24-bit scaling                  |
      |  fft_real32x32_ie     |  2 - 32-bit dynamic scaling          |    
      |                       |  3 - fixed scaling before each stage |   
      +-----------------------+--------------------------------------+
    
  NOTES:
  1. Bit-reversing reordering is done here.
  2. INPUT DATA MAY APPEAR DAMAGED after the call.
  3. FFT functions may use input and output buffers for temporal storage
     of intermediate 32-bit data, so FFT functions with 24-bit packed
     I/O (Nx3-byte data) require that the buffers are large enough to 
     keep Nx4-byte data.
  4. FFT of size N may be supplied with constant data (twiddle factors) 
     of a larger-sized FFT = N*twdstep.

  Precision:
  16x16_ie      16-bit input/outputs, 16-bit data, 16-bit twiddles
  24x24_ie      24-bit input/outputs, 24-bit data, 24-bit twiddles
  32x16_ie      32-bit input/outputs, 32-bit data, 16-bit twiddles
  32x32_ie      32-bit input/outputs, 32-bit data, 32-bit twiddles
  24x24_ie_24p  24-bit packed input/outputs, 24-bit data, 24-bit twiddles
  32x16_ie_24p  24-bit packed input/outputs, 32-bit data, 16-bit twiddles
  f_ie          floating point

  Input:
  x                     input signal: 
  --------------+----------+-----------------+----------+
  Function      |   Size   |  Allocated Size |  type    |
  --------------+----------+-----------------+-----------
  16x16_ie      |     N    |      N          |  int16_t |
  24x24_ie      |     N    |      N          |   f24    |
  32x16_ie      |     N    |      N          |  int32_t |
  32x32_ie      |     N    |      N          |  int32_t |
  24x24_ie_24p  |     3*N  |      4*N+8      |  uint8_t |
  32x16_ie_24p  |     3*N  |      4*N+8      |  uint8_t |
  --------------+----------+-----------------+----------+

  twd[N*twdstep*3/4]    twiddle factor table of a complex-valued 
                        FFT of size N*twdstep
  N                     FFT size
  twdstep               twiddle step
  scalingOpt            scaling option (see table above), not applicable 
                        to the floating point function

  Output:
  y                     output spectrum. Real and imaginary data 
                        are interleaved and real data goes first:
  --------------+----------+-----------------+---------------+
  Function      |   Size   |  Allocated Size |  type         |
  --------------+----------+-----------------+----------------
  16x16_ie      |   N/2+1  |      N/2+1      |complex_fract16|
  24x24_ie      |   N/2+1  |      N/2+1      |complex_fract32|
  32x16_ie      |   N/2+1  |      N/2+1      |complex_fract32|
  32x32_ie      |   N/2+1  |      N/2+1      |complex_fract32|
  24x24_ie_24p  |  3*(N+2) |      4*N+8      |  uint8_t      |
  32x16_ie_24p  |  3*(N+2) |      4*N+8      |  uint8_t      |
  f_ie          |   N/2+1  |      N/2+1      |complex_float  |
  --------------+----------+-----------------+---------------+

  Returned value: total number of right shifts occurred during scaling
  procedure

  Restrictions:
  x,y     should not overlap
  x,y     aligned on 8-bytes boundary
-------------------------------------------------------------------------*/
 
int fft_real24x24_ie    (complex_fract32* y, f24      * x, const complex_fract32* twd, int twdstep, int N, int scalingOpt)
{
    int scale = 0;
    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT(scalingOpt==3);
    NASSERT(N==256||N==512||N==1024);

    scale = fft_cplx24x24_ie(y, (complex_fract32*)x, twd, twdstep*2, N/2, scalingOpt);  
    splitPart_x2_24x24((int32_t*)y, N, (const int32_t*)twd, twdstep);

  return scale;
} /* fft_real24x24_ie() */
