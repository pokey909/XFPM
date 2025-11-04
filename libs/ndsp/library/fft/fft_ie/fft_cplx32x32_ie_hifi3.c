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
#include "NatureDSP_Signal_fft.h"
#include "NatureDSP_Signal_vector.h"
#include "common.h"
#include "fft_twiddles32x32.h"

#define SWAP_PTR(_x, _y) {complex_fract32 *tmp = _x; _x = _y ; _y = tmp; } 

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

int fft_cplx32x32_ie(complex_fract32* y, complex_fract32* x, const complex_fract32* twd, int twdstep, int N, int scalingOpt)
{
    int bexp = 0;
    int v = 1;
    int shift = 0;

    complex_fract32 *pdest = y; 
    int log2N = 30 - NSA(N); 

    const fft_cplx32x32_stage_t first_stg_fn = (scalingOpt == 2) ? fft_stageS2_DFT4_first_32x32 : fft_stageS3_DFT4_first_32x32;
    const fft_cplx32x32_stage_t stg_fn = (scalingOpt == 2) ? fft_stageS2_DFT4_32x32 : fft_stageS3_DFT4_32x32;
    const fft_cplx32x32_stage_t last_stg_fn = 
                                            (log2N & 1) ? 
                                            ((scalingOpt == 2) ? fft_stageS2_DFT2_last_32x32 : fft_stageS3_DFT2_last_32x32) :
                                            ((scalingOpt == 2) ? fft_stageS2_DFT4_last_32x32 : fft_stageS3_DFT4_last_32x32) ;
    NASSERT_ALIGN8(x); 
    NASSERT_ALIGN8(y);
    NASSERT_ALIGN8(twd);
    NASSERT( (N&(N-1))==0 );
    NASSERT(x != y);
    NASSERT(scalingOpt == 2 || scalingOpt == 3); 

    if (scalingOpt==2)
    {
        bexp = vec_bexp32((int32_t*)x, 2*N); 
    }

    shift += first_stg_fn((const int32_t*)twd, (int32_t*)x, (int32_t*)y, N, &v, twdstep, &bexp);
    SWAP_PTR(x, y);
    log2N -= 2;
    twdstep *= 4;

    while ( log2N > 2 )
    {
        shift += stg_fn((const int32_t*)twd, (int32_t*)x, (int32_t*)y, N, &v, twdstep, &bexp);
        SWAP_PTR(x, y); 
        log2N -= 2; 
        twdstep *= 4;
    }
   
    if (y != pdest)
    {  
        /* Execute the last stage inplace */
        y = x;
    }

    /* Last stage */
    shift += last_stg_fn(NULL, (int32_t*)x, (int32_t*)y, N, &v, 0, &bexp);
    return shift;
}




