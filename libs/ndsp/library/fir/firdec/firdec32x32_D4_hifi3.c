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
  Decimating block real FIR filter, 32x32-bit
*/

/*-----------------------------------------------------------------------------
* Data processing function of a particular decimating filter. Stores a
* block of input samples to the circular delay line buffer and computes
* decimating FIR filter's response.
* Input:
*   delayLine - circular delay line buffer start address
*   delayLen  - Delay line buffer length
*   wrIx    - next position in the buffer to be filled with an input sample
*   x[N*D]  - input samples
*   h[]     - decimating FIR filter coefficients, array layout varies
* Output:
*   y[N]    - output samples
*   retval  - updated index of the oldest sample
* Notes and restrictions:
*   1. Most of data processing functions feature a single, hard-coded
*      decimation factor, so they expect a determined value for parameter D.
*   2. All pointers with the exception of y[N] must be aligned on an 8-bytes
*      boundary.
*   3. N - must be a multiple of 4.
*   4. M - must be a multiple of 4.
-----------------------------------------------------------------------------*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"

#include "firdec32x32_common.h"

#define sz_i32    sizeof(int32_t)

 int fir32x32_D4_proc( int32_t * restrict y,
                                  int32_t * delayLine,
                                  int32_t delayLen,
                          const int32_t * restrict x,
                          const int32_t * restrict h,
                          int wrIx, int D, int N, int M )
{
        ae_int32x2 * restrict D_wr;
  const ae_int32x2 *          D_rd0;
  const ae_int32x2 *          X;
        ae_f32     * restrict Y;
  const ae_int32x2 *          C;

  ae_int32x2 t0, t1, t2, t3, t4, t5, t6, t7;
  ae_int32x2 c0;
  ae_f64   q0, q1, q2, q3;

  int m, n;

  NASSERT( y && delayLine && x && h && (D > 0) && (N > 0) && (M > 0) );
  NASSERT( !( N & 7 ) && !( M & 3 ) );
  NASSERT( IS_ALIGN(delayLine) && IS_ALIGN(x) && IS_ALIGN(h) );

  //
  // Setup pointers and circular delay line buffer.
  //

  D_wr = (      ae_int32x2 *)(delayLine+wrIx);
  X    = (const ae_int32x2 *)x;
  Y    = (      ae_f32     *)y;

  WUR_AE_CBEGIN0( (uintptr_t)(delayLine            ) );
  WUR_AE_CEND0  ( (uintptr_t)(delayLine + delayLen ) );

  //
  // Break the input signal into 4*D-samples blocks. For each block, store
  // 4*D samples to the delay line buffer, and compute 4 samples of decimated
  // response signal.
  //
  __Pragma("loop_count min=1")
  for ( n=0; n<(N>>2); n++ )
  {
    // Reset the coefficients pointer. Now it looks at the tap corresponding
    // to the oldest sample in the delay line.
    C = (const ae_int32x2*)h;

    AE_L32X2_IP(t0, X, +8);
    AE_L32X2_IP(t1, X, +8);
    AE_L32X2_IP(t2, X, +8);
    AE_L32X2_IP(t3, X, +8);
    AE_L32X2_IP(t4, X, +8);
    AE_L32X2_IP(t5, X, +8);
    AE_L32X2_IP(t6, X, +8);
    AE_L32X2_IP(t7, X, +8);

    AE_S32X2_XC(t0, D_wr, +8);
    AE_S32X2_XC(t1, D_wr, +8);
    AE_S32X2_XC(t2, D_wr, +8);
    AE_S32X2_XC(t3, D_wr, +8);
    AE_S32X2_XC(t4, D_wr, +8);
    AE_S32X2_XC(t5, D_wr, +8);
    AE_S32X2_XC(t6, D_wr, +8);
    AE_S32X2_XC(t7, D_wr, +8);

    D_rd0 = (const ae_int32x2 *)D_wr;

    // Zero the accumulators.
    q0 = q1 = q2 = q3 = AE_ZERO64();

    AE_L32X2_XC(t0, D_rd0, +8);
    AE_L32X2_XC(t1, D_rd0, +8);
    AE_L32X2_XC(t2, D_rd0, +8);              
    AE_L32X2_XC(t3, D_rd0, +8);
    AE_L32X2_XC(t4, D_rd0, +8);
    AE_L32X2_XC(t5, D_rd0, +8);           
    
    __Pragma("loop_count min=2")
    for (m = 0; m<(M >> 1) + 2; m++)
    {
        // Load the next 2 tap coefficients.
        // Q31
        AE_L32X2_IP(c0, C, +8);

        AE_L32X2_XC(t6, D_rd0, +8);
   
#if (defined(AE_MULAAFD32R_HH_LL))
      AE_MULAAFD32R_HH_LL(q0, t0, c0);

      AE_MULAAFD32R_HH_LL(q1, t2, c0);

      AE_MULAAFD32R_HH_LL(q2, t4, c0);

      AE_MULAAFD32R_HH_LL(q3, t6, c0);
#else
        AE_MULAF32R_HH(q0, t0, c0);
        AE_MULAF32R_LL(q0, t0, c0);

        AE_MULAF32R_HH(q1, t2, c0);
        AE_MULAF32R_LL(q1, t2, c0);

        AE_MULAF32R_HH(q2, t4, c0);
        AE_MULAF32R_LL(q2, t4, c0);

        AE_MULAF32R_HH(q3, t6, c0);
        AE_MULAF32R_LL(q3, t6, c0);
#endif

        t0 = t1; t1 = t2;
        t2 = t3; t3 = t4;
        t4 = t5; t5 = t6;
    }
                                               
    // Convert and store filter outputs.
    // Q31 <- Q16.47 - 16 w/ rounding and saturation.
    AE_S32RA64S_IP( q0, Y, +4 );
    AE_S32RA64S_IP( q1, Y, +4 );
    AE_S32RA64S_IP( q2, Y, +4 );
    AE_S32RA64S_IP( q3, Y, +4 );
  }
  return (int)((int32_t *)D_wr - delayLine) ;
} /* fir32x32_D4_proc() */
