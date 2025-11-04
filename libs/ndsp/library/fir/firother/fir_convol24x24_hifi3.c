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
    Real data circular convolution, 24x24-bit
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Circular Convolution
  Performs circular convolution between vectors x (of length N) and y (of 
  length M)  resulting in vector r of length N.

  Precision: 
  16x16     16x16-bit data, 16-bit outputs
  24x24     24x24-bit data, 24-bit outputs
  32x16     32x16-bit data, 32-bit outputs 
  32x32     32x32-bit data, 32-bit outputs
  f         floating point

  Input:
  x[N]      input data, Q15, Q31 or floating point
  y[M]      input data, Q15, Q31 or floating point
  N         length of x
  M         length of y
  Output:
  r[N]      output data, Q15, Q31 or floating point

  Restriction:
  x,y,r     should not overlap
  x,y,r     aligned on an 8-bytes boundary
  N,M       multiples of 4 and >0
-------------------------------------------------------------------------*/
/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

void fir_convol24x24( int32_t * restrict r,
                const f24     * restrict x,
                const f24     * restrict y,
                int N, int M )
{
  //
  // Circular convolution algorithm:
  //
  //   r[n] = sum( x[mod(n-m,N)]*y[m] )
  //        m=0..M-1
  //
  //   where n = 0..N-1
  //

  const int32_t  *          xn;
  const ae_f24x2 *          X;
  const ae_f24x2 *          Y;
        ae_f24x2 * restrict R;

  ae_f64   q0, q1, q2, q3, q4, q5, q6, q7;
  ae_f24x2 x0, x1, x2, x3, x4, x5;
  ae_f24x2 y0, y1;

  int n, m;

  ASSERT( r && x && y && (N > 0) && (M > 0) );

  ASSERT( !( N & 3 ) && !( M & 3 ) );

  ASSERT( IS_ALIGN( x ) && IS_ALIGN( y ) && IS_ALIGN( r ) );

  //
  // Setup pointers and circular addressing for array x[N].
  //

  xn = x + 6;

  R = (ae_f24x2 *)r;

  WUR_AE_CBEGIN0( (uintptr_t)( x + 0 ) );
  WUR_AE_CEND0  ( (uintptr_t)( x + N ) );

  //
  // Compute (N&~7) convolution samples, 8 at a time.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Set X to &x[8*n+6]
    X = (const ae_f24x2 *)xn;

    xn += 8;

    // Load x[8*n+7]..x[8*n+0]
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_RIC( x5, X );
    AE_L32X2F24_RIC( x4, X );
    AE_L32X2F24_RIC( x3, X );
    AE_L32X2F24_RIC( x2, X );

    //
    // Inner loop prologue: process first 4 entries of y[M] for 8 accumulators.
    //

    // Load x[8*n-1]..x[8*n-4].
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_RIC( x1, X );
    AE_L32X2F24_RIC( x0, X );

    // Set Y to &y[0]
    Y = (const ae_f24x2 *)y;

    // Load y[0]..y[3].
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( y0, Y, +8 );
    AE_L32X2F24_IP( y1, Y, +8 );

    // Q16.47 <- [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
    AE_MULFD24X2_FIR_H( q1, q0, x2, x1, y0 );
    AE_MULFD24X2_FIR_H( q3, q2, x3, x2, y0 );
    AE_MULFD24X2_FIR_H( q5, q4, x4, x3, y0 );
    AE_MULFD24X2_FIR_H( q7, q6, x5, x4, y0 );

    // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
    AE_MULAFD24X2_FIR_H( q1, q0, x1, x0, y1 );
    AE_MULAFD24X2_FIR_H( q3, q2, x2, x1, y1 );
    AE_MULAFD24X2_FIR_H( q5, q4, x3, x2, y1 );
    AE_MULAFD24X2_FIR_H( q7, q6, x4, x3, y1 );

    // 4 convolution products are done, move out x[8*n+7]..x[8*n+4].
    x5 = x3; x4 = x2; 
    x3 = x1; x2 = x0;

    //
    // Inner loop kernel: process next 4 entries of y[M] for 8 accumulators.
    //

    for ( m=0; m<(M>>2)-1; m++ )
    {
      // Load x[8*n-4*m-5]..x[8*n-4*m-8].
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_RIC( x1, X );
      AE_L32X2F24_RIC( x0, X );

      // Load y[4*m+4]..y[4*m+7].
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( y0, Y, +8 );
      AE_L32X2F24_IP( y1, Y, +8 );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q1, q0, x2, x1, y0 );
      AE_MULAFD24X2_FIR_H( q3, q2, x3, x2, y0 );
      AE_MULAFD24X2_FIR_H( q5, q4, x4, x3, y0 );
      AE_MULAFD24X2_FIR_H( q7, q6, x5, x4, y0 );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q1, q0, x1, x0, y1 );
      AE_MULAFD24X2_FIR_H( q3, q2, x2, x1, y1 );
      AE_MULAFD24X2_FIR_H( q5, q4, x3, x2, y1 );
      AE_MULAFD24X2_FIR_H( q7, q6, x4, x3, y1 );

      // 4 convolution products are done, move out x[8*n-4*m+3]..x[8*n-4*m+0].
      x5 = x3; x4 = x2;
      x3 = x1; x2 = x0;
    }

    // Store r[8*n+0]..r[8*n+7].
    // Q(23+8) <- [ Q16.47 - 24 w/ rounding and saturation ] + 8
    AE_S24X2RA64S_IP( q0, q1, R );
    AE_S24X2RA64S_IP( q2, q3, R );
    AE_S24X2RA64S_IP( q4, q5, R );
    AE_S24X2RA64S_IP( q6, q7, R );
  }

  //
  // If N is not a multiple of 8, compute the last 4 convolution samples.
  //

  if ( N & 4 )
  {
    // Set X to &x[(N&~7)+2].
    X = (const ae_f24x2 *)( xn - 4 );

    // Load x[(N&~7)+3]..x[(N&~7)+0].
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_RIC( x3, X );
    AE_L32X2F24_RIC( x2, X );

    //
    // Inner loop prologue: process first 4 entries of y[M] for 4 accumulators.
    //

    // Load x[(N&~7)-1]..x[(N&~7)-4].
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_RIC( x1, X );
    AE_L32X2F24_RIC( x0, X );

    // Set Y to &y[0].
    Y = (const ae_f24x2 *)y;

    // Load y[0]..y[3].
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( y0, Y, +8 );
    AE_L32X2F24_IP( y1, Y, +8 );

    // Q16.47 <- [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
    AE_MULFD24X2_FIR_H( q1, q0, x2, x1, y0 );
    AE_MULFD24X2_FIR_H( q3, q2, x3, x2, y0 );

    // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
    AE_MULAFD24X2_FIR_H( q1, q0, x1, x0, y1 );
    AE_MULAFD24X2_FIR_H( q3, q2, x2, x1, y1 );

    // 4 convolution products are done, move out x[(N&~7)+3]..x[(N&~7)+0].
    x3 = x1; x2 = x0;

    for ( m=0; m<(M>>2)-1; m++ )
    {
      // Load x[(N&~7)-4*m-5]..x[(N&~7)-4*m-8].
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_RIC( x1, X );
      AE_L32X2F24_RIC( x0, X );

      // Load y[4*m+4]..y[4*m+7].
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( y0, Y, +8 );
      AE_L32X2F24_IP( y1, Y, +8 );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q1, q0, x2, x1, y0 );
      AE_MULAFD24X2_FIR_H( q3, q2, x3, x2, y0 );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q1, q0, x1, x0, y1 );
      AE_MULAFD24X2_FIR_H( q3, q2, x2, x1, y1 );

      // 4 products are done, move out x[(N&~7)-4*m-1]..x[(N&~7)-4*m-4].
      x3 = x1; x2 = x0;
    }

    // Store r[N-4]..r[N-1].
    // Q(23+8) <- [ Q16.47 - 24 w/ rounding and saturation ] + 8
    AE_S24X2RA64S_IP( q0, q1, R );
    AE_S24X2RA64S_IP( q2, q3, R );
  }

} /* fir_convol24x24() */
