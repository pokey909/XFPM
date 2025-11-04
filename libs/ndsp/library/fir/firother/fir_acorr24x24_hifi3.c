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
    Real data circular auto-correlation, 24x24-bit
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Circular Autocorrelation 
  Estimates the auto-correlation of vector x. Returns autocorrelation of 
  length N.

  Precision: 
  16x16     16-bit data, 16-bit outputs
  24x24     24-bit data, 24-bit outputs
  32x32     32-bit data, 32-bit outputs
  f         floating point

  Input:
  x[N]      input data, Q15, Q31 or floating point
  N         length of x
  Output:
  r[N]      output data, Q15, Q31 or floating point

  Restrictions:
  x,r       should not overlap
  x,r       aligned on an 8-bytes boundary
  N         multiple of 4 and >0
-------------------------------------------------------------------------*/
/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

void fir_acorr24x24( f24 * restrict r,
               const f24 * restrict x,
               int N )
{
  //
  // Circular autocorrelation algorithm:
  //
  //   r[n] = sum( x[mod(n+m,N)]*x[m] )
  //        m=0..N-1
  //
  // where n = 0..N-1. Throughout the algorithm implementation, x[mod(n+m,N)],
  // m = 0..N-1 is referred to as the left-hand array, while x[m] is called
  // the right-hand array, and it is accessed through the Y pointer.
  //

  const ae_f24x2 *          X;
  const ae_f24x2 *          S;
  const ae_f24x2 *          Y;
        ae_f24x2 * restrict R;

  ae_f64   q0, q1, q2, q3, q4, q5, q6, q7;
  ae_f24x2 x0, x1, x2, x3, x4, x5;
  ae_f24x2 y0, y1;

  int n, m;

  NASSERT(x);
  NASSERT(r);
  NASSERT_ALIGN(x,8);
  NASSERT_ALIGN(r,8);
  NASSERT((N>0) && ((N%4)==0));

  //
  // Setup pointers and circular addressing for array x[].
  //

  X = (const ae_f24x2 *)x;
  R = (      ae_f24x2 *)r;

  WUR_AE_CBEGIN0( (uintptr_t)( x + 0 ) );
  WUR_AE_CEND0  ( (uintptr_t)( x + N ) );

  //
  // Compute (N&~7) autocorrelation samples, 8 at a time.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 8 left-most unprocessed left-hand array entries.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_XC( x0, X, +8 );
    AE_L32X2F24_XC( x1, X, +8 );
    AE_L32X2F24_XC( x2, X, +8 );
    AE_L32X2F24_XC( x3, X, +8 );

    // Use the shuttle pointer for the inner loop; preserve X for the next
    // iteration.
    S = X;

    //
    // Inner loop prologue: process first 4 entries of the right-hand array
    // for 8 accumulators.
    //

    // Load 4 more left-hand array entries. Altogether we have 12 entries 
    // residing in 6 AE registers.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_XC( x4, S, +8 );
    AE_L32X2F24_XC( x5, S, +8 );

    // Reset the right-hand array pointer.
    Y = (const ae_f24x2 *)x;

    // Load first 4 right-hand array entries.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( y0, Y, +8 );
    AE_L32X2F24_IP( y1, Y, +8 );

    // Q16.47 <- [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
    AE_MULFD24X2_FIR_H( q0, q1, x0, x1, y0 );
    AE_MULFD24X2_FIR_H( q2, q3, x1, x2, y0 );
    AE_MULFD24X2_FIR_H( q4, q5, x2, x3, y0 );
    AE_MULFD24X2_FIR_H( q6, q7, x3, x4, y0 );

    // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
    AE_MULAFD24X2_FIR_H( q0, q1, x1, x2, y1 );
    AE_MULAFD24X2_FIR_H( q2, q3, x2, x3, y1 );
    AE_MULAFD24X2_FIR_H( q4, q5, x3, x4, y1 );
    AE_MULAFD24X2_FIR_H( q6, q7, x4, x5, y1 );

    // 4 left-hand array entries are done, shift them out of the registers.
    x0 = x2; x1 = x3; 
    x2 = x4; x3 = x5;

    //
    // Inner loop kernel: process 4 entries of the right-hand array
    // for 8 accumulators. Registers-based delay line is updated on each 
    // iteration in the same way as for the loop prologue.
    //

    for ( m=0; m<(N>>2)-1; m++ )
    {
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_XC( x4, S, +8 );
      AE_L32X2F24_XC( x5, S, +8 );

      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( y0, Y, +8 );
      AE_L32X2F24_IP( y1, Y, +8 );

      // Q16.47 <- [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, x0, x1, y0 );
      AE_MULAFD24X2_FIR_H( q2, q3, x1, x2, y0 );
      AE_MULAFD24X2_FIR_H( q4, q5, x2, x3, y0 );
      AE_MULAFD24X2_FIR_H( q6, q7, x3, x4, y0 );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, x1, x2, y1 );
      AE_MULAFD24X2_FIR_H( q2, q3, x2, x3, y1 );
      AE_MULAFD24X2_FIR_H( q4, q5, x3, x4, y1 );
      AE_MULAFD24X2_FIR_H( q6, q7, x4, x5, y1 );

      // 4 left-hand array entries are done, shift them out of the registers.
      x0 = x2; x1 = x3; 
      x2 = x4; x3 = x5;
    }

    //
    // Format and save 8 autocorrelation results.
    //

    // Store r[8*n+0]..r[8*n+7].
    // Q(23+8) <- [ Q16.47 - 24 w/ rounding and saturation ] + 8
    AE_S24X2RA64S_IP( q0, q1, R );
    AE_S24X2RA64S_IP( q2, q3, R );
    AE_S24X2RA64S_IP( q4, q5, R );
    AE_S24X2RA64S_IP( q6, q7, R );
  }

  //
  // If N is not a multiple of 8, compute the last 4 autocorrelation results.
  //

  if ( N & 4 )
  {
    // Load 4 left-most unprocessed left-hand array entries.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_XC( x0, X, +8 );
    AE_L32X2F24_XC( x1, X, +8 );

    //
    // Inner loop prologue: process first 4 entries of the right-hand array
    // for 4 accumulators.
    //

    // Load 4 more left-hand array entries. Altogether we have 8 entries 
    // residing in 4 AE registers.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_XC( x2, X, +8 );
    AE_L32X2F24_XC( x3, X, +8 );

    // Reset the right-hand array pointer.
    Y = (const ae_f24x2 *)x;

    // Load first 4 right-hand array entries.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( y0, Y, +8 );
    AE_L32X2F24_IP( y1, Y, +8 );

    // Q16.47 <- [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
    AE_MULFD24X2_FIR_H( q0, q1, x0, x1, y0 );
    AE_MULFD24X2_FIR_H( q2, q3, x1, x2, y0 );

    // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
    AE_MULAFD24X2_FIR_H( q0, q1, x1, x2, y1 );
    AE_MULAFD24X2_FIR_H( q2, q3, x2, x3, y1 );

    // 4 left-hand array entries are done, shift them out of the registers.
    x0 = x2; x1 = x3;

    //
    // Inner loop kernel: process 4 entries of the right-hand array
    // for 4 accumulators. Registers-based delay line is updated on each 
    // iteration in the same way as for the loop prologue.
    //

    for ( m=0; m<(N>>2)-1; m++ )
    {
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_XC( x2, X, +8 );
      AE_L32X2F24_XC( x3, X, +8 );

      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( y0, Y, +8 );
      AE_L32X2F24_IP( y1, Y, +8 );

      // Q16.47 <- [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, x0, x1, y0 );
      AE_MULAFD24X2_FIR_H( q2, q3, x1, x2, y0 );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, x1, x2, y1 );
      AE_MULAFD24X2_FIR_H( q2, q3, x2, x3, y1 );

      // 4 left-hand array entries are done, shift them out of the registers.
      x0 = x2; x1 = x3;
    }

    //
    // Format and save 4 autocorrelation results.
    //

     // Store r[N-4]..r[N-1].
    // Q(23+8) <- [ Q16.47 - 24 w/ rounding and saturation ] + 8
    AE_S24X2RA64S_IP( q0, q1, R );
    AE_S24X2RA64S_IP( q2, q3, R );
  }

} /* fir_acorr24x24() */
