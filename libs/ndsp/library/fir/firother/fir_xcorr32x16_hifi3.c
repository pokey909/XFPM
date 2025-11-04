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
    Real data circular cross-correlation, 32x16-bit
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Circular Correlation
  Estimates the circular cross-correlation between vectors x (of length N) 
  and y (of length M)  resulting in vector r of length N. It is a similar 
  to correlation but x is read in opposite direction.

  Precision: 
  16x16     16x16-bit data, 16-bit outputs
  24x24     24x24-bit data, 24-bit outputs
  32x16     32x16-bit data, 32-bit outputs
  32x32     32x32-bit data, 32-bit outputs
  f         floating point 


  Input:
  x[N]      input data Q15, Q31 or floating point
  y[M]      input data Q15, Q31 or floating point
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

void fir_xcorr32x16( int32_t * restrict r,
               const int32_t * restrict x,
               const int16_t * restrict y,
               int N, int M )
{
  //
  // Circular cross-correlation algorithm:
  //
  //   r[n] = sum( x[mod(n+m,N)]*y[m] )
  //        m=0..M-1
  //
  //   where n = 0..N-1
  //

  const ae_int32x2 *          X;
  const ae_int32x2 *          S;
  const ae_f16x4   *          Y;
        ae_int32x2 * restrict R;

  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_int32x2 t0, t1, t2, t3, t4, t5;
  ae_f32x2   x0, x1, x2, x3, x4, x5;
  ae_f16x4   y0;

  int n, m;

  NASSERT(r);
  NASSERT(x);
  NASSERT(y);
  NASSERT_ALIGN(r, 8);
  NASSERT_ALIGN(x, 8);
  NASSERT_ALIGN(y, 8);
  NASSERT((M>0) && ((M%4)==0));
  NASSERT((N>0) && ((N%4)==0));

  //
  // Setup pointers and circular addressing for array x[N].
  //

  X = (const ae_int32x2 *)x;
  R = (      ae_int32x2 *)r;

  WUR_AE_CBEGIN0( (uintptr_t)( x + 0 ) );
  WUR_AE_CEND0  ( (uintptr_t)( x + N ) );

  //
  // Compute (N&~7) correlation samples, 8 at a time.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 8 left-most unprocessed x[] entries.
    // Q31
    AE_L32X2_XC( t0, X, +8 );
    AE_L32X2_XC( t1, X, +8 );
    AE_L32X2_XC( t2, X, +8 );
    AE_L32X2_XC( t3, X, +8 );

    x0 = ( t0 );
    x1 = ( t1 );
    x2 = ( t2 );
    x3 = ( t3 );

    // Use the shuttle pointer for the inner loop; preserve X for the next
    // iteration.
    S = X;

    //
    // Inner loop prologue: process first 4 entries of y[M] for 8 accumulators.
    //

    // Load next 4 entries of x[]. Altogether we have 12 x[] entries residing in
    // 6 AE registers.
    // Q31
    AE_L32X2_XC( t4, S, +8 );
    AE_L32X2_XC( t5, S, +8 );

    x4 = ( t4 );
    x5 = ( t5 );

    Y = (const ae_f16x4 *)y;

    // Load y[0]..y[3].
    // Q15
    ae_f16x4_loadip( y0, Y, +8 );

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH( q0, q1, x0, x1, y0 );
    AE_MULFD32X16X2_FIR_HH( q2, q3, x1, x2, y0 );
    AE_MULFD32X16X2_FIR_HH( q4, q5, x2, x3, y0 );
    AE_MULFD32X16X2_FIR_HH( q6, q7, x3, x4, y0 );

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL( q0, q1, x1, x2, y0 );
    AE_MULAFD32X16X2_FIR_HL( q2, q3, x2, x3, y0 );
    AE_MULAFD32X16X2_FIR_HL( q4, q5, x3, x4, y0 );
    AE_MULAFD32X16X2_FIR_HL( q6, q7, x4, x5, y0 );

    // 4 x[] entries are done. Move them out of registers.
    x0 = x2; x1 = x3; x2 = x4; x3 = x5;

    //
    // Inner loop kernel: process next 4 entries of y[M] for 8 accumulators.
    //

    for ( m=0; m<(M>>2)-1; m++ )
    {
      // Load next 4 entries of x[].
      // Q31
      AE_L32X2_XC( t4, S, +8 );
      AE_L32X2_XC( t5, S, +8 );

      x4 = ( t4 );
      x5 = ( t5 );

      // Load y[4*m+4]..y[4*m+7].
      // Q15
      ae_f16x4_loadip( y0, Y, +8 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH( q0, q1, x0, x1, y0 );
      AE_MULAFD32X16X2_FIR_HH( q2, q3, x1, x2, y0 );
      AE_MULAFD32X16X2_FIR_HH( q4, q5, x2, x3, y0 );
      AE_MULAFD32X16X2_FIR_HH( q6, q7, x3, x4, y0 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, x1, x2, y0 );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, x2, x3, y0 );
      AE_MULAFD32X16X2_FIR_HL( q4, q5, x3, x4, y0 );
      AE_MULAFD32X16X2_FIR_HL( q6, q7, x4, x5, y0 );

      // 4 x[] entries are done. Move them out of registers.
      x0 = x2; x1 = x3; 
      x2 = x4; x3 = x5;
    }


    // Store r[8*n+0]..r[8*n+7].
    // Q31 <- Q16.47 - 16 w/ rounding and saturation
    AE_S32X2RA64S_IP(q0, q1, R);
    AE_S32X2RA64S_IP(q2, q3, R);
    AE_S32X2RA64S_IP(q4, q5, R);
    AE_S32X2RA64S_IP(q6, q7, R);
  }

  //
  // If N is not a multiple of 8, compute the last 4 correlation samples.
  //

  if ( N & 4 )
  {
    // Load 4 left-most unprocessed x[] entries.
    // Q31
    AE_L32X2_XC(t0, X, +8);
    AE_L32X2_XC(t1, X, +8);

    x0 = t0;
    x1 = t1;

    // Use the shuttle pointer for the inner loop; preserve X for the next
    // iteration.
    S = X;

    //
    // Inner loop prologue: process first 4 entries of y[M] for 4 accumulators.
    //

    // Load next 4 entries of x[]. Altogether we have 8 x[] entries residing in
    // 4 AE registers.
    // Q31
    AE_L32X2_XC(t2, S, +8);
    AE_L32X2_XC(t3, S, +8);

    x2 = t2;
    x3 = t3;

    Y = (const ae_f16x4 *)y;

    // Load y[0]..y[3].
    // Q15
    ae_f16x4_loadip(y0, Y, +8);

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH(q0, q1, x0, x1, y0);
    AE_MULFD32X16X2_FIR_HH(q2, q3, x1, x2, y0);

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL(q0, q1, x1, x2, y0);
    AE_MULAFD32X16X2_FIR_HL(q2, q3, x2, x3, y0);

    // 4 x[] entries are done. Move them out of registers.
    x0 = x2;
    x1 = x3;

    for ( m=0; m<(M>>2)-1; m++ )
    {
      // Load next 4 entries of x[].
      // Q31
      AE_L32X2_XC(t2, S, +8);
      AE_L32X2_XC(t3, S, +8);

      x2 = t2;
      x3 = t3;

      // Load y[4*m+4]..y[4*m+7].
      // Q15
      ae_f16x4_loadip(y0, Y, +8);

      // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH(q0, q1, x0, x1, y0);
      AE_MULAFD32X16X2_FIR_HH(q2, q3, x1, x2, y0);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL(q0, q1, x1, x2, y0);
      AE_MULAFD32X16X2_FIR_HL(q2, q3, x2, x3, y0);

      // 4 x[] entries are done. Move them out of registers.
      x0 = x2; x1 = x3;
    }

    // Store r[N-4]..r[N-1].
    // Q31 <- Q16.47 - 16 w/ rounding and saturation
    AE_S32X2RA64S_IP(q0, q1, R);
    AE_S32X2RA64S_IP(q2, q3, R);
  }

} /* fir_xcorr32x16() */
