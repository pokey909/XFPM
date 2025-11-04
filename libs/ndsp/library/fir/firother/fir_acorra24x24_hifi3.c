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
    Real data circular auto-correlation, 24x24-bit, no requirements on vectors
    length and alignment.
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Circular Autocorrelation 
  Estimates the auto-correlation of vector x. Returns autocorrelation of 
  length N.
  These functions implement the circular autocorrelation algorithm with no 
  limitations on x vector length and alignment at the cost of increased 
  processing complexity. In addition, this implementation variant requires
  scratch memory area.

  Precision: 
  16x16    16-bit data, 16-bit outputs
  24x24    24-bit data, 24-bit outputs
  32x32    32-bit data, 32-bit outputs
  f        floating point

  Input:
  s[]     scratch area of
            FIR_ACORRA16X16_SCRATCH_SIZE( N ) or
            FIR_ACORRA24X24_SCRATCH_SIZE( N ) or
            FIR_ACORRA32X32_SCRATCH_SIZE( N ) or
            FIR_ACORRAF_SCRATCH_SIZE( N ) bytes
              
  x[N]    input data Q15, Q31 or floating point
  N       length of x
  Output:
  r[N]    output data, Q15, Q31 or floating point

  Restrictions:
  x,r,s   should not overlap
  N       must be non-zero
  s       aligned on an 8-bytes boundary
-------------------------------------------------------------------------*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
      (void*)( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

void fir_acorra24x24( void * restrict s,
                      f24  * restrict r,
                const f24  * restrict x,
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

  void * s_ptr;
  f24  * x_buf;
  f24  * y_buf;

  const ae_f24x2 *          S;
        ae_f24x2 *          D;
        ae_f24   *          D_s;
  const ae_f24x2 *          X;
  const ae_f24x2 *          Y;
        ae_f24x2 * restrict R;
        ae_f24   * restrict R_s;

  ae_valign S_va, D_va, Y_va, R_va;

  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_f24x2   x0, x1, x2, x3, x4, x5;
  ae_f24x2   y0, y1;
  ae_f24x2   r0, r1, r2, r3;
  ae_f24x2   z_f24x2;

  int M;
  int n, m;

  NASSERT(s);
  NASSERT(r);
  NASSERT(x);
  NASSERT_ALIGN(s,8);
  NASSERT(N>0);

  //
  // Partition the scratch memory.
  //

  s_ptr = s;
  x_buf = (f24*)ALIGNED_ADDR( s_ptr, 8 );
  s_ptr = x_buf + 2*N + 3;

  ASSERT( (int8_t *)s_ptr - (int8_t *)s <= (int)FIR_ACORRA24X24_SCRATCH_SIZE( N ) );

  y_buf = x_buf + N;

  //----------------------------------------------------------------------------
  // Copy x[N] into the scratch are in a way that simplifies the 
  // autocorrelation computation:
  //   x[0]..x[N-1] x[0]..x[N-1] 0 0 0

  //
  // First copy of x[N].
  //

  S = (const ae_f24x2 *)x;
  D = (      ae_f24x2 *)x_buf;

  S_va = AE_LA64_PP( S );

#ifdef COMPILER_XTENSA
  #pragma concurrent
#endif
  for ( n=0; n<((N+1)>>1); n++ )
  {
    // Q23 <- Q(23+8) - 8
    AE_LA32X2F24_IP( x0, S_va, S );
    // Q(23+8) <- Q23 + 8
    AE_S32X2F24_IP( x0, D, +8 );
  }

  //
  // Second copy of x[N].
  //

  S = (const ae_f24x2 *)x;
  D = (      ae_f24x2 *)( (f24 *)D - ( N & 1 ) );

  S_va = AE_LA64_PP( S );
  D_va = AE_ZALIGN64();

#ifdef COMPILER_XTENSA
  #pragma concurrent
#endif
  for ( n=0; n<((N+1)>>1); n++ )
  {
    // Q23 <- Q(23+8) - 8
    AE_LA32X2F24_IP( x0, S_va, S );
    // Q(23+8) <- Q23 + 8
    AE_SA32X2F24_IP( x0, D_va, D );
  }

  AE_SA64POS_FP( D_va, D );

  //
  // Append 3 zeros to allow block-4 processing for the inner loop of the 
  // autocorrelation algorithm.
  //

  D_s = (ae_f24 *)( (f24 *)D - ( N & 1 ) );

  z_f24x2 = 0;

  for ( n=0; n<3; n++ )
  {
    AE_S32F24_L_IP( z_f24x2, D_s, +4 );
  }

  //----------------------------------------------------------------------------
  // Compute first (N&~7) autocorrelation results.

  //
  // Setup pointers and complement the inner loop trip count to the next 
  // multiple of 4.
  //

  X = (const ae_f24x2 *)x_buf;
  R = (      ae_f24x2 *)r;

  R_va = AE_ZALIGN64();

  M = N + ( -N & 3 );

  //
  // Compute (N&~7) autocorrelation samples, 8 at a time.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 8 left-most unprocessed left-hand array entries.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( x0, X, +8 );
    AE_L32X2F24_IP( x1, X, +8 );
    AE_L32X2F24_IP( x2, X, +8 );
    AE_L32X2F24_IP( x3, X, +8 );

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
    AE_L32X2F24_IP( x4, S, +8 );
    AE_L32X2F24_IP( x5, S, +8 );

    // Reset the right-hand array pointer. The second copy of x[N] that is used
    // for the right-hand array starts at an unaligned address whenever N is odd
    Y = (const ae_f24x2 *)y_buf;

    Y_va = AE_LA64_PP( Y );

    // Load first 4 right-hand array entries.
    // Q23 <- Q(23+8) - 8
    AE_LA32X2F24_IP( y0, Y_va, Y );
    AE_LA32X2F24_IP( y1, Y_va, Y );

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
    x0 = x2; x1 = x3; x2 = x4; x3 = x5;

    //
    // Inner loop kernel: process 4 entries of the right-hand array
    // for 8 accumulators. Registers-based delay line is updated on each 
    // iteration in the same way as for the loop prologue.
    //

    for ( m=0; m<(M>>2)-1; m++ )
    {
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( x4, S, +8 );
      AE_L32X2F24_IP( x5, S, +8 );

      // Q23 <- Q(23+8) - 8
      AE_LA32X2F24_IP( y0, Y_va, Y );
      AE_LA32X2F24_IP( y1, Y_va, Y );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
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
      x0 = x2; x1 = x3; x2 = x4; x3 = x5;
    }

    //
    // Format and save 8 autocorrelation results.
    //

    // Q23 <- Q16.47 - 24 w/ rounding and saturation.
    r0 = AE_ROUND24X2F48SASYM( q0, q1 );
    r1 = AE_ROUND24X2F48SASYM( q2, q3 );
    r2 = AE_ROUND24X2F48SASYM( q4, q5 );
    r3 = AE_ROUND24X2F48SASYM( q6, q7 );

    // Q(23+8) <- Q23 + 8
    AE_SA32X2F24_IP( r0, R_va, R );
    AE_SA32X2F24_IP( r1, R_va, R );
    AE_SA32X2F24_IP( r2, R_va, R );
    AE_SA32X2F24_IP( r3, R_va, R );
  }

  //----------------------------------------------------------------------------
  // (N&~7) autocorrelation results are done by now. If (N&4) != 0, then compute
  // a block of 4 results.

  if ( N & 4 )
  {
    // Load 4 left-most unprocessed left-hand array entries.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( x0, X, +8 );
    AE_L32X2F24_IP( x1, X, +8 );

    // Use the shuttle pointer for the inner loop; preserve X for the next
    // iteration.
    S = X;

    //
    // Inner loop prologue: process first 4 entries of the right-hand array
    // for 4 accumulators.
    //

    // Load 4 more left-hand array entries. Altogether we have 8 entries 
    // residing in 4 AE registers.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( x2, S, +8 );
    AE_L32X2F24_IP( x3, S, +8 );

    // Reset the right-hand array pointer.
    Y = (const ae_f24x2 *)y_buf;

    Y_va = AE_LA64_PP( Y );

    // Load first 4 right-hand array entries.
    // Q23 <- Q(23+8) - 8
    AE_LA32X2F24_IP( y0, Y_va, Y );
    AE_LA32X2F24_IP( y1, Y_va, Y );

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

    for ( m=0; m<(M>>2)-1; m++ )
    {
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( x2, S, +8 );
      AE_L32X2F24_IP( x3, S, +8 );

      // Q23 <- Q(23+8) - 8
      AE_LA32X2F24_IP( y0, Y_va, Y );
      AE_LA32X2F24_IP( y1, Y_va, Y );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
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

    // Q23 <- Q16.47 - 24 w/ rounding and saturation.
    r0 = AE_ROUND24X2F48SASYM( q0, q1 );
    r1 = AE_ROUND24X2F48SASYM( q2, q3 );

    // Store r[N-4]..r[N-1].
    // Q(23+8) <- Q23 + 8
    AE_SA32X2F24_IP( r0, R_va, R );
    AE_SA32X2F24_IP( r1, R_va, R );
  }

  AE_SA64POS_FP( R_va, R );

  //----------------------------------------------------------------------------
  // (N&~3) autocorrelation results are done by now. If N is not a multiple
  // of 4, compute the last (N&3) results.

  if ( N & 3 )
  {
    // Load 4 left-most unprocessed left-hand array entries.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( x0, X, +8 );
    AE_L32X2F24_IP( x1, X, +8 );

    //
    // Inner loop prologue: process first 4 entries of the right-hand array
    // for 4 accumulators.
    //

    // Load 4 more left-hand array entries. Altogether we have 8 entries 
    // residing in 4 AE registers.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( x2, X, +8 );
    AE_L32X2F24_IP( x3, X, +8 );

    // Reset the right-hand array pointer.
    Y = (const ae_f24x2 *)y_buf;

    Y_va = AE_LA64_PP( Y );

    // Load first 4 right-hand array entries.
    // Q23 <- Q(23+8) - 8
    AE_LA32X2F24_IP( y0, Y_va, Y );
    AE_LA32X2F24_IP( y1, Y_va, Y );

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

    for ( m=0; m<(M>>2)-1; m++ )
    {
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( x2, X, +8 );
      AE_L32X2F24_IP( x3, X, +8 );

      // Q23 <- Q(23+8) - 8
      AE_LA32X2F24_IP( y0, Y_va, Y );
      AE_LA32X2F24_IP( y1, Y_va, Y );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, x0, x1, y0 );
      AE_MULAFD24X2_FIR_H( q2, q3, x1, x2, y0 );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, x1, x2, y1 );
      AE_MULAFD24X2_FIR_H( q2, q3, x2, x3, y1 );

      // 4 left-hand array entries are done, shift them out of the registers.
      x0 = x2; x1 = x3;
    }

    //
    // Format and save (N&3) autocorrelation results.
    //

    R_s = (ae_f24 *)( r + N-1 );

    switch ( N & 3 )
    {
    // Q(23+8) <- [ Q16.47 - 24 w/ rounding and saturation ] + 8
    case 3: AE_S24RA64S_IP( q2, R_s, -4 );
    case 2: AE_S24RA64S_IP( q1, R_s, -4 );
    }

    // Q(23+8) <- [ Q16.47 - 24 w/ rounding and saturation ] + 8
    AE_S24RA64S_IP( q0, R_s, -4 );
  }

} /* fir_acorra24x24() */
