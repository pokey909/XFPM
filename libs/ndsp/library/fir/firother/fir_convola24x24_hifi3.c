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
    Real data circular convolution, 24x24-bit, no requirements on vectors 
    length and alignment.
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Circular Convolution
  Performs circular convolution between vectors x (of length N) and y (of 
  length M) resulting in vector r of length N.
  These functions implement the circular convolution algorithm with no 
  limitations on x and y vectors length and alignment at the cost of 
  increased processing complexity. In addition, this implementation variant
  requires scratch memory area.

  Precision: 
  16x16     16x16-bit data, 16-bit outputs
  24x24     24x24-bit data, 24-bit outputs
  32x16     32x16-bit data, 32-bit outputs 
  32x32     32x32-bit data, 32-bit outputs
  f         floating point

  Input:
  s[]       scratch area, 
              FIR_CONVOLA16X16_SCRATCH_SIZE(N,M) or
              FIR_CONVOLA24X24_SCRATCH_SIZE(N,M) or
              FIR_CONVOLA32X16_SCRATCH_SIZE(N,M) or
              CXFIR_CONVOLA32X16_SCRATCH_SIZE(N,M) or
              FIR_CONVOLA32X32_SCRATCH_SIZE(N,M) or
              FIR_CONVOLAF_SCRATCH_SIZE(N,M) bytes

  x[N]      input data Q15, Q31 or floating point
  y[M]      input data Q15, Q31 or floating point
  N         length of x
  M         length of y
  Output:
  r[N]      output data, Q15, Q31 or floating point

  Restriction:
  x,y,r,s   should not overlap
  s         must be aligned on an 8-bytes boundary
  N,M       must be >0
  N >= M-1  minimum allowed length of vector x is the length of y minus one
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

void fir_convola24x24( void * restrict s,
                       f24  * restrict r,
                 const f24  * restrict x,
                 const f24  * restrict y,
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

  void * s_ptr;
  f24  * x_buf;
  int32_t * x_tmp;
  f24  * y_buf;

  const ae_f24x2 *          S;
        ae_f24x2 * restrict D;
        ae_f24   * restrict D_s;
  const ae_f24x2 *          X;
  const ae_f24x2 *          Y;
        ae_f24x2 * restrict R;
        ae_f24   * restrict R_s;

  ae_valign S_va, D_va, R_va;

  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_f24x2   x0, x1, x2, x3, x4, x5;
  ae_f24x2   y0, y1;
  ae_f24x2   r0, r1, r2, r3;
  ae_int24x2 z_24x2;
  ae_f24x2   z_f24x2;

  int _M;
  int n, m;

  ASSERT( s && r && x && y && (N > 0) && (M > 0) && (N >= (M-1)) );

  ASSERT( IS_ALIGN( s ) );

  //----------------------------------------------------------------------------
  // Partition the scratch memory area.

  s_ptr = s;
  x_buf = (f24*)ALIGNED_ADDR( s_ptr, 8 );
  s_ptr = x_buf + M-1 + N + 3;
  y_buf = (f24*)ALIGNED_ADDR( s_ptr, 8 );
  s_ptr = y_buf + M + 3;

  ASSERT( (int8_t *)s_ptr - (int8_t *)s <= (int)FIR_CONVOLA24X24_SCRATCH_SIZE(N, M) );
  x_tmp = x_buf + M-1 + N + 3 - 4;
  x_tmp[0]=0;
  x_tmp[1]=0;
  x_tmp[2]=0;
  x_tmp[3]=0;
  x_tmp[4]=0;
  x_tmp[5]=0;
  x_tmp[6]=0;
  x_tmp[7]=0;

  //----------------------------------------------------------------------------
  // Copy x[N] data into the aligned buffer in a way that simplifies the
  // convolution calculation:
  //  x_buf[M-1+N+3] = { x[N-(M-1)]..x[N-1] x[0] x[1]..x[N-1] X X X }
  // Three X locations are reserved to allow for block-4 processing.

  //
  // Copy last M-1 entries of x{N].
  //

  S = (const ae_f24x2 *)( x + N - (M-1) );
  D = (      ae_f24x2 *)x_buf;

  S_va = AE_LA64_PP( S );

  for ( m=0; m<(M>>1); m++ )
  {
    // Q23 <- Q(23+8) - 8
    AE_LA32X2F24_IP( x0, S_va, S );
    // Q(23+8) <- Q23 + 8
    AE_S32X2F24_IP( x0, D, +8 );
  }

  //
  // Copy x[N].
  //

  S = (const ae_f24x2 *)x;
  D = (      ae_f24x2 *)( (f24 *)D - !( M & 1 ) );

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

  //----------------------------------------------------------------------------
  // Copy reverted y[M] data into the aligned buffer and append 3 zeros:
  //  y_buf[M+3] = { y[M-1]..y[0] 0 0 0 }

  //
  // Copy y[M] in reverted order.
  //

  S = (const ae_f24x2 *)( y + M-1 );
  D = (      ae_f24x2 *)y_buf;

  S_va = AE_LA64_PP( S );

  for ( m=0; m<((M+1)>>1); m++ )
  {
    // Q23 <- Q(23+8) - 8
    AE_LA32X2F24_RIP( y0, S_va, S );
    // Q(23+8) <- Q23 + 8
    AE_S32X2F24_IP( y0, D, +8 );
  }

  //
  // Append three zeros to allow for block-4 processing.
  //

  D_s = (ae_f24 *)( (f24 *)D - ( M & 1 ) );

  z_24x2  = AE_ZEROP48();
  z_f24x2 = ( z_24x2 );

  for ( m=0; m<3; m++ )
  {
    // Q(23+8) <- Q23 + 8
    AE_S32F24_L_IP( z_f24x2, D_s, +4 );
  }

  //----------------------------------------------------------------------------
  // Compute (N&~7) convolution results.

  X = (const ae_f24x2 *)x_buf;
  R = (      ae_f24x2 *)r;

  R_va = AE_ZALIGN64();

  _M = M + ((-M)&3);

  //
  // Process vector x data in 8-entries blocks.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 8 left-most unprocessed x[] entries, the first is x[8*n-(M-1)].
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( x0, X, +8 );
    AE_L32X2F24_IP( x1, X, +8 );
    AE_L32X2F24_IP( x2, X, +8 );
    AE_L32X2F24_IP( x3, X, +8 );

    // Use the shuttle pointer when computing the convolution. Preserve the X
    // pointer for the next iteration.
    S = X;

    //
    // Inner loop prologue: process first 4 y[] entries for 8 accumulators.
    //

    // Load 4 more x[] entries. Altogether we have 12 x[] entries stored in
    // 6 AE registers.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( x4, S, +8 );
    AE_L32X2F24_IP( x5, S, +8 );

    Y = (const ae_f24x2 *)y_buf;

    // Load y[M-1]..y[M-1-3].
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( y1, Y, +8 );
    AE_L32X2F24_IP( y0, Y, +8 );

    // Q16.47 <- [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
    AE_MULFD24X2_FIR_H( q0, q1, x0, x1, y1 );
    AE_MULFD24X2_FIR_H( q2, q3, x1, x2, y1 );
    AE_MULFD24X2_FIR_H( q4, q5, x2, x3, y1 );
    AE_MULFD24X2_FIR_H( q6, q7, x3, x4, y1 );

    // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
    AE_MULAFD24X2_FIR_H( q0, q1, x1, x2, y0 );
    AE_MULAFD24X2_FIR_H( q2, q3, x2, x3, y0 );
    AE_MULAFD24X2_FIR_H( q4, q5, x3, x4, y0 );
    AE_MULAFD24X2_FIR_H( q6, q7, x4, x5, y0 );

    // 4 x[] entries are done, shift them out of the registers.
    x0 = x2; x1 = x3; 
    x2 = x4; x3 = x5;

    //
    // Inner loop kernel: process 4 y[] entries for 8 accumulators. 12-entries 
    // register delay line is updated similarly to the loop prologue.
    //

    for ( m=0; m<(_M>>2)-1; m++ )
    {
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( x4, S, +8 );
      AE_L32X2F24_IP( x5, S, +8 );

      // Load y[(M-1)-4*(m+1)]..y[(M-1)-4*(m+1)-3].
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( y1, Y, +8 );
      AE_L32X2F24_IP( y0, Y, +8 );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, x0, x1, y1 );
      AE_MULAFD24X2_FIR_H( q2, q3, x1, x2, y1 );
      AE_MULAFD24X2_FIR_H( q4, q5, x2, x3, y1 );
      AE_MULAFD24X2_FIR_H( q6, q7, x3, x4, y1 );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, x1, x2, y0 );
      AE_MULAFD24X2_FIR_H( q2, q3, x2, x3, y0 );
      AE_MULAFD24X2_FIR_H( q4, q5, x3, x4, y0 );
      AE_MULAFD24X2_FIR_H( q6, q7, x4, x5, y0 );

      // 4 x[] entries are done, shift them out of the registers.
      x0 = x2; x1 = x3;
      x2 = x4; x3 = x5;
    }

    //
    // Convert and save 8 convolution results.
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
  // (N&~7) convolution results are done by now. If N - (N&~7) >= 4, then
  // compute a block of 4 results.

  if ( N & 4 )
  {
    // Load 4 left-most unprocessed x[] entries.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( x0, X, +8 );
    AE_L32X2F24_IP( x1, X, +8 );

    // Use the shuttle pointer when computing the convolution. Preserve the X
    // pointer for the last block of (N&3) results.
    S = X;

    //
    // Inner loop prologue: process first 4 y[] entries for 4 accumulators.
    //

    // Load 4 more x[] entries. Altogether we have 8 x[] entries residing in
    // 4 AE registers.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( x2, S, +8 );
    AE_L32X2F24_IP( x3, S, +8 );

    Y = (const ae_f24x2 *)y_buf;

    // Load y[M-1]..y[M-1-3].
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( y1, Y, +8 );
    AE_L32X2F24_IP( y0, Y, +8 );

    // Q16.47 <- [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
    AE_MULFD24X2_FIR_H( q0, q1, x0, x1, y1 );
    AE_MULFD24X2_FIR_H( q2, q3, x1, x2, y1 );

    // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
    AE_MULAFD24X2_FIR_H( q0, q1, x1, x2, y0 );
    AE_MULAFD24X2_FIR_H( q2, q3, x2, x3, y0 );

    // 4 x[] entries are done, shift them out of the registers.
    x0 = x2; x1 = x3;

    //
    // Inner loop kernel: process 4 y[] entries for 4 accumulators. 8-entries 
    // register delay line is updated similarly to the loop prologue.
    //

    for ( m=0; m<(_M>>2)-1; m++ )
    {
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( x2, S, +8 );
      AE_L32X2F24_IP( x3, S, +8 );

      // Load y[(M-1)-4*(m+1)]..y[(M-1)-4*(m+1)-3].
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( y1, Y, +8 );
      AE_L32X2F24_IP( y0, Y, +8 );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, x0, x1, y1 );
      AE_MULAFD24X2_FIR_H( q2, q3, x1, x2, y1 );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, x1, x2, y0 );
      AE_MULAFD24X2_FIR_H( q2, q3, x2, x3, y0 );

      // 4 x[] entries are done, shift them out of the registers.
      x0 = x2; x1 = x3;
    }

    //
    // Convert and save 4 convolution results.
    //

    // Q23 <- Q16.47 - 24 w/ rounding and saturation.
    r0 = AE_ROUND24X2F48SASYM( q0, q1 );
    r1 = AE_ROUND24X2F48SASYM( q2, q3 );

    // Q(23+8) <- Q23 + 8
    AE_SA32X2F24_IP( r0, R_va, R );
    AE_SA32X2F24_IP( r1, R_va, R );
  }

  AE_SA64POS_FP( R_va, R );

  //----------------------------------------------------------------------------
  // (N&~3) convolution results are done by now. If N is not a multiple of 4,
  // compute the last (N&3) results.

  if ( N & 3 )
  {
    // Load 4 left-most unprocessed x[] entries.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( x0, X, +8 );
    AE_L32X2F24_IP( x1, X, +8 );

    //
    // Inner loop prologue: process first 4 y[] entries for 4 accumulators.
    //

    // Load 4 more x[] entries. Altogether we have 8 x[] entries residing in
    // 4 AE registers.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( x2, X, +8 );
    AE_L32X2F24_IP( x3, X, +8 );

    Y = (const ae_f24x2 *)y_buf;

    // Load y[M-1]..y[M-1-3].
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( y1, Y, +8 );
    AE_L32X2F24_IP( y0, Y, +8 );

    // Q16.47 <- [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
    AE_MULFD24X2_FIR_H( q0, q1, x0, x1, y1 );
    AE_MULFD24X2_FIR_H( q2, q3, x1, x2, y1 );

    // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
    AE_MULAFD24X2_FIR_H( q0, q1, x1, x2, y0 );
    AE_MULAFD24X2_FIR_H( q2, q3, x2, x3, y0 );

    // 4 x[] entries are done, shift them out of the registers.
    x0 = x2; x1 = x3;

    //
    // Inner loop kernel: process 4 y[] entries for 4 accumulators. 8-entries 
    // register delay line is updated similarly to the loop prologue.
    //

    for ( m=0; m<(_M>>2)-1; m++ )
    {
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( x2, X, +8 );
      AE_L32X2F24_IP( x3, X, +8 );

      // Load y[(M-1)-4*(m+1)]..y[(M-1)-4*(m+1)-3].
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( y1, Y, +8 );
      AE_L32X2F24_IP( y0, Y, +8 );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, x0, x1, y1 );
      AE_MULAFD24X2_FIR_H( q2, q3, x1, x2, y1 );

      // Q16.47 <- Q16.47 + [ Q23*Q23 + 1 ] + [ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, x1, x2, y0 );
      AE_MULAFD24X2_FIR_H( q2, q3, x2, x3, y0 );

      // 4 x[] entries are done, shift them out of the registers.
      x0 = x2; x1 = x3;
    }

    //
    // Convert and save (N&3) convolution results.
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

} /* fir_convola24x24() */
