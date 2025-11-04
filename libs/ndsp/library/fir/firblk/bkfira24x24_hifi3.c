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
    Real block FIR filter, 24x24-bit, unaligned data and arbitrary M/N allowed
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Real FIR filter.
  Computes a real FIR filter (direct-form) using IR stored in vector h. The 
  real data input is stored in vector x. The filter output result is stored 
  in vector y. The filter calculates N output samples using M coefficients 
  and requires last M-1 samples in the delay line.
  These functions implement FIR filter with no limitation on size of data
  block, alignment and length of impulse response at the cost of increased
  processing complexity.
  NOTE: 
  User application is not responsible for management of delay lines.

  Precision: 
  16x16    16-bit data, 16-bit coefficients, 16-bit outputs
  24x24    24-bit data, 24-bit coefficients, 24-bit outputs
  32x16    32-bit data, 16-bit coefficients, 32-bit outputs
  32x32    32-bit data, 32-bit coefficients, 32-bit outputs
  f        floating point
  Input:
  x[N]     input samples, Q15, Q31, floating point
  h[M]     filter coefficients in normal order, Q15, Q31, floating point
  N        length of sample block
  M        length of filter
  Output:
  y[N]     input samples, Q15, Q31, floating point 

  Restrictions:
  x,y      should not be overlapping
-------------------------------------------------------------------------*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

/* Instance pointer validation number. */
#define MAGIC     0xe96e9f66

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
      ( (size_t)(size) + (align) - 1 )

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
      (void*)( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

#define sz_f24   sizeof(f24)

#define MIN(a,b)    ( (a) < (b) ? (a) : (b) )

/* Filter instance structure. */
typedef struct tag_bkfira24x24_t
{
  uint32_t    magic;     // Instance pointer validation number
  int         M;         // Number of filter coefficients
  const f24 * coef;      // M filter coefficients, reverted and aligned
  f24 *       delayLine; // Delay line for samples
  int         delayLen;  // Delay line length, in samples
  f24 *       delayPos;  // Delay line slot to be filled next

} bkfira24x24_t, *bkfira24x24_ptr_t;

/* Calculate the memory block size for an FIR filter with given attributes. */
size_t bkfira24x24_alloc( int M )
{
  ASSERT( M > 0 );

  return ( ALIGNED_SIZE( sizeof( bkfira24x24_t ), 4 )
           + // Delay line
           ALIGNED_SIZE( ( M + (-M&3) + 8 )*sz_f24, 8 )
           + // Filter coefficients
           ALIGNED_SIZE( ( M + (-M&3) + 4 )*sz_f24, 8 ) );

} /* bkfira24x24_alloc() */

/* Initialize the filter structure. The delay line is zeroed. */
bkfira24x24_handle_t bkfira24x24_init( void *         objmem, 
                                       int            M, 
                                 const f24 * restrict h )
{
  bkfira24x24_ptr_t bkfir;
  void *            ptr;
  f24 *             coef;
  f24 *             p_coef;
  f24 *             delLine;

  int m, _M;

  ASSERT( objmem && (M>0) && h );

  //
  // Partition the memory block
  //

  // Complement the filter length to the next multiple of 4.
  _M = M + (-M&3);

  ptr     = objmem;
  bkfir   = (bkfira24x24_ptr_t)ALIGNED_ADDR( ptr, 4 );
  ptr     = bkfir + 1;
  delLine = (f24 *)ALIGNED_ADDR( ptr, 8 );
  ptr     = delLine + _M + 8;
  coef    = (f24 *)ALIGNED_ADDR( ptr, 8 );
  ptr     = coef + _M + 4;

  NASSERT( (int8_t*)ptr - (int8_t*)objmem <= (int)bkfira24x24_alloc( M ) );

  //
  // Convert and copy the filter coefficients. The conversion may be described
  // in three steps:
  //  1. Filter length is extended to the next multiple of 4 by appending a
  //     number of zero coefficients after the last tap, which matches the
  //     oldest sample. This is done to ensure that the delay line buffer size
  //     is always a multiple of 8 bytes. New filter length is _M.
  //  2. Pad the impulse response with 4 zeros. 3 zeros go before the first tap,
  //     and 1 zero follows the last tap. The intention here is to eliminate a 
  //     one-sample delay that would otherwise be introduced into the filter
  //     response, because we use _M+8 delay elements instead of _M+8-1 (8 is
  //     the basic samples block length).
  //  3. Finally, the coefficients order is reverted, and they are copied into
  //     the internal storage.
  //

  p_coef = coef;

  for ( m=0; m<(-M&3)+1; m++ )
  {
    *p_coef++ = 0;
  }

  for ( m=0; m<M; m++ )
  {
    *p_coef++ = h[M-1-m];
  }

  for ( m=0; m<3; m++ )
  {
    *p_coef++ = 0;
  }

  //
  // Zero the delay line.
  //

  for ( m=0; m<_M+8; m++ )
  {
    delLine[m] = 0;
  }

  //
  // Initialize the filter instance.
  //

  bkfir->magic     = MAGIC;
  bkfir->M         = _M;
  bkfir->coef      = coef;
  bkfir->delayLine = delLine;
  bkfir->delayLen  = _M + 8;
  bkfir->delayPos  = delLine;

  return (bkfir);

} /* bkfira24x24_init() */

/* Put a chunk of input signal into the delay line and compute the filter
 * response. */
void bkfira24x24_process( bkfira24x24_handle_t _bkfir, 
                          f24 * restrict       y,
                    const f24 * restrict       x, int N )
{
  bkfira24x24_ptr_t bkfir = (bkfira24x24_ptr_t)_bkfir;

  const ae_f24x2 *          C;
  const ae_f24x2 *          D_rd;
  const ae_f24   *          D_s_rd;
        ae_f24x2 * restrict D_wr;
        ae_f24   * restrict D_s_wr;
  const ae_f24x2 *          X;
  const ae_f24   *          X_s;
        ae_f24x2 * restrict Y;
        ae_f24   * restrict Y_s;

  ae_valign D_va, X_va, Y_va;

  ae_f64     q0, q1, q2, q3;
  ae_f64     q4, q5, q6, q7;
  ae_f24x2   d0, d1, d2, d3;
  ae_f24x2   d4, d5, d6;
  ae_f24x2   c0, c1;

  int M;
  int m, n;

  ASSERT( bkfir && bkfir->magic == MAGIC && y && x );

  M = bkfir->M;

  NASSERT_ALIGN(( bkfir->delayLine                   ),8);
  NASSERT_ALIGN(( bkfir->delayLine + bkfir->delayLen ),8);
  NASSERT_ALIGN(( bkfir->coef                        ),8);
  if(N<=0) return;

  //
  // Setup pointers and circular delay line buffer.
  //

  D_wr = (      ae_f24x2*)bkfir->delayPos;
  X    = (const ae_f24x2*)x;
  Y    = (      ae_f24x2*)y;

  WUR_AE_CBEGIN0( (uintptr_t)( bkfir->delayLine                   ) );
  WUR_AE_CEND0  ( (uintptr_t)( bkfir->delayLine + bkfir->delayLen ) );

  X_va = AE_LA64_PP( X );
  Y_va = AE_ZALIGN64();

  //
  // Break the input signal into 8-samples blocks. For each block, store 8
  // samples to the delay line and compute the filter response.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 8 input samples.
    // Q23 <- Q(23+8) - 8
    AE_LA32X2F24_IP( d0, X_va, X );
    AE_LA32X2F24_IP( d1, X_va, X );
    AE_LA32X2F24_IP( d2, X_va, X );
    AE_LA32X2F24_IP( d3, X_va, X );

    D_va = AE_ZALIGN64();

    // Store 8 samples to the delay line buffer with circular address update.
    // Q(23+8) <- Q23 + 8
    AE_SA32X2F24_IC( d0, D_va, D_wr );
    AE_SA32X2F24_IC( d1, D_va, D_wr );
    AE_SA32X2F24_IC( d2, D_va, D_wr );
    AE_SA32X2F24_IC( d3, D_va, D_wr );

    AE_SA64POS_FP( D_va, D_wr );

    // Circular buffer pointer looks at the oldest sample: M+8 samples back from
    // the newest one.
    D_rd = D_wr;

    AE_LA32X2F24POS_PC( D_va, D_rd );

    // Load 8 oldest samples from the delay line.
    // Q23 <- Q(23+8) - 8
    AE_LA32X2F24_IC( d0, D_va, D_rd );
    AE_LA32X2F24_IC( d1, D_va, D_rd );
    AE_LA32X2F24_IC( d2, D_va, D_rd );
    AE_LA32X2F24_IC( d3, D_va, D_rd );

    //
    // Inner loop prologue: process the first 4 taps for each of 8 accumulators.
    //

    // Load next 4 samples. Altogether we have 12 samples residing in 6 AE
    // registers.
    // Q23 <- Q(23+8) - 8
    AE_LA32X2F24_IC( d4, D_va, D_rd );
    AE_LA32X2F24_IC( d5, D_va, D_rd );

    // Reset the coefficients pointer. Now it looks at the tap corresponding
    // to the oldest sample in the delay line.
    C = (const ae_f24x2*)bkfir->coef;

    // Load 4 tap coefficients.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( c0, C, +8 );
    AE_L32X2F24_IP( c1, C, +8 );

    // 2xQ16.47 <- 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
    AE_MULFD24X2_FIR_H( q0, q1, d0, d1, c0 );
    AE_MULFD24X2_FIR_H( q2, q3, d1, d2, c0 );
    AE_MULFD24X2_FIR_H( q4, q5, d2, d3, c0 );
    AE_MULFD24X2_FIR_H( q6, q7, d3, d4, c0 );

    // 2xQ16.47 <- 2xQ16.47 + 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
    AE_MULAFD24X2_FIR_H( q0, q1, d1, d2, c1 );
    AE_MULAFD24X2_FIR_H( q2, q3, d2, d3, c1 );
    AE_MULAFD24X2_FIR_H( q4, q5, d3, d4, c1 );
    AE_MULAFD24X2_FIR_H( q6, q7, d4, d5, c1 );

    // First 4 taps are done. Move 4 input samples out of the registers file.
    d0 = d2; d1 = d3; d2 = d4; d3 = d5;

    //
    // Inner loop: process 4 taps for 8 accumulators on each trip. Totally we 
    // perform M+4 MACs for each accumulator, 4 of which fall on zero taps
    // inserted into the impulse response during initialization.
    //
    NASSERT(M>=4);
    __Pragma("loop_count min=1")
    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples from the delay line. Altogether we have 12 samples
      // residing in 6 AE registers.
      // Q23 <- Q(23+8) - 8
      AE_LA32X2F24_IC( d4, D_va, D_rd );
      AE_LA32X2F24_IC( d5, D_va, D_rd );

      // Load the next 4 tap coefficients.
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( c0, C, +8 );
      AE_L32X2F24_IP( c1, C, +8 );

      // 2xQ16.47 <- 2xQ16.47 + 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, d0, d1, c0 );
      AE_MULAFD24X2_FIR_H( q2, q3, d1, d2, c0 );
      AE_MULAFD24X2_FIR_H( q4, q5, d2, d3, c0 );
      AE_MULAFD24X2_FIR_H( q6, q7, d3, d4, c0 );

      // 2xQ16.47 <- 2xQ16.47 + 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, d1, d2, c1 );
      AE_MULAFD24X2_FIR_H( q2, q3, d2, d3, c1 );
      AE_MULAFD24X2_FIR_H( q4, q5, d3, d4, c1 );
      AE_MULAFD24X2_FIR_H( q6, q7, d4, d5, c1 );

      // 4 taps are done. Move 4 input samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    // 2xQ23 <- 2xQ16.47 - 24 w/ rounding and saturation.
    d0 = AE_ROUND24X2F48SASYM( q0, q1 );
    d1 = AE_ROUND24X2F48SASYM( q2, q3 );
    d2 = AE_ROUND24X2F48SASYM( q4, q5 );
    d3 = AE_ROUND24X2F48SASYM( q6, q7 );

    // Store 8 filter outputs.
    // 2xQ(23+8) <- 2xQ23 + 8
    AE_SA32X2F24_IP( d0, Y_va, Y );
    AE_SA32X2F24_IP( d1, Y_va, Y );
    AE_SA32X2F24_IP( d2, Y_va, Y );
    AE_SA32X2F24_IP( d3, Y_va, Y );
  }

  AE_SA64POS_FP( Y_va, Y );

  //
  // Process the last N&7 samples. 
  //

  X_s    = (const ae_f24*)X;
  D_s_wr = (      ae_f24*)D_wr;

  if ( N&7 )
  {
    // Insert 1..7 input samples into the delay line one-by-one.
    for ( n=0; n<(N&7); n++ )
    {
      AE_L32F24_IP ( d0, X_s   , +4 );
      AE_S32F24_L_XC( d0, D_s_wr, +4 );
    }

    D_s_rd = D_s_wr;

    // Perform dummy reads to skip 8-(N&7) oldest samples.
    for ( n=0; n<(-N&7); n++ )
    {
      AE_L32F24_XC( d0, D_s_rd, +4 );
    }

    D_rd = (const ae_f24x2*)D_s_rd;

    AE_LA32X2F24POS_PC( D_va, D_rd );

    // Load 8 oldest samples from the delay line.
    // Q23 <- Q(23+8) - 8
    AE_LA32X2F24_IC( d0, D_va, D_rd );
    AE_LA32X2F24_IC( d1, D_va, D_rd );
    AE_LA32X2F24_IC( d2, D_va, D_rd );
    AE_LA32X2F24_IC( d3, D_va, D_rd );

    //
    // Inner loop prologue: process the first 4 taps for each of 8 accumulators.
    //

    // Load next 4 samples. Altogether we have 12 samples residing in 6 AE
    // registers.
    // Q23 <- Q(23+8) - 8
    AE_LA32X2F24_IC( d4, D_va, D_rd );
    AE_LA32X2F24_IC( d5, D_va, D_rd );

    // Reset the coefficients pointer. Now it looks at the tap corresponding
    // to the oldest sample in the delay line.
    C = (const ae_f24x2*)bkfir->coef;

    // Load 4 tap coefficients.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( c0, C, +8 );
    AE_L32X2F24_IP( c1, C, +8 );

    // 2xQ16.47 <- 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
    AE_MULFD24X2_FIR_H( q0, q1, d0, d1, c0 );
    AE_MULFD24X2_FIR_H( q2, q3, d1, d2, c0 );
    AE_MULFD24X2_FIR_H( q4, q5, d2, d3, c0 );
    AE_MULFD24X2_FIR_H( q6, q7, d3, d4, c0 );

    // 2xQ16.47 <- 2xQ16.47 + 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
    AE_MULAFD24X2_FIR_H( q0, q1, d1, d2, c1 );
    AE_MULAFD24X2_FIR_H( q2, q3, d2, d3, c1 );
    AE_MULAFD24X2_FIR_H( q4, q5, d3, d4, c1 );
    AE_MULAFD24X2_FIR_H( q6, q7, d4, d5, c1 );

    // First 4 taps are done. Move 4 input samples out of the registers file.
    d0 = d2; d1 = d3; d2 = d4; d3 = d5;

    //
    // Inner loop: process 4 taps for 8 accumulators on each trip. Totally we 
    // perform M+4 MACs for each accumulator, 4 of which fall on zero taps
    // inserted into the impulse response during initialization.
    //
    NASSERT(M>=4);
    __Pragma("loop_count min=1")
    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples from the delay line. Altogether we have 12 samples
      // residing in 6 AE registers.
      // Q23 <- Q(23+8) - 8
      AE_LA32X2F24_IC( d4, D_va, D_rd );
      AE_LA32X2F24_IC( d5, D_va, D_rd );

      // Load the next 4 tap coefficients.
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( c0, C, +8 );
      AE_L32X2F24_IP( c1, C, +8 );

      // 2xQ16.47 <- 2xQ16.47 + 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, d0, d1, c0 );
      AE_MULAFD24X2_FIR_H( q2, q3, d1, d2, c0 );
      AE_MULAFD24X2_FIR_H( q4, q5, d2, d3, c0 );
      AE_MULAFD24X2_FIR_H( q6, q7, d3, d4, c0 );

      // 2xQ16.47 <- 2xQ16.47 + 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, d1, d2, c1 );
      AE_MULAFD24X2_FIR_H( q2, q3, d2, d3, c1 );
      AE_MULAFD24X2_FIR_H( q4, q5, d3, d4, c1 );
      AE_MULAFD24X2_FIR_H( q6, q7, d4, d5, c1 );

      // 4 taps are done. Move 4 input samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    // 2xQ23 <- 2xQ16.47 - 24 w/ rounding and saturation.
    d0 = AE_ROUND24X2F48SASYM( q0, q0 );
    d1 = AE_ROUND24X2F48SASYM( q1, q1 );
    d2 = AE_ROUND24X2F48SASYM( q2, q2 );
    d3 = AE_ROUND24X2F48SASYM( q3, q3 );
    d4 = AE_ROUND24X2F48SASYM( q4, q4 );
    d5 = AE_ROUND24X2F48SASYM( q5, q5 );
    d6 = AE_ROUND24X2F48SASYM( q6, q6 );

    Y_s = (ae_f24*)( y + N-1 );

    // Store 1..7 filter outputs.
    // 2xQ(23+8) <- 2xQ23 + 8
    switch ( N & 7 )
    {
    case 7: AE_S32F24_L_IP( d6, Y_s, -4 );
    case 6: AE_S32F24_L_IP( d5, Y_s, -4 );
    case 5: AE_S32F24_L_IP( d4, Y_s, -4 );
    case 4: AE_S32F24_L_IP( d3, Y_s, -4 );
    case 3: AE_S32F24_L_IP( d2, Y_s, -4 );
    case 2: AE_S32F24_L_IP( d1, Y_s, -4 );
    }

    AE_S32F24_L_IP( d0, Y_s, -4 );
  }

  bkfir->delayPos = (f24*)D_s_wr;

} /* bkfira24x24_process() */
