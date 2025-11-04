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
    Real block FIR filter, 32x16-bit, unaligned data and arbitrary M/N allowed
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
#define MAGIC     0x463c5eca

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
      ( (size_t)(size) + (align) - 1 )

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
      (void*)( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

#define sz_i16   sizeof(int16_t)
#define sz_i32   sizeof(int32_t)

/* Filter instance structure. */
typedef struct tag_bkfira32x16_t
{
  uint32_t        magic;     // Instance pointer validation number
  int             M;         // Number of filter coefficients
  const int16_t * coef;      // M filter coefficients, aligned
  int32_t *       delayLine; // Delay line for samples, aligned
  int             delayLen;  // Delay line length, in samples
  int32_t *       delayPos;  // Delay line slot to be filled next

} bkfira32x16_t, *bkfira32x16_ptr_t;

/* Calculate the memory block size for an FIR filter with given attributes. */
size_t bkfira32x16_alloc( int M )
{
  ASSERT( M > 0 );

  return ( ALIGNED_SIZE( sizeof( bkfira32x16_t ), 4 )
           + // Delay line
           ALIGNED_SIZE( ( M + (-M&3) + 8 )*sz_i32, 8 )
           + // Filter coefficients
           ALIGNED_SIZE( ( M + (-M&3) + 4 )*sz_i16, 8 ) );

} /* bkfira32x16_alloc() */

/* Initialize the filter structure. The delay line is zeroed. */
bkfira32x16_handle_t bkfira32x16_init( void *              objmem, 
                                       int                 M, 
                                 const int16_t * restrict  h )
{
  bkfira32x16_ptr_t bkfir;
  void *            ptr;
  int32_t *         delLine;
  int16_t *         coef;
  int16_t *         p_coef;

  int m, _M;

  ASSERT( objmem && M > 0 && h );

  //
  // Partition the memory block
  //

  // Complement the filter length to the next multiple of 4.
  _M = M + (-M&3);

  ptr     = objmem;
  bkfir   = (bkfira32x16_ptr_t)ALIGNED_ADDR( ptr, 4 );
  ptr     = bkfir + 1;
  delLine = (int32_t *)ALIGNED_ADDR( ptr, 8 );
  ptr     = delLine + _M + 8;
  coef    = (int16_t *)ALIGNED_ADDR( ptr, 8 );
  ptr     = coef + _M + 4;

  ASSERT( (int8_t*)ptr - (int8_t*)objmem <= (int)bkfira32x16_alloc( M ) );

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
  bkfir->delayLen  = _M+8;
  bkfir->delayPos  = delLine;

  return (bkfir);

} /* bkfira32x16_init() */

/* Put a chunk of input signal into the delay line and compute the filter
 * response. */
void bkfira32x16_process( bkfira32x16_handle_t _bkfir, 
                          int32_t * restrict   y,
                    const int32_t * restrict   x, int N )
{
  bkfira32x16_ptr_t bkfir = (bkfira32x16_ptr_t)_bkfir;

  const ae_f16x4   *          C;
        ae_int32x2 * restrict D_wr;
        ae_int32   * restrict D_s_wr;
  const ae_int32x2 *          D_rd;
  const ae_int32   *          D_s_rd;
  const ae_int32x2 *          X;
  const ae_int32   *          X_s;
        ae_f32     * restrict Y;

  ae_valign D_va, X_va;

  ae_int32x2 t0, t1, t2, t3, t4, t5;
  ae_f32x2   d0, d1, d2, d3, d4, d5;
  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_f16x4   c;

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

  D_wr = (      ae_int32x2*)bkfir->delayPos;
  X    = (const ae_int32x2*)x;
  Y    = (      ae_f32    *)y;

  WUR_AE_CBEGIN0( (uintptr_t)( bkfir->delayLine                   ) );
  WUR_AE_CEND0  ( (uintptr_t)( bkfir->delayLine + bkfir->delayLen ) );

  X_va = AE_LA64_PP( X );

  //
  // Break the input signal into 8-samples blocks. For each block, store 8
  // samples to the delay line and compute the filter response.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 8 input samples.
    // Q31
    AE_LA32X2_IP( t0, X_va, X );
    AE_LA32X2_IP( t1, X_va, X );
    AE_LA32X2_IP( t2, X_va, X );
    AE_LA32X2_IP( t3, X_va, X );

    D_va = AE_ZALIGN64();

    // Store 8 samples to the delay line buffer with circular address update.
    // Q31
    AE_SA32X2_IC( t0, D_va, D_wr );
    AE_SA32X2_IC( t1, D_va, D_wr );
    AE_SA32X2_IC( t2, D_va, D_wr );
    AE_SA32X2_IC( t3, D_va, D_wr );

    AE_SA64POS_FP( D_va, D_wr );

    // Circular buffer pointer looks at the oldest sample: M+8 samples back from
    // the newest one.
    D_rd = D_wr;

    AE_LA32X2POS_PC( D_va, D_rd );

    // Load 8 oldest samples from the delay line.
    // Q31
    AE_LA32X2_IC( t0, D_va, D_rd );
    AE_LA32X2_IC( t1, D_va, D_rd );
    AE_LA32X2_IC( t2, D_va, D_rd );
    AE_LA32X2_IC( t3, D_va, D_rd );

    d0 = ( t0 );
    d1 = ( t1 );
    d2 = ( t2 );
    d3 = ( t3 );

    //
    // Inner loop prologue: process the first 4 taps for each of 8 accumulators.
    //

    // Load next 4 samples. Altogether we have 12 samples residing in 6 AE
    // registers.
    // Q31
    AE_LA32X2_IC( t4, D_va, D_rd );
    AE_LA32X2_IC( t5, D_va, D_rd );

    d4 = ( t4 );
    d5 = ( t5 );

    // Reset the coefficients pointer. Now it looks at the tap corresponding
    // to the oldest sample in the delay line.
    C = (const ae_f16x4*)bkfir->coef;

    // Load 4 tap coefficients.
    // Q15
    ae_f16x4_loadip( c, C, +8 );

    // 2xQ16.47 <- 2x[ Q31*Q15 + 1 ] + 2x[ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
    AE_MULFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
    AE_MULFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
    AE_MULFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

    // 2xQ16.47 <- 2xQ16.47 + 2x[ Q31*Q15 + 1 ] + 2x[ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
    AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
    AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
    AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

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
      // Q31
      AE_LA32X2_IC( t4, D_va, D_rd );
      AE_LA32X2_IC( t5, D_va, D_rd );

      d4 = ( t4 );
      d5 = ( t5 );

      // Load the next 4 tap coefficients.
      ae_f16x4_loadip( c, C, +8 );

      // 2xQ16.47 <- 2xQ16.47 + 2x[ Q31*Q15 + 1 ] + 2x[ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
      AE_MULAFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

      // 2xQ16.47 <- 2xQ16.47 + 2x[ Q31*Q15 + 1 ] + 2x[ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
      AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

      // 4 taps are done. Move 4 input samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    // Convert and save 8 outputs.
    // 2xQ31 <- 2xQ16.47 - 16 w/ rounding and saturation.
    AE_S32RA64S_IP( q0, Y, +4 );
    AE_S32RA64S_IP( q1, Y, +4 );
    AE_S32RA64S_IP( q2, Y, +4 );
    AE_S32RA64S_IP( q3, Y, +4 );
    AE_S32RA64S_IP( q4, Y, +4 );
    AE_S32RA64S_IP( q5, Y, +4 );
    AE_S32RA64S_IP( q6, Y, +4 );
    AE_S32RA64S_IP( q7, Y, +4 );
  }

  //
  // Process the last N&7 samples. 
  //

  X_s    = (const ae_int32*)X;
  D_s_wr = (      ae_int32*)D_wr;

  if ( N&7 )
  {
    // Insert 0..7 input samples into the delay line one-by-one.
    for ( n=0; n<(N&7); n++ )
    {
      AE_L32_IP ( t0, X_s   , +4 );
      AE_S32_L_XC( t0, D_s_wr, +4 );
    }

    D_s_rd = D_s_wr;

    // Perform dummy reads to skip 8-(N&7) oldest samples.
    for ( n=0; n<(-N&7); n++ )
    {
      AE_L32_XC( t0, D_s_rd, +4 );
    }

    D_rd = (const ae_int32x2*)D_s_rd;

    AE_LA32X2POS_PC( D_va, D_rd );

    // Load 8 oldest samples from the delay line.
    // Q31
    AE_LA32X2_IC( t0, D_va, D_rd );
    AE_LA32X2_IC( t1, D_va, D_rd );
    AE_LA32X2_IC( t2, D_va, D_rd );
    AE_LA32X2_IC( t3, D_va, D_rd );

    d0 = ( t0 );
    d1 = ( t1 );
    d2 = ( t2 );
    d3 = ( t3 );

    //
    // Inner loop prologue: process the first 4 taps for each of 8 accumulators.
    //

    // Load next 4 samples. Altogether we have 12 samples residing in 6 AE
    // registers.
    // Q31
    AE_LA32X2_IC( t4, D_va, D_rd );
    AE_LA32X2_IC( t5, D_va, D_rd );

    d4 = ( t4 );
    d5 = ( t5 );

    // Reset the coefficients pointer. Now it looks at the tap corresponding
    // to the oldest sample in the delay line.
    C = (const ae_f16x4*)bkfir->coef;

    // Load 4 tap coefficients.
    // Q15
    ae_f16x4_loadip( c, C, +8 );

    // 2xQ16.47 <- 2x[ Q31*Q15 + 1 ] + 2x[ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
    AE_MULFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
    AE_MULFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
    AE_MULFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

    // 2xQ16.47 <- 2xQ16.47 + 2x[ Q31*Q15 + 1 ] + 2x[ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
    AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
    AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
    AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

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
      // Q31
      AE_LA32X2_IC( t4, D_va, D_rd );
      AE_LA32X2_IC( t5, D_va, D_rd );

      d4 = ( t4 );
      d5 = ( t5 );

      // Load the next 4 tap coefficients.
      ae_f16x4_loadip( c, C, +8 );

      // 2xQ16.47 <- 2xQ16.47 + 2x[ Q31*Q15 + 1 ] + 2x[ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
      AE_MULAFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

      // 2xQ16.47 <- 2xQ16.47 + 2x[ Q31*Q15 + 1 ] + 2x[ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
      AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

      // 4 taps are done. Move 4 input samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    Y = (ae_f32*)( y + N-1 );

    // Convert and save 1..7 outputs.
    // 2xQ31 <- 2xQ16.47 - 16 w/ rounding and saturation.
    switch ( N & 7 )
    {
    case 7: AE_S32RA64S_IP( q6, Y, -4 );
    case 6: AE_S32RA64S_IP( q5, Y, -4 );
    case 5: AE_S32RA64S_IP( q4, Y, -4 );
    case 4: AE_S32RA64S_IP( q3, Y, -4 );
    case 3: AE_S32RA64S_IP( q2, Y, -4 );
    case 2: AE_S32RA64S_IP( q1, Y, -4 );
    }

    AE_S32RA64S_IP( q0, Y, -4 );
  }

  bkfir->delayPos = (int32_t*)D_s_wr;

} /* bkfira32x16_process() */
