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
    Real block FIR filter, 24x24-bit
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Real FIR filter.
  Computes a real FIR filter (direct-form) using IR stored in vector h. The 
  real data input is stored in vector x. The filter output result is stored 
  in vector y. The filter calculates N output samples using M coefficients 
  and requires last M-1 samples in the delay line.
  NOTE: 
  1. User application is not responsible for management of delay lines
  2. User has an option to set IR externally or copy from original location 
     (i.e. from the slower constant memory). In the first case, user is 
     responsible for right alignment, ordering and zero padding of filter 
     coefficients - usually array is composed from zeroes (left padding), 
     reverted IR and right zero padding.


  Precision: 
  16x16    16-bit data, 16-bit coefficients, 16-bit outputs
  24x24    24-bit data, 24-bit coefficients, 24-bit outputs
  24x24p   use 24-bit data packing for internal delay line buffer
           and internal coefficients storage
  32x16    32-bit data, 16-bit coefficients, 32-bit outputs
  32x32    32-bit data, 32-bit coefficients, 32-bit outputs
  f        floating point

  Input:
  x[N]     input samples, Q31, Q15, floating point
  h[M]     filter coefficients in normal order, Q31, Q15, floating point
  N        length of sample block, should be a multiple of 4
  M        length of filter, should be a multiple of 4
  extIR    if zero, IR is copied from original location, otherwise not
           but user should keep alignment, order of coefficients 
           and zero padding requirements shown below
  Output:
  y[N]     output samples, Q31, Q15, floating point

  Alignment, ordering and zero padding for external IR  (extIR!=0)
  ------------------------+----------+--------------+--------------+----------------
  Function                |Alignment,|Left zero     |   Coefficient| Right zero 
                          | bytes    |padding, bytes|   order      | padding, bytes
  ------------------------+----------+--------------+--------------+----------------
  bkfir16x16_init         |     8    |      2       |  inverted    |  6
  bkfir24x24_init         |     8    |      4       |  inverted    |  12
  bkfir24x24p_init        |     8    |((-M&4)+1)*3  |  inverted    |  7
  bkfir32x16_init (M>32)  |     8    |     10       |  inverted    |  6
  bkfir32x16_init (M<=32) |     8    |      2       |  inverted    |  6
  bkfir32x32_init         |     8    |      4       |  inverted    |  12
  bkfirf_init             |     8    |      0       |  direct      |  0
  ------------------------+----------+--------------+--------------+----------------

  Restrictions:
  x, y     should not be overlapping
  x, h     aligned on a 8-bytes boundary
  N, M     multiples of 4 
-------------------------------------------------------------------------*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

/* Instance pointer validation number. */
#define MAGIC     0x9305f5b1

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
      ( (size_t)(size) + (align) - 1 )

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
      (void*)( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

#define sz_f24    sizeof(f24)

/* Filter instance structure. */
typedef struct tag_bkfir24x24_t
{
  uint32_t    magic;     // Instance pointer validation number
  int         M;         // Number of filter coefficients
  const f24 * coef;      // M+4 filter coefficients, reverted and aligned
  f24 *       delayLine; // Delay line for samples
  int         delayLen;  // Delay line length, in samples
  f24 *       delayPos;  // Delay line slot to be filled next

} bkfir24x24_t, *bkfir24x24_ptr_t;

/* Calculate the memory block size for an FIR filter with given attributes. */
size_t bkfir24x24_alloc( int M, int extIR )
{
  NASSERT( M > 0 );
  NASSERT( M%4==0 );

  return ( ALIGNED_SIZE( sizeof( bkfir24x24_t ), 4 )
           + // Delay line
           ALIGNED_SIZE( ( M + 8 )*sz_f24, 8 )
           + // Filter coefficients
           (extIR?0:ALIGNED_SIZE( ( M + 4 )*sz_f24, 8 ) ));

} // bkfir24x24_alloc()

/* Initialize the filter structure. The delay line is zeroed. */
bkfir24x24_handle_t bkfir24x24_init( void *         objmem, 
                                     int            M,
                                     int extIR,
                               const f24 * restrict h )
{
  bkfir24x24_ptr_t bkfir;
  void *           ptr;
  f24 *            coef;
  f24 *            delLine;

  int m;

  NASSERT( objmem && M > 0 && h );
  NASSERT(  M%4==0 && IS_ALIGN( h ) );

  //
  // Partition the memory block
  //

  ptr     = objmem;
  bkfir   = (bkfir24x24_ptr_t)ALIGNED_ADDR( ptr, 4 );
  ptr     = bkfir + 1;
  delLine = (f24 *)ALIGNED_ADDR( ptr, 8 );
  ptr     = delLine + M + 8;
  if (extIR)
  {
      coef = (f24 *)h;
  }
  else
  {
      coef    = (f24 *)ALIGNED_ADDR( ptr, 8 );
      ptr     = coef + M + 4;
  }

  NASSERT( (int8_t*)ptr - (int8_t*)objmem <= (int)bkfir24x24_alloc( M, extIR ) );

  //
  // Copy the filter coefficients and zero the delay line. Original impulse
  // response is padded with zeros: three zeros go before the first tap
  // (corresponds to the newest sample), one zero follows the last tap,
  // which matches the oldest sample. After that the order of filter
  // coefficients is reverted.
  //
  if (extIR == 0)
  {
      coef[0] = 0;

      for ( m=1; m<M+1; m++ )
      {
        coef[m] = h[M-m];
      }

      for ( ; m<M+4; m++ )
      {
        coef[m] = 0;
      }
  }
  for ( m=0; m<M+8; m++ )
  {
    delLine[m] = 0;
  }

  //
  // Initialize the filter instance.
  //

  bkfir->magic     = MAGIC;
  bkfir->M         = M;
  bkfir->coef      = coef;
  bkfir->delayLine = delLine;
  bkfir->delayLen  = M + 8;
  bkfir->delayPos  = delLine;

  return (bkfir);

} // bkfir24x24_init()

/* Put a chunk of input signal into the delay line and compute the filter
 * response. */
void bkfir24x24_process( bkfir24x24_handle_t _bkfir, 
                         f24 * restrict      y,
                   const f24 * restrict      x, int N )
{
  bkfir24x24_ptr_t bkfir = (bkfir24x24_ptr_t)_bkfir;

  const ae_f24x2 *          C;
  const ae_f24x2 *          D_rd;
        ae_f24x2 * restrict D_wr;
  const ae_f24x2 *          X;
        ae_f24x2 * restrict Y;

  ae_valign Y_va;

  ae_f64   q0, q1, q2, q3;
  ae_f64   q4, q5, q6, q7;
  ae_f24x2 d0, d1, d2, d3;
  ae_f24x2 d4, d5;
  ae_f24x2 c0, c1;

  int M;
  int m, n;

  NASSERT( bkfir && bkfir->magic == MAGIC && y && x );
  NASSERT( IS_ALIGN( x ) );

  M = bkfir->M;
  NASSERT_ALIGN((bkfir->delayLine                  ),8);
  NASSERT_ALIGN((bkfir->delayLine + bkfir->delayLen),8);
  NASSERT_ALIGN((bkfir->delayPos                   ),8);
  NASSERT_ALIGN((bkfir->coef                       ),8);
  NASSERT(M%4==0);
  NASSERT(N%4==0);
  if(N<0) return;
  //
  // Setup pointers and circular delay line buffer.
  //

  D_wr = (      ae_f24x2*)bkfir->delayPos;
  X    = (const ae_f24x2*)x;
  Y    = (      ae_f24x2*)y;

  WUR_AE_CBEGIN0( (uintptr_t)( bkfir->delayLine                   ) );
  WUR_AE_CEND0  ( (uintptr_t)( bkfir->delayLine + bkfir->delayLen ) );

  Y_va = AE_ZALIGN64();

  //
  // Break the input signal into 8-samples blocks. For each block, store 8
  // samples to the delay line and compute the filter response.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 8 input samples.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( d0, X, +8 );
    AE_L32X2F24_IP( d1, X, +8 );
    AE_L32X2F24_IP( d2, X, +8 );
    AE_L32X2F24_IP( d3, X, +8 );

    // Store 8 samples to the delay line buffer with circular address update.
    // Q(23+8) <- Q23 + 8
    AE_S32X2F24_XC( d0, D_wr, +8 );
    AE_S32X2F24_XC( d1, D_wr, +8 );
    AE_S32X2F24_XC( d2, D_wr, +8 );
    AE_S32X2F24_XC( d3, D_wr, +8 );

    // Circular buffer pointer looks at the oldest sample: M+8 samples back from
    // the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples from the delay line.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_XC( d0, D_rd, +8 );
    AE_L32X2F24_XC( d1, D_rd, +8 ); 
    AE_L32X2F24_XC( d2, D_rd, +8 );
    AE_L32X2F24_XC( d3, D_rd, +8 );

    //
    // Inner loop prologue: process the first 4 taps for each of 8 accumulators.
    //

    // Load next 4 samples. Altogether we have 12 samples residing in 6 AE
    // registers.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_XC( d4, D_rd, +8 );
    AE_L32X2F24_XC( d5, D_rd, +8 );

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

    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples from the delay line. Altogether we have 12 samples
      // residing in 6 AE registers.
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_XC( d4, D_rd, +8 );
      AE_L32X2F24_XC( d5, D_rd, +8 );

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

  //
  // If the signal chunk length N is not a multiple of 8, process the last
  // quadruple: store it to the delay line and compute the filter response.
  //

  if ( N & 4 )
  {
    // Load 4 input samples.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_IP( d0, X, +8 );
    AE_L32X2F24_IP( d1, X, +8 );

    // Store 4 samples to the delay line buffer with circular address update.
    // Q(23+8) <- Q23 + 8
    AE_S32X2F24_XC( d0, D_wr, +8 );
    AE_S32X2F24_XC( d1, D_wr, +8 );

    // Circular buffer pointer looks at the oldest sample: M+8 samples back from
    // the newest one. 
    D_rd = D_wr;

    // Perform dummy reads to jump over 4 oldest 32-bit entries with circular
    // address update. Now the pointer is M+4 samples apart from the newest one.
    AE_L32X2F24_XC( d0, D_rd, +8 );
    AE_L32X2F24_XC( d1, D_rd, +8 );

    // Load 4 samples from the delay line.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_XC( d0, D_rd, +8 );
    AE_L32X2F24_XC( d1, D_rd, +8 );

    //
    // Inner loop prologue: process the first 4 taps for each of 4 accumulators.
    //

    // Load next 4 samples. Altogether we have 8 samples residing in 4 AE 
    // registers.
    // Q23 <- Q(23+8) - 8
    AE_L32X2F24_XC( d2, D_rd, +8 );
    AE_L32X2F24_XC( d3, D_rd, +8 );

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

    // 2xQ16.47 <- 2xQ16.47 + 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
    AE_MULAFD24X2_FIR_H( q0, q1, d1, d2, c1 );
    AE_MULAFD24X2_FIR_H( q2, q3, d2, d3, c1 );

    // First 4 taps are done. Move 4 input samples out of the registers file.
    d0 = d2; d1 = d3;

    //
    // Inner loop: process 4 taps for 4 accumulators on each trip. Totally we 
    // perform M+4 MACs for each accumulator, 4 of which fall on zero taps
    // inserted into the impulse response during initialization.
    //

    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples from the delay line. Altogether we have 8 samples 
      // residing in 4 AE registers.
      AE_L32X2F24_XC( d2, D_rd, +8 );
      AE_L32X2F24_XC( d3, D_rd, +8 );

      // Load the next 4 tap coefficients.
      // Q23 <- Q(23+8) - 8
      AE_L32X2F24_IP( c0, C, +8 );
      AE_L32X2F24_IP( c1, C, +8 );

      // 2xQ16.47 <- 2xQ16.47 + 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, d0, d1, c0 );
      AE_MULAFD24X2_FIR_H( q2, q3, d1, d2, c0 );

      // 2xQ16.47 <- 2xQ16.47 + 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
      AE_MULAFD24X2_FIR_H( q0, q1, d1, d2, c1 );
      AE_MULAFD24X2_FIR_H( q2, q3, d2, d3, c1 );

      // 4 taps are done. Move 4 input samples out of the registers file.
      d0 = d2; d1 = d3;
    }

    // 2xQ23 <- 2xQ16.47 - 24 w/ rounding and saturation.
    d0 = AE_ROUND24X2F48SASYM( q0, q1 );
    d1 = AE_ROUND24X2F48SASYM( q2, q3 );

    // Store 4 filter outputs.
    // 2xQ(23+8) <- 2xQ23 + 8
    AE_SA32X2F24_IP( d0, Y_va, Y );
    AE_SA32X2F24_IP( d1, Y_va, Y );
  }

  AE_SA64POS_FP( Y_va, Y );

  bkfir->delayPos = (f24*)D_wr;

} // bkfir24x24_process()
