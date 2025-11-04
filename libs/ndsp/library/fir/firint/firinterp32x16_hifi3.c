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
    Interpolating block real FIR filter, 32x16-bit
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Interpolating Block Real FIR Filter
  Computes a real FIR filter (direct-form) with interpolation using IR stored 
  in vector h. The real data input is stored in vector x. The filter output 
  result is stored in vector y. The filter calculates N*D output samples 
  using M*D coefficients and requires last N+M*D-1 samples on the delay line.
  NOTE:
  user application is not responsible for management of delay lines

  Precision: 
  16x16     16-bit data, 16-bit coefficients, 16-bit outputs
  24x24     24-bit data, 24-bit coefficients, 24-bit outputs
  32x16     32-bit data, 16-bit coefficients, 32-bit outputs
  32x32     32-bit data, 32-bit coefficients, 32-bit outputs
  f         floating point

  Input:
  h[M*D]    filter coefficients; h[0] is to be multiplied with the 
            newest sample, Q15, Q31, floating point
  D         interpolation ratio
  N         length of input sample block
  M         length of subfilter. Total length of filter is M*D
  x[N]      input samples, Q15, Q31, floating point
  Output:
  y[N*D]    output samples, Q15, Q31, floating point

  Restrictions:
  x,h,y     should not overlap
  x,h       aligned on an 8-bytes boundary
  N         multiple of 8
  M         multiple of 4
  D         should be >1

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  D - 2, 3 or 4
-------------------------------------------------------------------------*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

/* Instance pointer validation number. */
#define MAGIC     0x54352160

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
      ( (size_t)(size) + (align) - 1 )

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
      (void*)( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

#define sz_i16    sizeof(int16_t)
#define sz_i32    sizeof(int32_t)

/* Data processing function of a particular interpolating filter. Stores a
 * block of N input samples to the circular delay line buffer and computes
 * N*D samples of interpolating FIR filter's response.
 * Input:
 *   cbegin  - circular delay line buffer start address
 *   cend    - points at the next to the last byte of the circular buffer
 *   cpos    - next position in the buffer to be filled with an input sample
 *   x[N]    - input samples
 *   h[]     - decimating FIR filter coefficients, array layout varies
 * Output:
 *   y[N*D]  - output samples
 *   retval  - updated circular buffer's next position pointer
 * Notes and restrictions:
 *   1. Most of data processing functions feature a single, hard-coded 
 *      interpolation factor, so they expect a determined value for parameter D.
 *   2. All pointers with the exception of y[N] must be aligned on an 8-bytes
 *      boundary.
 *   3. N - must be a multiple of 8.
 *   4. M - must be a multiple of 4. */
typedef int32_t * (proc_fxn_t)( int32_t * restrict y,
                                int32_t *          cbegin,
                                int32_t *          cend,
                                int32_t * restrict cpos,
                          const int32_t * restrict x,
                          const int16_t * restrict h,
                          int D, int N, int M );

static proc_fxn_t firinterpX_proc;
static proc_fxn_t firinterp2_proc;
static proc_fxn_t firinterp3_proc;
static proc_fxn_t firinterp4_proc;

/* Interpolator instance structure. */
typedef struct tag_firinterp32x16_t
{
  uint32_t        magic;     // Instance pointer validation number
  int             D;         // Interpolation factor
  int             M;         // Number of filter coefficients
  proc_fxn_t *    procFxn;   // Filter data processing function
  const int16_t * coef;      // Filter coefficients
  int32_t *       delayLine; // Delay line for complex samples
  int             delayLen;  // Delay line length, in complex samples
  int32_t *       delayPos;  // Delay line slot to be filled next

} firinterp32x16_t, *firinterp32x16_ptr_t;

/* Calculate the memory block size for an interpolator with given 
 * attributes. */
size_t firinterp32x16_alloc( int D, int M )
{
  int delayLen, coefNum;

  ASSERT( D > 1  && M > 0 );

  ASSERT(!( M & 3 ) );

  delayLen = ( M + 8 );
  coefNum  = ( M + 4 )*D;

  return ( ALIGNED_SIZE( sizeof( firinterp32x16_t ), 4 )
           + // Delay line
           ALIGNED_SIZE( delayLen*sz_i32, 8 )
           + // Coefficients
           ALIGNED_SIZE( coefNum*sz_i16, 8 ) );

} // firinterp32x16_alloc()

/* Initialize the interpolator structure. The delay line is zeroed. */
firinterp32x16_handle_t firinterp32x16_init( void * objmem, 
                                             int D, int M, 
                                             const int16_t * restrict h )
{
  firinterp32x16_ptr_t firinterp;
  void *               ptr;
  int32_t *            delayLine;
  int                  delayLen;
  int16_t *            coefBuf;
  int16_t *            coefBank;
  int                  coefNum;
  proc_fxn_t *         procFxn;

  int d, m;

  NASSERT( objmem && D > 1 && M > 0 && h );

  NASSERT(  !( M & 3 ) && IS_ALIGN( h ) );

  //
  // Select the processing function, delay line length and coefficients
  // block layout.
  //

  delayLen = ( M + 8 );
  coefNum  = ( M + 4 )*D;

  procFxn = ( D == 2 ? &firinterp2_proc :
              D == 3 ? &firinterp3_proc :
              D == 4 ? &firinterp4_proc :
                       &firinterpX_proc );

  //
  // Partition the memory block.
  //
  ptr       = objmem;
  firinterp = (firinterp32x16_ptr_t)(ALIGNED_ADDR( ptr, 4 ));
  ptr       = firinterp + 1;
  delayLine = (int32_t*)(ALIGNED_ADDR( ptr, 8 ));
  ptr       = delayLine + delayLen;
  coefBuf   = (int16_t*)(ALIGNED_ADDR( ptr, 8 ));
  ptr       = coefBuf + coefNum;

  NASSERT( (int8_t*)ptr - (int8_t*)objmem <= 
          (int)firinterp32x16_alloc( D, M ) );

  //
  // Break the impulse response into D coefficients banks and copy them in
  // reverted order.
  //

  for ( d=0; d<D; d++ )
  {
    coefBank = coefBuf + d*(M+4);

    // To avoid a 1-sample delay, insert a zero coefficient that will match the
    // oldest sample in the delay line. To keep the filter length a multiple of
    // 4, append 3 zeros after the last coefficient.
    coefBank[0] = coefBank[M+1] = coefBank[M+2] = coefBank[M+3] = 0;

    // Copy bank's coefficients in reverted order.
    for ( m=1; m<=M; m++ )
    {
      coefBank[m] = h[D*(M-m)+d];
    }
  }

  //
  // Zero the delay line.
  //

  for ( m=0; m<delayLen; m++ )
  {
    delayLine[m] = 0;
  }

  //
  // Initialize the interpolator instance.
  //

  firinterp->magic     = MAGIC;
  firinterp->D         = D;
  firinterp->M         = M;
  firinterp->procFxn   = procFxn;
  firinterp->coef      = coefBuf;
  firinterp->delayLine = delayLine;
  firinterp->delayLen  = delayLen;
  firinterp->delayPos  = delayLine;

  return (firinterp);

} // firinterp32x16_init()

/* Put a chunk of input signal into the delay line and compute the filter
 * response. */
void firinterp32x16_process( firinterp32x16_handle_t _firinterp, 
                             int32_t * restrict      y,
                       const int32_t * restrict      x, int N )
{
  firinterp32x16_ptr_t firinterp = (firinterp32x16_ptr_t)_firinterp;

  NASSERT( firinterp && firinterp->magic == MAGIC && y && x );

  NASSERT( IS_ALIGN( x ) );

  //
  // Call filter's data processing function. It will store the block of input 
  // samples to the delay line, and compute the filter response. Returns the
  // updated next position pointer into the delay line buffer.
  //

  NASSERT( firinterp->procFxn );

  if(N>0)
  firinterp->delayPos = (*firinterp->procFxn)( 
                                   y,
                                   firinterp->delayLine,
                                   firinterp->delayLine + firinterp->delayLen,
                                   firinterp->delayPos,
                                   x,
                                   firinterp->coef,
                                   firinterp->D,
                                   N,
                                   firinterp->M );

} // firinterp32x16_process()

/* Data processing function for a factor 2 interpolating FIR filter. */
static int32_t * firinterp2_proc( int32_t * restrict y,
                                  int32_t *          cbegin,
                                  int32_t *          cend,
                                  int32_t * restrict cpos,
                          const int32_t * restrict x,
                          const int16_t * restrict h,
                          int D, int N, int M )
{
        ae_int32x2 * restrict D_wr;
  const ae_int32x2 *          D_rd;
  const ae_int32x2 *          X;
        ae_f32     * restrict Y;
  const ae_f16x4   *          C;

  ae_f32x2   d0, d1, d2, d3, d4, d5;
  ae_int32x2 t0, t1, t2, t3, t4, t5;
  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_int64   u0, u1, u2, u3, u4, u5, u6, u7;
  ae_f16x4   c;

  int n, m;

  NASSERT( y && cbegin && cend && cpos && x && h && (D == 2) && (N > 0) && (M > 0) );

  NASSERT( !( N & 7 ) && !( M & 3 ) );

  NASSERT( IS_ALIGN( cbegin ) &&
          IS_ALIGN( cend   ) &&
          IS_ALIGN( cpos   ) &&
          IS_ALIGN( x      ) &&
          IS_ALIGN( h      ) );

  //
  // Setup pointers and circular delay line buffer.
  //

  D_wr = (      ae_int32x2 *)cpos;
  X    = (const ae_int32x2 *)x;
  Y    = (      ae_f32     *)y;

  WUR_AE_CBEGIN0( (uintptr_t)cbegin );
  WUR_AE_CEND0  ( (uintptr_t)cend   );

  //
  // Break input signal into 8-samples blocks. For each block, store 8 samples
  // to the delay line buffer and compute 16 samples of interpolating FIR filter
  // response.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 8 input samples.
    // Q31
    AE_L32X2_IP( t0, X, +8 );
    AE_L32X2_IP( t1, X, +8 );
    AE_L32X2_IP( t2, X, +8 );
    AE_L32X2_IP( t3, X, +8 );

    // Store input samples into the circular delay line.
    // Q31
    AE_S32X2_XC( t0, D_wr, +8 );
    AE_S32X2_XC( t1, D_wr, +8 );
    AE_S32X2_XC( t2, D_wr, +8 );
    AE_S32X2_XC( t3, D_wr, +8 );

    //--------------------------------------------------------------------------
    // Bank 0
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q31
    AE_L32X2_XC( t0, D_rd, +8 );
    AE_L32X2_XC( t1, D_rd, +8 );
    AE_L32X2_XC( t2, D_rd, +8 );
    AE_L32X2_XC( t3, D_rd, +8 );

    d0 = ( t0 );
    d1 = ( t1 );
    d2 = ( t2 );
    d3 = ( t3 );

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q31
    AE_L32X2_XC( t4, D_rd, +8 );
    AE_L32X2_XC( t5, D_rd, +8 );

    d4 = ( t4 );
    d5 = ( t5 );

    // Reset the coefficients pointer. Now it looks at the oldest sample
    // coefficient, bank #0.
    C = (const ae_f16x4 *)h;

    // Load first 4 tap coefficients of bank #0.
    // Q15
    ae_f16x4_loadip( c, C, +8 );

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
    AE_MULFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
    AE_MULFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
    AE_MULFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
    AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
    AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
    AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3; d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //
    NASSERT(M>=4);
    __Pragma("loop_count min=1")
    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q31
      AE_L32X2_XC( t4, D_rd, +8 );
      AE_L32X2_XC( t5, D_rd, +8 );

      d4 = ( t4 );
      d5 = ( t5 );

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip( c, C, +8 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
      AE_MULAFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
      AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    //
    // Scale outputs by the interpolation factor 2.
    //

    // Q16.47 <- Q16.47 + Q16.47 w/ saturation
    u0 = ( q0 );
    u0 = AE_ADD64S( u0, u0 );
    q0 = ( u0  );

    u1 = ( q1 );
    u1 = AE_ADD64S( u1, u1 );
    q1 = ( u1  );

    u2 = ( q2 );
    u2 = AE_ADD64S( u2, u2 );
    q2 = ( u2  );

    u3 = ( q3 );
    u3 = AE_ADD64S( u3, u3 );
    q3 = ( u3  );

    u4 = ( q4 );
    u4 = AE_ADD64S( u4, u4 );
    q4 = ( u4  );

    u5 = ( q5 );
    u5 = AE_ADD64S( u5, u5 );
    q5 = ( u5  );

    u6 = ( q6 );
    u6 = AE_ADD64S( u6, u6 );
    q6 = ( u6  );

    u7 = ( q7 );
    u7 = AE_ADD64S( u7, u7 );
    q7 = ( u7  );

    //
    // Store 8 outputs to the output array with 2-samples stride.
    //

    // Q31 <- Q16.47 - 16 w/ rounding and saturation
    AE_S32RA64S_IP( q0, Y,  +2*4 ); //  0
    AE_S32RA64S_IP( q1, Y,  +2*4 ); //  2
    AE_S32RA64S_IP( q2, Y,  +2*4 ); //  4
    AE_S32RA64S_IP( q3, Y,  +2*4 ); //  6
    AE_S32RA64S_IP( q4, Y,  +2*4 ); //  8
    AE_S32RA64S_IP( q5, Y,  +2*4 ); // 10
    AE_S32RA64S_IP( q6, Y,  +2*4 ); // 12
    AE_S32RA64S_XP( q7, Y, -13*4 ); // 14

    //--------------------------------------------------------------------------
    // Bank 1
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q31
    AE_L32X2_XC( t0, D_rd, +8 );
    AE_L32X2_XC( t1, D_rd, +8 );
    AE_L32X2_XC( t2, D_rd, +8 );
    AE_L32X2_XC( t3, D_rd, +8 );

    d0 = ( t0 );
    d1 = ( t1 );
    d2 = ( t2 );
    d3 = ( t3 );

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q31
    AE_L32X2_XC( t4, D_rd, +8 );
    AE_L32X2_XC( t5, D_rd, +8 );

    d4 = ( t4 );
    d5 = ( t5 );

    // Load first 4 tap coefficients of bank #1.
    // Q15
    ae_f16x4_loadip( c, C, +8 );

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
    AE_MULFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
    AE_MULFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
    AE_MULFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
    AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
    AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
    AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3; d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //
    NASSERT(M>=4);
    __Pragma("loop_count min=1")
    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q31
      AE_L32X2_XC( t4, D_rd, +8 );
      AE_L32X2_XC( t5, D_rd, +8 );

      d4 = ( t4 );
      d5 = ( t5 );

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip( c, C, +8 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
      AE_MULAFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
      AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    //
    // Scale outputs by the interpolation factor 2.
    //

    // Q16.47 <- Q16.47 + Q16.47 w/ saturation
    u0 = ( q0 );
    u0 = AE_ADD64S( u0, u0 );
    q0 = ( u0  );

    u1 = ( q1 );
    u1 = AE_ADD64S( u1, u1 );
    q1 = ( u1  );

    u2 = ( q2 );
    u2 = AE_ADD64S( u2, u2 );
    q2 = ( u2  );

    u3 = ( q3 );
    u3 = AE_ADD64S( u3, u3 );
    q3 = ( u3  );

    u4 = ( q4 );
    u4 = AE_ADD64S( u4, u4 );
    q4 = ( u4  );

    u5 = ( q5 );
    u5 = AE_ADD64S( u5, u5 );
    q5 = ( u5  );

    u6 = ( q6 );
    u6 = AE_ADD64S( u6, u6 );
    q6 = ( u6  );

    u7 = ( q7 );
    u7 = AE_ADD64S( u7, u7 );
    q7 = ( u7  );

    //
    // Store 8 outputs to the output array with 2-samples stride.
    //

    // Q31 <- Q16.47 - 16 w/ rounding and saturation
    AE_S32RA64S_IP( q0, Y, +2*4 ); //  1
    AE_S32RA64S_IP( q1, Y, +2*4 ); //  3
    AE_S32RA64S_IP( q2, Y, +2*4 ); //  5
    AE_S32RA64S_IP( q3, Y, +2*4 ); //  7
    AE_S32RA64S_IP( q4, Y, +2*4 ); //  9
    AE_S32RA64S_IP( q5, Y, +2*4 ); // 11
    AE_S32RA64S_IP( q6, Y, +2*4 ); // 13
    AE_S32RA64S_IP( q7, Y, +1*4 ); // 15
  }

  return ( (int32_t*)D_wr );

} // firinterp2_proc()

/* Data processing function for a factor 3 interpolating FIR filter. */
static int32_t * firinterp3_proc( int32_t * restrict y,
                                  int32_t *          cbegin,
                                  int32_t *          cend,
                                  int32_t * restrict cpos,
                          const int32_t * restrict x,
                          const int16_t * restrict h,
                          int D, int N, int M )
{
        ae_int32x2 * restrict D_wr;
  const ae_int32x2 *          D_rd;
  const ae_int32x2 *          X;
        ae_f32     * restrict Y;
  const ae_f16x4   *          C;

  ae_f32x2   d0, d1, d2, d3, d4, d5;
  ae_int32x2 t0, t1, t2, t3, t4, t5;
  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_int64   u0, u1, u2, u3, u4, u5, u6, u7;
  ae_int64   v0, v1, v2, v3, v4, v5, v6, v7;
  ae_f16x4   c;

  int n, m;

  NASSERT( y && cbegin && cend && cpos && x && h 
          &&
          D == 3 && N > 0 && M > 0 );

  NASSERT( !( N & 7 ) && !( M & 3 ) );

  NASSERT( IS_ALIGN( cbegin ) &&
          IS_ALIGN( cend   ) &&
          IS_ALIGN( cpos   ) &&
          IS_ALIGN( x      ) &&
          IS_ALIGN( h      ) );

  //
  // Setup pointers and circular delay line buffer.
  //

  D_wr = (      ae_int32x2 *)cpos;
  X    = (const ae_int32x2 *)x;
  Y    = (      ae_f32     *)y;

  WUR_AE_CBEGIN0( (uintptr_t)cbegin );
  WUR_AE_CEND0  ( (uintptr_t)cend   );

  //
  // Break input signal into 8-samples blocks. For each block, store 8 samples
  // to the delay line buffer and compute 16 samples of interpolating FIR filter
  // response.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 8 input samples.
    // Q31
    AE_L32X2_IP( t0, X, +8 );
    AE_L32X2_IP( t1, X, +8 );
    AE_L32X2_IP( t2, X, +8 );
    AE_L32X2_IP( t3, X, +8 );

    // Store input samples into the circular delay line.
    // Q31
    AE_S32X2_XC( t0, D_wr, +8 );
    AE_S32X2_XC( t1, D_wr, +8 );
    AE_S32X2_XC( t2, D_wr, +8 );
    AE_S32X2_XC( t3, D_wr, +8 );

    //--------------------------------------------------------------------------
    // Bank 0
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q31
    AE_L32X2_XC( t0, D_rd, +8 );
    AE_L32X2_XC( t1, D_rd, +8 );
    AE_L32X2_XC( t2, D_rd, +8 );
    AE_L32X2_XC( t3, D_rd, +8 );

    d0 = ( t0 );
    d1 = ( t1 );
    d2 = ( t2 );
    d3 = ( t3 );

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q31
    AE_L32X2_XC( t4, D_rd, +8 );
    AE_L32X2_XC( t5, D_rd, +8 );

    d4 = ( t4 );
    d5 = ( t5 );

    // Reset the coefficients pointer. Now it looks at the oldest sample
    // coefficient, bank #0.
    C = (const ae_f16x4 *)h;

    // Load first 4 tap coefficients of bank #0.
    // Q15
    ae_f16x4_loadip( c, C, +8 );

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
    AE_MULFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
    AE_MULFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
    AE_MULFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
    AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
    AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
    AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3; d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //
    NASSERT(M>=4);
    __Pragma("loop_count min=1")
    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q31
      AE_L32X2_XC( t4, D_rd, +8 );
      AE_L32X2_XC( t5, D_rd, +8 );

      d4 = ( t4 );
      d5 = ( t5 );

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip( c, C, +8 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
      AE_MULAFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
      AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    //
    // Scale outputs by the interpolation factor 3.
    //

    // Q16.47 <- Q16.47 + Q16.47 + Q16.47 w/ saturation
    u0 = ( q0 );
    v0 = AE_ADD64S( u0, u0 );
    v0 = AE_ADD64S( u0, v0 );
    q0 = ( v0 );

    u1 = ( q1 );
    v1 = AE_ADD64S( u1, u1 );
    v1 = AE_ADD64S( u1, v1 );
    q1 = ( v1 );

    u2 = ( q2 );
    v2 = AE_ADD64S( u2, u2 );
    v2 = AE_ADD64S( u2, v2 );
    q2 = ( v2 );

    u3 = ( q3 );
    v3 = AE_ADD64S( u3, u3 );
    v3 = AE_ADD64S( u3, v3 );
    q3 = ( v3 );

    u4 = ( q4 );
    v4 = AE_ADD64S( u4, u4 );
    v4 = AE_ADD64S( u4, v4 );
    q4 = ( v4 );

    u5 = ( q5 );
    v5 = AE_ADD64S( u5, u5 );
    v5 = AE_ADD64S( u5, v5 );
    q5 = ( v5 );

    u6 = ( q6 );
    v6 = AE_ADD64S( u6, u6 );
    v6 = AE_ADD64S( u6, v6 );
    q6 = ( v6 );

    u7 = ( q7 );
    v7 = AE_ADD64S( u7, u7 );
    v7 = AE_ADD64S( u7, v7 );
    q7 = ( v7 );

    //
    // Store 8 outputs to the output array with 3-samples stride.
    //

    // Q31 <- Q16.47 - 16 w/ rounding and saturation
    AE_S32RA64S_IP( q0, Y,  +3*4 ); //  0
    AE_S32RA64S_IP( q1, Y,  +3*4 ); //  3
    AE_S32RA64S_IP( q2, Y,  +3*4 ); //  6
    AE_S32RA64S_IP( q3, Y,  +3*4 ); //  9
    AE_S32RA64S_IP( q4, Y,  +3*4 ); // 12
    AE_S32RA64S_IP( q5, Y,  +3*4 ); // 15
    AE_S32RA64S_IP( q6, Y,  +3*4 ); // 18
    AE_S32RA64S_XP( q7, Y, -20*4 ); // 21

    //--------------------------------------------------------------------------
    // Bank 1
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q31
    AE_L32X2_XC( t0, D_rd, +8 );
    AE_L32X2_XC( t1, D_rd, +8 );
    AE_L32X2_XC( t2, D_rd, +8 );
    AE_L32X2_XC( t3, D_rd, +8 );

    d0 = ( t0 );
    d1 = ( t1 );
    d2 = ( t2 );
    d3 = ( t3 );

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q31
    AE_L32X2_XC( t4, D_rd, +8 );
    AE_L32X2_XC( t5, D_rd, +8 );

    d4 = ( t4 );
    d5 = ( t5 );

    // Load first 4 tap coefficients of bank #1.
    // Q15
    ae_f16x4_loadip( c, C, +8 );

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
    AE_MULFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
    AE_MULFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
    AE_MULFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
    AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
    AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
    AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3; d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //
    NASSERT(M>=4);
    __Pragma("loop_count min=1")
    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q31
      AE_L32X2_XC( t4, D_rd, +8 );
      AE_L32X2_XC( t5, D_rd, +8 );

      d4 = ( t4 );
      d5 = ( t5 );

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip( c, C, +8 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
      AE_MULAFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
      AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    //
    // Scale outputs by the interpolation factor 3.
    //

    // Q16.47 <- Q16.47 + Q16.47 + Q16.47 w/ saturation
    u0 = ( q0 );
    v0 = AE_ADD64S( u0, u0 );
    v0 = AE_ADD64S( u0, v0 );
    q0 = ( v0 );

    u1 = ( q1 );
    v1 = AE_ADD64S( u1, u1 );
    v1 = AE_ADD64S( u1, v1 );
    q1 = ( v1 );

    u2 = ( q2 );
    v2 = AE_ADD64S( u2, u2 );
    v2 = AE_ADD64S( u2, v2 );
    q2 = ( v2 );

    u3 = ( q3 );
    v3 = AE_ADD64S( u3, u3 );
    v3 = AE_ADD64S( u3, v3 );
    q3 = ( v3 );

    u4 = ( q4 );
    v4 = AE_ADD64S( u4, u4 );
    v4 = AE_ADD64S( u4, v4 );
    q4 = ( v4 );

    u5 = ( q5 );
    v5 = AE_ADD64S( u5, u5 );
    v5 = AE_ADD64S( u5, v5 );
    q5 = ( v5 );

    u6 = ( q6 );
    v6 = AE_ADD64S( u6, u6 );
    v6 = AE_ADD64S( u6, v6 );
    q6 = ( v6 );

    u7 = ( q7 );
    v7 = AE_ADD64S( u7, u7 );
    v7 = AE_ADD64S( u7, v7 );
    q7 = ( v7 );

    //
    // Store 8 outputs to the output array with 3-samples stride.
    //

    // Q31 <- Q16.47 - 16 w/ rounding and saturation
    AE_S32RA64S_IP( q0, Y,  +3*4 ); //  1
    AE_S32RA64S_IP( q1, Y,  +3*4 ); //  4
    AE_S32RA64S_IP( q2, Y,  +3*4 ); //  7
    AE_S32RA64S_IP( q3, Y,  +3*4 ); // 10
    AE_S32RA64S_IP( q4, Y,  +3*4 ); // 13
    AE_S32RA64S_IP( q5, Y,  +3*4 ); // 16
    AE_S32RA64S_IP( q6, Y,  +3*4 ); // 19
    AE_S32RA64S_XP( q7, Y, -20*4 ); // 22

    //--------------------------------------------------------------------------
    // Bank 2
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q31
    AE_L32X2_XC( t0, D_rd, +8 );
    AE_L32X2_XC( t1, D_rd, +8 );
    AE_L32X2_XC( t2, D_rd, +8 );
    AE_L32X2_XC( t3, D_rd, +8 );

    d0 = ( t0 );
    d1 = ( t1 );
    d2 = ( t2 );
    d3 = ( t3 );

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q31
    AE_L32X2_XC( t4, D_rd, +8 );
    AE_L32X2_XC( t5, D_rd, +8 );

    d4 = ( t4 );
    d5 = ( t5 );

    // Load first 4 tap coefficients of bank #2.
    // Q15
    ae_f16x4_loadip( c, C, +8 );

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
    AE_MULFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
    AE_MULFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
    AE_MULFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
    AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
    AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
    AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3; d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //
    NASSERT(M>=4);
    __Pragma("loop_count min=1")
    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q31
      AE_L32X2_XC( t4, D_rd, +8 );
      AE_L32X2_XC( t5, D_rd, +8 );

      d4 = ( t4 );
      d5 = ( t5 );

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip( c, C, +8 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
      AE_MULAFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
      AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    //
    // Scale outputs by the interpolation factor 3.
    //

    // Q16.47 <- Q16.47 + Q16.47 + Q16.47 w/ saturation
    u0 = ( q0 );
    v0 = AE_ADD64S( u0, u0 );
    v0 = AE_ADD64S( u0, v0 );
    q0 = ( v0 );

    u1 = ( q1 );
    v1 = AE_ADD64S( u1, u1 );
    v1 = AE_ADD64S( u1, v1 );
    q1 = ( v1 );

    u2 = ( q2 );
    v2 = AE_ADD64S( u2, u2 );
    v2 = AE_ADD64S( u2, v2 );
    q2 = ( v2 );

    u3 = ( q3 );
    v3 = AE_ADD64S( u3, u3 );
    v3 = AE_ADD64S( u3, v3 );
    q3 = ( v3 );

    u4 = ( q4 );
    v4 = AE_ADD64S( u4, u4 );
    v4 = AE_ADD64S( u4, v4 );
    q4 = ( v4 );

    u5 = ( q5 );
    v5 = AE_ADD64S( u5, u5 );
    v5 = AE_ADD64S( u5, v5 );
    q5 = ( v5 );

    u6 = ( q6 );
    v6 = AE_ADD64S( u6, u6 );
    v6 = AE_ADD64S( u6, v6 );
    q6 = ( v6 );

    u7 = ( q7 );
    v7 = AE_ADD64S( u7, u7 );
    v7 = AE_ADD64S( u7, v7 );
    q7 = ( v7 );

    //
    // Store 8 outputs to the output array with 3-samples stride.
    //

    // Q31 <- Q16.47 - 16 w/ rounding and saturation
    AE_S32RA64S_IP( q0, Y, +3*4 ); //  2
    AE_S32RA64S_IP( q1, Y, +3*4 ); //  5
    AE_S32RA64S_IP( q2, Y, +3*4 ); //  8
    AE_S32RA64S_IP( q3, Y, +3*4 ); // 11
    AE_S32RA64S_IP( q4, Y, +3*4 ); // 14
    AE_S32RA64S_IP( q5, Y, +3*4 ); // 17
    AE_S32RA64S_IP( q6, Y, +3*4 ); // 20
    AE_S32RA64S_IP( q7, Y, +1*4 ); // 23
  }

  return ( (int32_t*)D_wr );

} // firinterp3_proc()

/* Data processing function for a factor 4 interpolating FIR filter. */
static int32_t * firinterp4_proc( int32_t * restrict y,
                                  int32_t *          cbegin,
                                  int32_t *          cend,
                                  int32_t * restrict cpos,
                          const int32_t * restrict x,
                          const int16_t * restrict h,
                          int D, int N, int M )
{
        ae_int32x2 * restrict D_wr;
  const ae_int32x2 *          D_rd;
  const ae_int32x2 *          X;
        ae_f32     * restrict Y;
  const ae_f16x4   *          C;

  ae_f32x2   d0, d1, d2, d3, d4, d5;
  ae_int32x2 t0, t1, t2, t3, t4, t5;
  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_int64   u0, u1, u2, u3, u4, u5, u6, u7;
  ae_f16x4   c;

  int n, m;

  NASSERT( y && cbegin && cend && cpos && x && h 
          &&
          D == 4 && N > 0 && M > 0 );

  NASSERT( !( N & 7 ) && !( M & 3 ) );

  NASSERT( IS_ALIGN( cbegin ) &&
          IS_ALIGN( cend   ) &&
          IS_ALIGN( cpos   ) &&
          IS_ALIGN( x      ) &&
          IS_ALIGN( h      ) );

  //
  // Setup pointers and circular delay line buffer.
  //

  D_wr = (      ae_int32x2 *)cpos;
  X    = (const ae_int32x2 *)x;
  Y    = (      ae_f32     *)y;

  WUR_AE_CBEGIN0( (uintptr_t)cbegin );
  WUR_AE_CEND0  ( (uintptr_t)cend   );

  //
  // Break input signal into 8-samples blocks. For each block, store 8 samples
  // to the delay line buffer and compute 16 samples of interpolating FIR filter
  // response.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 8 input samples.
    // Q31
    AE_L32X2_IP( t0, X, +8 );
    AE_L32X2_IP( t1, X, +8 );
    AE_L32X2_IP( t2, X, +8 );
    AE_L32X2_IP( t3, X, +8 );

    // Store input samples into the circular delay line.
    // Q31
    AE_S32X2_XC( t0, D_wr, +8 );
    AE_S32X2_XC( t1, D_wr, +8 );
    AE_S32X2_XC( t2, D_wr, +8 );
    AE_S32X2_XC( t3, D_wr, +8 );

    //--------------------------------------------------------------------------
    // Bank 0
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q31
    AE_L32X2_XC( t0, D_rd, +8 );
    AE_L32X2_XC( t1, D_rd, +8 );
    AE_L32X2_XC( t2, D_rd, +8 );
    AE_L32X2_XC( t3, D_rd, +8 );

    d0 = ( t0 );
    d1 = ( t1 );
    d2 = ( t2 );
    d3 = ( t3 );

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q31
    AE_L32X2_XC( t4, D_rd, +8 );
    AE_L32X2_XC( t5, D_rd, +8 );

    d4 = ( t4 );
    d5 = ( t5 );

    // Reset the coefficients pointer. Now it looks at the oldest sample
    // coefficient, bank #0.
    C = (const ae_f16x4 *)h;

    // Load first 4 tap coefficients of bank #0.
    // Q15
    ae_f16x4_loadip( c, C, +8 );

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
    AE_MULFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
    AE_MULFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
    AE_MULFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
    AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
    AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
    AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3; d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //
    NASSERT(M>=4);
    __Pragma("loop_count min=1")
    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q31
      AE_L32X2_XC( t4, D_rd, +8 );
      AE_L32X2_XC( t5, D_rd, +8 );

      d4 = ( t4 );
      d5 = ( t5 );

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip( c, C, +8 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
      AE_MULAFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
      AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    //
    // Scale outputs by the interpolation factor 4.
    //

    // Q16.47 <- Q16.47*Q0
    u0 = ( q0 );
    u0 = AE_SLAI64S( u0, 2 );
    q0 = ( u0  );

    u1 = ( q1 );
    u1 = AE_SLAI64S( u1, 2 );
    q1 = ( u1  );

    u2 = ( q2 );
    u2 = AE_SLAI64S( u2, 2 );
    q2 = ( u2  );

    u3 = ( q3 );
    u3 = AE_SLAI64S( u3, 2 );
    q3 = ( u3  );

    u4 = ( q4 );
    u4 = AE_SLAI64S( u4, 2 );
    q4 = ( u4  );

    u5 = ( q5 );
    u5 = AE_SLAI64S( u5, 2 );
    q5 = ( u5  );

    u6 = ( q6 );
    u6 = AE_SLAI64S( u6, 2 );
    q6 = ( u6  );

    u7 = ( q7 );
    u7 = AE_SLAI64S( u7, 2 );
    q7 = ( u7  );

    //
    // Store 8 outputs to the output array with 4-samples stride.
    //

    // Q31 <- Q16.47 - 16 w/ rounding and saturation
    AE_S32RA64S_IP( q0, Y,  +4*4 ); //  0
    AE_S32RA64S_IP( q1, Y,  +4*4 ); //  4
    AE_S32RA64S_IP( q2, Y,  +4*4 ); //  8
    AE_S32RA64S_IP( q3, Y,  +4*4 ); // 12
    AE_S32RA64S_IP( q4, Y,  +4*4 ); // 16
    AE_S32RA64S_IP( q5, Y,  +4*4 ); // 20
    AE_S32RA64S_IP( q6, Y,  +4*4 ); // 24
    AE_S32RA64S_XP( q7, Y, -27*4 ); // 28

    //--------------------------------------------------------------------------
    // Bank 1
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q31
    AE_L32X2_XC( t0, D_rd, +8 );
    AE_L32X2_XC( t1, D_rd, +8 );
    AE_L32X2_XC( t2, D_rd, +8 );
    AE_L32X2_XC( t3, D_rd, +8 );

    d0 = ( t0 );
    d1 = ( t1 );
    d2 = ( t2 );
    d3 = ( t3 );

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q31
    AE_L32X2_XC( t4, D_rd, +8 );
    AE_L32X2_XC( t5, D_rd, +8 );

    d4 = ( t4 );
    d5 = ( t5 );

    // Load first 4 tap coefficients of bank #1.
    // Q15
    ae_f16x4_loadip( c, C, +8 );

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
    AE_MULFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
    AE_MULFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
    AE_MULFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
    AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
    AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
    AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3; d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //
    NASSERT(M>=4);
    __Pragma("loop_count min=1")
    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q31
      AE_L32X2_XC( t4, D_rd, +8 );
      AE_L32X2_XC( t5, D_rd, +8 );

      d4 = ( t4 );
      d5 = ( t5 );

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip( c, C, +8 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
      AE_MULAFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
      AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    //
    // Scale outputs by the interpolation factor 4.
    //

    // Q16.47 <- Q16.47*Q0
    u0 = ( q0 );
    u0 = AE_SLAI64S( u0, 2 );
    q0 = ( u0  );

    u1 = ( q1 );
    u1 = AE_SLAI64S( u1, 2 );
    q1 = ( u1  );

    u2 = ( q2 );
    u2 = AE_SLAI64S( u2, 2 );
    q2 = ( u2  );

    u3 = ( q3 );
    u3 = AE_SLAI64S( u3, 2 );
    q3 = ( u3  );

    u4 = ( q4 );
    u4 = AE_SLAI64S( u4, 2 );
    q4 = ( u4  );

    u5 = ( q5 );
    u5 = AE_SLAI64S( u5, 2 );
    q5 = ( u5  );

    u6 = ( q6 );
    u6 = AE_SLAI64S( u6, 2 );
    q6 = ( u6  );

    u7 = ( q7 );
    u7 = AE_SLAI64S( u7, 2 );
    q7 = ( u7  );

    //
    // Store 8 outputs to the output array with 4-samples stride.
    //

    // Q31 <- Q16.47 - 16 w/ rounding and saturation
    AE_S32RA64S_IP( q0, Y,  +4*4 ); //  1
    AE_S32RA64S_IP( q1, Y,  +4*4 ); //  5
    AE_S32RA64S_IP( q2, Y,  +4*4 ); //  9
    AE_S32RA64S_IP( q3, Y,  +4*4 ); // 13
    AE_S32RA64S_IP( q4, Y,  +4*4 ); // 17
    AE_S32RA64S_IP( q5, Y,  +4*4 ); // 21
    AE_S32RA64S_IP( q6, Y,  +4*4 ); // 25
    AE_S32RA64S_XP( q7, Y, -27*4 ); // 29

    //--------------------------------------------------------------------------
    // Bank 2
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q31
    AE_L32X2_XC( t0, D_rd, +8 );
    AE_L32X2_XC( t1, D_rd, +8 );
    AE_L32X2_XC( t2, D_rd, +8 );
    AE_L32X2_XC( t3, D_rd, +8 );

    d0 = ( t0 );
    d1 = ( t1 );
    d2 = ( t2 );
    d3 = ( t3 );

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q31
    AE_L32X2_XC( t4, D_rd, +8 );
    AE_L32X2_XC( t5, D_rd, +8 );

    d4 = ( t4 );
    d5 = ( t5 );

    // Load first 4 tap coefficients of bank #2.
    // Q15
    ae_f16x4_loadip( c, C, +8 );

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
    AE_MULFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
    AE_MULFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
    AE_MULFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
    AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
    AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
    AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3; d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //
    NASSERT(M>=4);
    __Pragma("loop_count min=1")
    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q31
      AE_L32X2_XC( t4, D_rd, +8 );
      AE_L32X2_XC( t5, D_rd, +8 );

      d4 = ( t4 );
      d5 = ( t5 );

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip( c, C, +8 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
      AE_MULAFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
      AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    //
    // Scale outputs by the interpolation factor 4.
    //

    // Q16.47 <- Q16.47*Q0
    u0 = ( q0 );
    u0 = AE_SLAI64S( u0, 2 );
    q0 = ( u0  );

    u1 = ( q1 );
    u1 = AE_SLAI64S( u1, 2 );
    q1 = ( u1  );

    u2 = ( q2 );
    u2 = AE_SLAI64S( u2, 2 );
    q2 = ( u2  );

    u3 = ( q3 );
    u3 = AE_SLAI64S( u3, 2 );
    q3 = ( u3  );

    u4 = ( q4 );
    u4 = AE_SLAI64S( u4, 2 );
    q4 = ( u4  );

    u5 = ( q5 );
    u5 = AE_SLAI64S( u5, 2 );
    q5 = ( u5  );

    u6 = ( q6 );
    u6 = AE_SLAI64S( u6, 2 );
    q6 = ( u6  );

    u7 = ( q7 );
    u7 = AE_SLAI64S( u7, 2 );
    q7 = ( u7  );

    //
    // Store 8 outputs to the output array with 4-samples stride.
    //

    // Q31 <- Q16.47 - 16 w/ rounding and saturation
    AE_S32RA64S_IP( q0, Y,  +4*4 ); //  2
    AE_S32RA64S_IP( q1, Y,  +4*4 ); //  6
    AE_S32RA64S_IP( q2, Y,  +4*4 ); // 10
    AE_S32RA64S_IP( q3, Y,  +4*4 ); // 14
    AE_S32RA64S_IP( q4, Y,  +4*4 ); // 18
    AE_S32RA64S_IP( q5, Y,  +4*4 ); // 22
    AE_S32RA64S_IP( q6, Y,  +4*4 ); // 26
    AE_S32RA64S_XP( q7, Y, -27*4 ); // 30

    //--------------------------------------------------------------------------
    // Bank 3
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q31
    AE_L32X2_XC( t0, D_rd, +8 );
    AE_L32X2_XC( t1, D_rd, +8 );
    AE_L32X2_XC( t2, D_rd, +8 );
    AE_L32X2_XC( t3, D_rd, +8 );

    d0 = ( t0 );
    d1 = ( t1 );
    d2 = ( t2 );
    d3 = ( t3 );

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q31
    AE_L32X2_XC( t4, D_rd, +8 );
    AE_L32X2_XC( t5, D_rd, +8 );

    d4 = ( t4 );
    d5 = ( t5 );

    // Load first 4 tap coefficients of bank #3.
    // Q15
    ae_f16x4_loadip( c, C, +8 );

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
    AE_MULFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
    AE_MULFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
    AE_MULFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
    AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
    AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
    AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3; d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //
    NASSERT(M>=4);
    __Pragma("loop_count min=1")
    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q31
      AE_L32X2_XC( t4, D_rd, +8 );
      AE_L32X2_XC( t5, D_rd, +8 );

      d4 = ( t4 );
      d5 = ( t5 );

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip( c, C, +8 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
      AE_MULAFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
      AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    //
    // Scale outputs by the interpolation factor 4.
    //

    // Q16.47 <- Q16.47*Q0
    u0 = ( q0 );
    u0 = AE_SLAI64S( u0, 2 );
    q0 = ( u0  );

    u1 = ( q1 );
    u1 = AE_SLAI64S( u1, 2 );
    q1 = ( u1  );

    u2 = ( q2 );
    u2 = AE_SLAI64S( u2, 2 );
    q2 = ( u2  );

    u3 = ( q3 );
    u3 = AE_SLAI64S( u3, 2 );
    q3 = ( u3  );

    u4 = ( q4 );
    u4 = AE_SLAI64S( u4, 2 );
    q4 = ( u4  );

    u5 = ( q5 );
    u5 = AE_SLAI64S( u5, 2 );
    q5 = ( u5  );

    u6 = ( q6 );
    u6 = AE_SLAI64S( u6, 2 );
    q6 = ( u6  );

    u7 = ( q7 );
    u7 = AE_SLAI64S( u7, 2 );
    q7 = ( u7  );

    //
    // Store 8 outputs to the output array with 4-samples stride.
    //

    // Q31 <- Q16.47 - 16 w/ rounding and saturation
    AE_S32RA64S_IP( q0, Y, +4*4 ); //  3
    AE_S32RA64S_IP( q1, Y, +4*4 ); //  7
    AE_S32RA64S_IP( q2, Y, +4*4 ); // 11
    AE_S32RA64S_IP( q3, Y, +4*4 ); // 15
    AE_S32RA64S_IP( q4, Y, +4*4 ); // 19
    AE_S32RA64S_IP( q5, Y, +4*4 ); // 23
    AE_S32RA64S_IP( q6, Y, +4*4 ); // 27
    AE_S32RA64S_IP( q7, Y, +1*4 ); // 31
  }

  return ( (int32_t*)D_wr );

} // firinterp4_proc()

/* Data processing function for a generic interpolating FIR filter with an
 * arbitrary interpolation factor. */
static int32_t * firinterpX_proc( int32_t * restrict y,
                                  int32_t *          cbegin,
                                  int32_t *          cend,
                                  int32_t * restrict cpos,
                          const int32_t * restrict x,
                          const int16_t * restrict h,
                          int D, int N, int M )
{
        ae_int32x2 * restrict D_wr;
  const ae_int32x2 *          D_rd;
  const ae_int32x2 *          X;
        ae_f32     * restrict Y;
  const ae_f16x4   *          C;

  ae_f32x2   d0, d1, d2, d3, d4, d5;
  ae_int32x2 t0, t1, t2, t3, t4, t5;
  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_int64   u0, u1, u2, u3, u4, u5, u6, u7;
  ae_f16x4   c;

  int        g_exp;
  ae_p24x2s g_frac;

  int d, n, m;

  NASSERT( y && cbegin && cend && cpos && x && h 
          &&
          D > 1 && N > 0 && M > 0 );

  NASSERT( !( N & 7 ) && !( M & 3 ) );

  NASSERT( IS_ALIGN( cbegin ) &&
          IS_ALIGN( cend   ) &&
          IS_ALIGN( cpos   ) &&
          IS_ALIGN( x      ) &&
          IS_ALIGN( h      ) );

  //
  // Calculate the scaling exponent and fraction.
  //

  // Q(15+8-g_exp) <- Q0 + 15 + 8 - g_exp
  g_exp  = 31 - AE_NSAZ32_L( D );
  {
  int q_frac_int =  D << ( 15 + 8 - g_exp );
  g_frac = *(ae_p24s *) &q_frac_int;
  }

  WUR_AE_SAR( g_exp );

  //
  // Setup pointers and circular delay line buffer.
  //

  D_wr = (      ae_int32x2 *)cpos;
  X    = (const ae_int32x2 *)x;
  Y    = (      ae_f32     *)y;

  WUR_AE_CBEGIN0( (uintptr_t)cbegin );
  WUR_AE_CEND0  ( (uintptr_t)cend   );

  //
  // Break input signal into 8-samples blocks. For each block, store 8 samples
  // to the delay line buffer and compute 8 outputs per filter bank,
  // in all 8*D samples
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 8 input samples.
    // Q31
    AE_L32X2_IP( t0, X, +8 );
    AE_L32X2_IP( t1, X, +8 );
    AE_L32X2_IP( t2, X, +8 );
    AE_L32X2_IP( t3, X, +8 );

    // Store input samples into the circular delay line.
    // Q31
    AE_S32X2_XC( t0, D_wr, +8 );
    AE_S32X2_XC( t1, D_wr, +8 );
    AE_S32X2_XC( t2, D_wr, +8 );
    AE_S32X2_XC( t3, D_wr, +8 );

    // Reset the coefficients pointer. Now it looks at the oldest sample
    // coefficient, bank #0.
    C = (const ae_f16x4 *)h;

    // Loop through the filter banks #0..#D-1
    for ( d=0; d<D; d++ )
    {
      // Start reading the delay line from the oldest sample, M+8 samples back
      // from the newest one.
      D_rd = D_wr;

      // Load 8 oldest samples.
      // Q31
      AE_L32X2_XC( t0, D_rd, +8 );
      AE_L32X2_XC( t1, D_rd, +8 );
      AE_L32X2_XC( t2, D_rd, +8 );
      AE_L32X2_XC( t3, D_rd, +8 );

      d0 = ( t0 );
      d1 = ( t1 );
      d2 = ( t2 );
      d3 = ( t3 );

      //
      // Inner loop prologue: process first 4 taps for 8 accumulators.
      //

      // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
      // registers to perform 4 taps processing for 8 accumulators.
      // Q31
      AE_L32X2_XC( t4, D_rd, +8 );
      AE_L32X2_XC( t5, D_rd, +8 );

      d4 = ( t4 );
      d5 = ( t5 );

      // Load first 4 tap coefficients of bank #d.
      // Q15
      ae_f16x4_loadip( c, C, +8 );

      // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
      AE_MULFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
      AE_MULFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
      AE_MULFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
      AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
      AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;

      //
      // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
      // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
      // coefficients appended to the impulse response to avoid a 1-sample
      // response delay.
      //
      NASSERT(M>=4);
      __Pragma("loop_count min=1")
      for ( m=0; m<(M>>2); m++ )
      {
        // Load 4 samples.
        // Q31
        AE_L32X2_XC( t4, D_rd, +8 );
        AE_L32X2_XC( t5, D_rd, +8 );

        d4 = ( t4 );
        d5 = ( t5 );

        // Load 4 tap coefficients.
        // Q15
        ae_f16x4_loadip( c, C, +8 );

        // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
        AE_MULAFD32X16X2_FIR_HH( q0, q1, d0, d1, c );
        AE_MULAFD32X16X2_FIR_HH( q2, q3, d1, d2, c );
        AE_MULAFD32X16X2_FIR_HH( q4, q5, d2, d3, c );
        AE_MULAFD32X16X2_FIR_HH( q6, q7, d3, d4, c );

        // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
        AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, c );
        AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, c );
        AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, c );
        AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, c );

        // 4 taps are done. Move 4 samples out of the registers file.
        d0 = d2; d1 = d3; d2 = d4; d3 = d5;
      }

      //
      // Scale outputs by the interpolation factor D.
      //

      u0 = ( q0 );
      // Q(31+16) <- sat32( Q16.47 - 16 ) + 16 w/ rounding
      u0 = AE_ROUNDSQ32F48ASYM( u0 );
      // Q(16.47-g_exp) <- [ Q(31+16) - 16 ] * [ Q(15+8-g_exp) - 8 ] + 1
      u0 = AE_MULF48Q32SP16S_L( u0, g_frac );
      // Q16.47 <- Q(16.47-g_exp) + g_exp
      u0 = AE_SLAS64S( u0 );
      q0 = ( u0 );

      u1 = ( q1 );
      // Q(31+16) <- sat32( Q16.47 - 16 ) + 16 w/ rounding
      u1 = AE_ROUNDSQ32F48ASYM( u1 );
      // Q(16.47-g_exp) <- [ Q(31+16) - 16 ] * [ Q(15+8-g_exp) - 8 ] + 1
      u1 = AE_MULF48Q32SP16S_L( u1, g_frac );
      // Q16.47 <- Q(16.47-g_exp) + g_exp
      u1 = AE_SLAS64S( u1 );
      q1 = ( u1 );

      u2 = ( q2 );
      // Q(31+16) <- sat32( Q16.47 - 16 ) + 16 w/ rounding
      u2 = AE_ROUNDSQ32F48ASYM( u2 );
      // Q(16.47-g_exp) <- [ Q(31+16) - 16 ] * [ Q(15+8-g_exp) - 8 ] + 1
      u2 = AE_MULF48Q32SP16S_L( u2, g_frac );
      // Q16.47 <- Q(16.47-g_exp) + g_exp
      u2 = AE_SLAS64S( u2 );
      q2 = ( u2 );

      u3 = ( q3 );
      // Q(31+16) <- sat32( Q16.47 - 16 ) + 16 w/ rounding
      u3 = AE_ROUNDSQ32F48ASYM( u3 );
      // Q(16.47-g_exp) <- [ Q(31+16) - 16 ] * [ Q(15+8-g_exp) - 8 ] + 1
      u3 = AE_MULF48Q32SP16S_L( u3, g_frac );
      // Q16.47 <- Q(16.47-g_exp) + g_exp
      u3 = AE_SLAS64S( u3 );
      q3 = ( u3 );

      u4 = ( q4 );
      // Q(31+16) <- sat32( Q16.47 - 16 ) + 16 w/ rounding
      u4 = AE_ROUNDSQ32F48ASYM( u4 );
      // Q(16.47-g_exp) <- [ Q(31+16) - 16 ] * [ Q(15+8-g_exp) - 8 ] + 1
      u4 = AE_MULF48Q32SP16S_L( u4, g_frac );
      // Q16.47 <- Q(16.47-g_exp) + g_exp
      u4 = AE_SLAS64S( u4 );
      q4 = ( u4 );

      u5 = ( q5 );
      // Q(31+16) <- sat32( Q16.47 - 16 ) + 16 w/ rounding
      u5 = AE_ROUNDSQ32F48ASYM( u5 );
      // Q(16.47-g_exp) <- [ Q(31+16) - 16 ] * [ Q(15+8-g_exp) - 8 ] + 1
      u5 = AE_MULF48Q32SP16S_L( u5, g_frac );
      // Q16.47 <- Q(16.47-g_exp) + g_exp
      u5 = AE_SLAS64S( u5 );
      q5 = ( u5 );

      u6 = ( q6 );
      // Q(31+16) <- sat32( Q16.47 - 16 ) + 16 w/ rounding
      u6 = AE_ROUNDSQ32F48ASYM( u6 );
      // Q(16.47-g_exp) <- [ Q(31+16) - 16 ] * [ Q(15+8-g_exp) - 8 ] + 1
      u6 = AE_MULF48Q32SP16S_L( u6, g_frac );
      // Q16.47 <- Q(16.47-g_exp) + g_exp
      u6 = AE_SLAS64S( u6 );
      q6 = ( u6 );

      u7 = ( q7 );
      // Q(31+16) <- sat32( Q16.47 - 16 ) + 16 w/ rounding
      u7 = AE_ROUNDSQ32F48ASYM( u7 );
      // Q(16.47-g_exp) <- [ Q(31+16) - 16 ] * [ Q(15+8-g_exp) - 8 ] + 1
      u7 = AE_MULF48Q32SP16S_L( u7, g_frac );
      // Q16.47 <- Q(16.47-g_exp) + g_exp
      u7 = AE_SLAS64S( u7 );
      q7 = ( u7 );

      //
      // Store 8 outputs to the output array with D-samples stride.
      //

      // Q31 <- Q16.47 - 16 w/ rounding and saturation
      AE_S32RA64S_XP( q0, Y, D*4 );       // d+0*D
      AE_S32RA64S_XP( q1, Y, D*4 );       // d+1*D
      AE_S32RA64S_XP( q2, Y, D*4 );       // d+2*D
      AE_S32RA64S_XP( q3, Y, D*4 );       // d+3*D
      AE_S32RA64S_XP( q4, Y, D*4 );       // d+4*D
      AE_S32RA64S_XP( q5, Y, D*4 );       // d+5*D
      AE_S32RA64S_XP( q6, Y, D*4 );       // d+6*D
      AE_S32RA64S_XP( q7, Y, (1-7*D)*4 ); // d+7*D
    }

    Y = (ae_f32*)( (uintptr_t)Y + 7*D*4 );
  }

  return ( (int32_t*)D_wr );

} // firinterpX_proc()
