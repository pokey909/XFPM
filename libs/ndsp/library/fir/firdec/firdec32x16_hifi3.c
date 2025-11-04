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
    Decimating block real FIR filter, 32x16-bit
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Decimating Block Real FIR Filter
  Computes a real FIR filter (direct-form) with decimation using IR stored 
  in vector h. The real data input is stored in vector x. The filter output 
  result is stored in vector r. The filter calculates N output samples using
  M coefficients and requires last D*N+M-1 samples on the delay line.
  NOTE:
  - To avoid aliasing IR should be synthesized in such a way to be narrower 
    than input sample rate divided to 2D.
  - user application is not responsible for management of delay lines

  Precision: 
  16x16     16-bit data, 16-bit coefficients, 16-bit outputs
  24x24     24-bit data, 24-bit coefficients, 24-bit outputs
  32x16     32-bit data, 16-bit coefficients, 32-bit outputs
  32x32     32-bit data, 32-bit coefficients, 32-bit outputs
  f         floating point

  Input:
  h[M]      filter coefficients; h[0] is to be multiplied with the newest 
            sample, Q15, Q31, floating point
  D         decimation factor 
  N         length of output sample block
  M         length of filter
  x[D*N]    input samples, Q15, Q31, floating point
  Output:
  y[N]      output samples, Q15, Q31, floating point

  Restriction:
  x,h,r     should not overlap
  x,h       aligned on an 8-bytes boundary
  N         multiple of 8
  D         should exceed 1

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
#define MAGIC     0xed1722f0

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
      ( (size_t)(size) + (align) - 1 )

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
      (void*)( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

#define sz_i16    sizeof(int16_t)
#define sz_i32    sizeof(int32_t)

/*
 * Data processing function of a particular decimating filter. Stores a
 * block of input samples to the circular delay line buffer and computes
 * decimating FIR filter's response.
 * Input:
 *   cbegin  - circular delay line buffer start address
 *   cend    - points at the next to the last byte of the circular buffer
 *   cpos    - next position in the buffer to be filled with an input sample
 *   x[N*D]  - input samples
 *   h[]     - decimating FIR filter coefficients, array layout varies
 * Output:
 *   y[N]    - output samples
 *   retval  - updated circular buffer's next position pointer
 * Notes and restrictions:
 *   1. Most of data processing functions feature a single, hard-coded 
 *      decimation factor, so they expect a determined value for parameter D.
 *   2. All pointers with the exception of y[N] must be aligned on an 8-bytes
 *      boundary.
 *   3. N - must be a multiple of 8.
 *   4. M - must be a multiple of 4. 
*/
typedef int32_t * (proc_fxn_t)( int32_t * restrict y,
                                int32_t *          cbegin,
                                int32_t *          cend,
                                int32_t * restrict cpos,
                          const int32_t * restrict x,
                          const int16_t * restrict h,
                          int D, int N, int M );

static proc_fxn_t firdec2_proc;
static proc_fxn_t firdec3_proc;
static proc_fxn_t firdec4_proc;
static proc_fxn_t firdecX_proc;

/* Decimator instance structure. */
typedef struct tag_firdec32x16_t
{
  uint32_t        magic;     // Instance pointer validation number
  int             D;         // Decimation factor
  int             M;         // Number of filter coefficients
  proc_fxn_t *    procFxn;   // Filter data processing function
  const int16_t * coef;      // Filter coefficients
  int32_t *       delayLine; // Delay line
  int             delayLen;  // Delay line length, in samples
  int32_t *       delayPos;  // Delay line slot to be filled next

} firdec32x16_t, *firdec32x16_ptr_t;

/* Calculate the memory block size for a decimator with given attributes. */
size_t firdec32x16_alloc( int D, int M )
{
  int bM, P;
  int delayLen, coefNum;

  ASSERT( (D > 1) &&  (M > 0) );

  if ( D == 2 )
  {
    P  = ( M + 1 )/2;    // Coefficients per bank
    bM = P + ( -P & 3 ); // Subfilter order

    coefNum  = 2*bM + 4;
    delayLen = 2*( bM + 8 );
  }
  else if ( D == 3 )
  {
    P  = ( M + 2 )/3;    // Coefficients per bank
    bM = P + ( -P & 3 ); // Subfilter order

    coefNum  = 3*bM + 4;
    delayLen = 3*( bM + 8 );
  }
  else if ( D == 4 )
  {
    P  = ( M + 3 ) / 4;  // Coefficients per bank
    bM = P + ( -P & 3 ); // Subfilter order

    coefNum  = 4*bM + 4;
    delayLen = 4*( bM + 8 );
  }
  else
  {
    bM = M + ( -M & 3 );

    coefNum  = bM + 4;
    delayLen = bM + 8*D;
  }

  return ( ALIGNED_SIZE( sizeof( firdec32x16_t ), 4 ) 
           + // Delay line
           ALIGNED_SIZE( delayLen*sz_i32, 8 )
           + // Coefficients
           ALIGNED_SIZE( coefNum*sz_i16, 8 ) );

} /* firdec32x16_alloc() */

/* Initialize the decimator structure. The delay line is zeroed. */
firdec32x16_handle_t firdec32x16_init( void * objmem, int D, 
                                       int M, const int16_t * restrict h )
{
  firdec32x16_ptr_t firdec;
  void *            ptr;
  int16_t *         coefBuf;
  int               coefNum;
  int32_t *         delayLine;
  int               delayLen;
  proc_fxn_t *      procFxn;

  int16_t * coef_a;
  int16_t * coef_b;
  int16_t * coef_c;
  int16_t * coef_d;

  int bM, P=0;
  int m;

  ASSERT( objmem && (D > 1) && (M > 0) && h );

  ASSERT( IS_ALIGN( h ) );

  //
  // Select the processing function, delay line length and coefficients
  // block layout.
  //

  if ( D == 2 )
  {
    P  = ( M + 1 )/2;    // Coefficients per bank
    bM = P + ( -P & 3 ); // Subfilter order

    coefNum  = 2*bM + 4;
    delayLen = 2*( bM + 8 );
    procFxn  = &firdec2_proc;
  }
  else if ( D == 3 )
  {
    P  = ( M + 2 )/3;    // Coefficients per bank
    bM = P + ( -P & 3 ); // Subfilter order

    coefNum  = 3*bM + 4;
    delayLen = 3*( bM + 8 );
    procFxn  = &firdec3_proc;
  }
  else if ( D == 4 )
  {
    P  = ( M + 3 ) / 4;  // Coefficients per bank
    bM = P + ( -P & 3 ); // Subfilter order

    coefNum  = 4*bM + 4;
    delayLen = 4*( bM + 8 );
    procFxn  = &firdec4_proc;
  }
  else
  {
    bM = M + ( -M & 3 );

    coefNum  = bM + 4;
    delayLen = bM + 8*D;
    procFxn  = &firdecX_proc;
  }

  //
  // Partition the memory block.
  //

  ptr       = objmem;
  firdec    = (firdec32x16_ptr_t)(ALIGNED_ADDR( ptr, 4 ));
  ptr       = firdec + 1;
  delayLine = (int32_t*)(ALIGNED_ADDR( ptr, 8 ));
  ptr       = delayLine + delayLen;
  coefBuf   = (int16_t*)(ALIGNED_ADDR( ptr, 8 ));
  ptr       = coefBuf + coefNum;

  ASSERT( (int8_t*)ptr - (int8_t*)objmem <= (int)firdec32x16_alloc( D, M ) );

  //
  // Copy the filter coefficients.
  //

  if ( D == 2 )
  {
    // Partition the coefficients area into two banks.
    coef_a = coefBuf + 0*bM;
    coef_b = coefBuf + 1*bM;

    // Insert zeros to make subfilter order a multiple of 4.
    for ( m=0; m<(bM-P); m++ )
    {
      *coef_a++ = 0;
      *coef_b++ = 0;
    }

    // For bank #1, pad the coefficients with 4 zeros to avoid a 1-sample delay.
    coef_b[0] = coef_b[P+1] = coef_b[P+2] = coef_b[P+3] = 0;

    // Align on coefficients pairs.
    switch ( M - 2*(P-1) )
    {
    case 1: 
      coef_a[0] = 0;
      coef_b[1] = h[2*(P-1)+0];
      break;
    case 2:
      coef_a[0] = h[2*(P-1)+1];
      coef_b[1] = h[2*(P-1)+0];
    }

    // Copy coefficient pairs in reverted order.
    for ( m=1; m<P; m++ )
    {
      coef_a[m  ] = h[2*(P-1-m)+1];
      coef_b[m+1] = h[2*(P-1-m)+0];
    }
  }
  else if ( D == 3 )
  {
    // Partition the coefficients area into three banks.
    coef_a = coefBuf + 0*bM;
    coef_b = coefBuf + 1*bM;
    coef_c = coefBuf + 2*bM;

    // Insert zeros to make subfilter order a multiple of 4.
    for ( m=0; m<(bM-P); m++ )
    {
      *coef_a++ = 0;
      *coef_b++ = 0;
      *coef_c++ = 0;
    }

    // For bank #2, pad the coefficients with 4 zeros to avoid a 1-sample delay.
    coef_c[0] = coef_c[P+1] = coef_c[P+2] = coef_c[P+3] = 0;

    // Align on coefficients triples.
    switch ( M - 3*(P-1) )
    {
    case 1: 
      coef_a[0] = 0;
      coef_b[0] = 0;
      coef_c[1] = h[3*(P-1)+0];
      break;
    case 2:
      coef_a[0] = 0;
      coef_b[0] = h[3*(P-1)+1];
      coef_c[1] = h[3*(P-1)+0];
      break;
    case 3:
      coef_a[0] = h[3*(P-1)+2];
      coef_b[0] = h[3*(P-1)+1];
      coef_c[1] = h[3*(P-1)+0];
    }

    // Copy coefficient triples in reverted order.
    for ( m=1; m<P; m++ )
    {
      coef_a[m  ] = h[3*(P-1-m)+2];
      coef_b[m  ] = h[3*(P-1-m)+1];
      coef_c[m+1] = h[3*(P-1-m)+0];
    }
  }
  else if ( D == 4 )
  {
    // Partition the coefficients area into four banks.
    coef_a = coefBuf + 0*bM;
    coef_b = coefBuf + 1*bM;
    coef_c = coefBuf + 2*bM;
    coef_d = coefBuf + 3*bM;

    // Insert zeros to make subfilter order a multiple of 4.
    for ( m=0; m<(bM-P); m++ )
    {
      *coef_a++ = 0;
      *coef_b++ = 0;
      *coef_c++ = 0;
      *coef_d++ = 0;
    }

    // For bank #3, pad the coefficients with 4 zeros to avoid a 1-sample delay.
    coef_d[0] = coef_d[P+1] = coef_d[P+2] = coef_d[P+3] = 0;

    // Align on coefficients quadruples.
    switch ( M - 4*(P-1) )
    {
    case 1: 
      coef_a[0] = 0;
      coef_b[0] = 0;
      coef_c[0] = 0;
      coef_d[1] = h[4*(P-1)+0];
      break;
    case 2:
      coef_a[0] = 0;
      coef_b[0] = 0;
      coef_c[0] = h[4*(P-1)+1];
      coef_d[1] = h[4*(P-1)+0];
      break;
    case 3:
      coef_a[0] = 0;
      coef_b[0] = h[4*(P-1)+2];
      coef_c[0] = h[4*(P-1)+1];
      coef_d[1] = h[4*(P-1)+0];
      break;
    case 4:
      coef_a[0] = h[4*(P-1)+3];
      coef_b[0] = h[4*(P-1)+2];
      coef_c[0] = h[4*(P-1)+1];
      coef_d[1] = h[4*(P-1)+0];
    }

    // Copy coefficient quadruples in reverted order.
    for ( m=1; m<P; m++ )
    {
      coef_a[m  ] = h[4*(P-1-m)+3];
      coef_b[m  ] = h[4*(P-1-m)+2];
      coef_c[m  ] = h[4*(P-1-m)+1];
      coef_d[m+1] = h[4*(P-1-m)+0];
    }
  }
  else
  {
    coef_a = coefBuf;

    // Begin by a few zeros to make the number of coefficients a multiple of 4.
    for ( m=M; m<bM; m++ )
    {
      *coef_a++ = 0;
    }

    // Pad the coefficients with 4 zeros (one at the beginning, three following 
    // the end) to avoid a 1-sample delay of filter response.
    coef_a[0] = coef_a[M+1] = coef_a[M+2] = coef_a[M+3] = 0;

    // Copy coefficients in reverted order.
    for ( m=1; m<=M; m++ )
    {
      coef_a[m] = h[M-m];
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

  firdec->magic     = MAGIC;
  firdec->D         = D;
  firdec->M         = bM;
  firdec->procFxn   = procFxn;
  firdec->coef      = coefBuf;
  firdec->delayLine = delayLine;
  firdec->delayLen  = delayLen;
  firdec->delayPos  = delayLine;

  return (firdec);

} /* firdec32x16_init() */

/* Put a chunk of input signal into the delay line and compute the filter
 * response. */
void firdec32x16_process( firdec32x16_handle_t _firdec, 
                          int32_t * restrict   y, 
                    const int32_t *            x, int N )
{
  firdec32x16_ptr_t firdec = (firdec32x16_ptr_t)_firdec;

  ASSERT( firdec && firdec->magic == MAGIC && y && x );

  ASSERT( IS_ALIGN( x ) );

  //
  // Call filter's data processing function. It will store the block of input 
  // samples to the delay line, and compute the filter response. Returns the
  // updated next position pointer into the delay line buffer.
  //

  ASSERT( firdec->procFxn );
  if(N>0)
  firdec->delayPos = (*firdec->procFxn)( y,
                                         firdec->delayLine,
                                         firdec->delayLine + firdec->delayLen,
                                         firdec->delayPos,
                                         x,
                                         firdec->coef,
                                         firdec->D,
                                         N,
                                         firdec->M );

} /* firdec32x16_process() */

/* Data processing function for a factor 2 decimating FIR filter. */
static int32_t * firdec2_proc( int32_t * restrict y,
                               int32_t *          cbegin,
                               int32_t *          cend,
                               int32_t * restrict cpos,
                         const int32_t * restrict x,
                         const int16_t * restrict h,
                         int D, int N, int M )
{
        ae_int32x2 * restrict D_wr;
  const ae_int32x2 *          D_rd;
  const ae_f16x4   *          C;
  const ae_int32x2 *          X;
        ae_f32     * restrict Y;

  int32_t * cbuf1;
  int       cbufSize;

  ae_f16x4   c;
  ae_f32x2   d0, d1, d2, d3, d4, d5;
  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_int32x2 t0, t1, t2, t3, t4, t5, t6, t7;
  ae_int32x2 r00, r01, r02, r03;
  ae_int32x2 r10, r11, r12, r13;

  int n, m;

  ASSERT( y && cbegin && cend && cpos && x && h && (D == 2) && (N > 0) && (M > 0) );

  ASSERT( !( N & 7 ) && !( M & 3 ) );

  ASSERT( IS_ALIGN( cbegin ) &&
          IS_ALIGN( cend   ) &&
          IS_ALIGN( cpos   ) &&
          IS_ALIGN( x      ) &&
          IS_ALIGN( h      ) );

  //
  // Setup pointers and 2-bank circular delay line buffers.
  //

  D_wr = (      ae_int32x2 *)cpos;
  X    = (const ae_int32x2 *)x;
  Y    = (      ae_f32     *)y;

  // Bank delay line size, in bytes. M holds the number of taps per bank.
  cbufSize = 4*( M + 8 );

  // cbegin points at the delay line of bank 0. Make pointer for bank #1.
  cbuf1 = (int32_t*)( (uintptr_t)cbegin + cbufSize );

  //
  // Break the input signal into 16-samples blocks. For each block, decimate the
  // samples with both possible offsets, and obtain 2 8-sample chunks. Store
  // each chunk to the delay line of the corresponding filter bank (#0, #1).
  // After that compute 2-banks filter response.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 16 input samples.
    // Q31
    AE_L32X2_IP( t0,  X, +8 );
    AE_L32X2_IP( t1,  X, +8 );
    AE_L32X2_IP( t2,  X, +8 );
    AE_L32X2_IP( t3,  X, +8 );
    AE_L32X2_IP( t4,  X, +8 );
    AE_L32X2_IP( t5,  X, +8 );
    AE_L32X2_IP( t6,  X, +8 );
    AE_L32X2_IP( t7,  X, +8 );

    //
    // To preserve bit-exactness with the reference code, we have to use the 
    // decimation phase 1. Thus the first sample in each triple goes to
    // bank #1. To keep both filter banks in-phase, the response of bank #0
    // is implicitly delayed by one sample.
    //

    // Input samples for bank #1.
    r10 = AE_SEL32_HH( t0, t1 ); //  0,  2
    r11 = AE_SEL32_HH( t2, t3 ); //  4,  6
    r12 = AE_SEL32_HH( t4, t5 ); //  8, 10
    r13 = AE_SEL32_HH( t6, t7 ); // 12, 14

    // Input samples for bank #0
    r00 = AE_SEL32_LL( t0, t1 ); //  1,  3
    r01 = AE_SEL32_LL( t2, t3 ); //  5,  7
    r02 = AE_SEL32_LL( t4, t5 ); //  9, 11
    r03 = AE_SEL32_LL( t6, t7 ); // 13, 15

    // Setup circular addressing for bank #0. D_wr always looks in bank #0.
    WUR_AE_CBEGIN0( (uintptr_t)cbegin );
    WUR_AE_CEND0  ( (uintptr_t)cbuf1  );

    // Q31
    AE_S32X2_X( r10, D_wr, cbufSize ); // #1 
    AE_S32X2_XC( r00, D_wr, +8       ); // #0

    // Q31
    AE_S32X2_X( r11, D_wr, cbufSize ); // #1
    AE_S32X2_XC( r01, D_wr, +8       ); // #0

    // Q31
    AE_S32X2_X( r12, D_wr, cbufSize ); // #1
    AE_S32X2_XC( r02, D_wr, +8       ); // #0

    // Q31
    AE_S32X2_X( r13, D_wr, cbufSize ); // #1
    AE_S32X2_XC( r03, D_wr, +8       ); // #0

    // Reset the coefficients pointer. Now it looks at the tap corresponding
    // to the oldest sample in the delay line of bank #0. The structure of the
    // coefficients block (2*M+4 entries):
    //   <M coefs of #0> 0 <M coefs of #1> 0 0 0
    C = (const ae_f16x4*)h;

    //--------------------------------------------------------------------------
    // Bank #0
    //

    // Circular buffer pointer looks at the oldest sample in bank #0:
    // M+8 samples back from the newest one.
    D_rd = (const ae_int32x2 *)D_wr;

    // Load 8 oldest samples from the delay line.
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
    // Inner loop prologue: process the first 4 taps for each of 8 accumulators.
    //

    // Load next 4 samples. Altogether we have 12 samples residing in 6 AE
    // registers.
    // Q31
    AE_L32X2_XC( t4, D_rd, +8 );
    AE_L32X2_XC( t5, D_rd, +8 );

    d4 = ( t4 );
    d5 = ( t5 );

    // Load 4 tap coefficients.
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

    // First 4 taps are done. Move 2 input samples out of the registers file.
    d0 = d2; d1 = d3; d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 2 taps for 8 accumulators on each trip. 
    // Filter's response is implicitly delayed by 1 sample.
    //

    for ( m=0; m<(M>>2)-1; m++ )
    {
      // Load next 4 samples. 
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

      // 4 taps are done. Move 2 input samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    //--------------------------------------------------------------------------
    // Bank #1
    //

    // Setup circular addressing for bank #1.
    WUR_AE_CBEGIN0( (uintptr_t)cbuf1 );
    WUR_AE_CEND0  ( (uintptr_t)cend  );

    // Circular buffer pointer looks at the oldest sample in bank #1:
    // M+8 samples back from the newest one.
    D_rd = (const ae_int32x2 *)( (uintptr_t)D_wr + cbufSize );

    // Load 8 oldest samples from the delay line.
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
    // Inner loop: process 4 taps for 8 accumulators on each trip, add the 
    // result to the output samples of bank #0. There is no delay  introduced
    // to the response of filter bank #1. To accomplish this zero-delay we 
    // actually perform M+4 MACs for each accumulator, with 3 MACs falling
    // on zero coefficients inserted in the impulse response of bank #0 during
    // initialization.
    //

    for ( m=0; m<(M>>2)+1; m++ )
    {
      // Load next 4 samples. 
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

      // 4 taps are done. Move 2 input samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    // Convert and store filter outputs.
    // Q31 <- Q16.47 - 16 w/ rounding and saturation.
    AE_S32RA64S_IP( q0, Y, +4 );
    AE_S32RA64S_IP( q1, Y, +4 );
    AE_S32RA64S_IP( q2, Y, +4 );
    AE_S32RA64S_IP( q3, Y, +4 );
    AE_S32RA64S_IP( q4, Y, +4 );
    AE_S32RA64S_IP( q5, Y, +4 );
    AE_S32RA64S_IP( q6, Y, +4 );
    AE_S32RA64S_IP( q7, Y, +4 );
  }

  return ( (int32_t*)D_wr );

} /* firdec2_proc() */

/* Data processing function for a factor 3 decimating FIR filter. */
static int32_t * firdec3_proc( int32_t * restrict y,
                               int32_t *          cbegin,
                               int32_t *          cend,
                               int32_t * restrict cpos,
                         const int32_t * restrict x,
                         const int16_t * restrict h,
                         int D, int N, int M )
{
        ae_int32x2 * restrict D_wr;
  const ae_int32x2 *          D_rd;
  const ae_f16x4   *          C;
  const ae_int32x2 *          X;
        ae_f32     * restrict Y;

  int32_t * cbuf1;
  int32_t * cbuf2;
  int       cbufSize;

  ae_f16x4   c;
  ae_f32x2   d0, d1, d2, d3, d4, d5;
  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_int32x2 t0, t1, t2, t3, t4, t5;
  ae_int32x2 t6, t7, t8, t9, t10, t11;
  ae_int32x2 r00, r01, r02, r03;
  ae_int32x2 r10, r11, r12, r13;
  ae_int32x2 r20, r21, r22, r23;

  int n, m;

  ASSERT( y && cbegin && cend && cpos && x && h
          && 
          (D == 3) && (N > 0) && (M > 0) );

  ASSERT( !( N & 7 ) && !( M & 3 ) );

  ASSERT( IS_ALIGN( cbegin ) &&
          IS_ALIGN( cend   ) &&
          IS_ALIGN( cpos   ) &&
          IS_ALIGN( x      ) &&
          IS_ALIGN( h      ) );

  //
  // Setup pointers and 3-bank circular delay line buffers.
  //

  D_wr = (      ae_int32x2 *)cpos;
  X    = (const ae_int32x2 *)x;
  Y    = (      ae_f32     *)y;

  // Bank delay line size, in bytes. M holds the number of taps per bank.
  cbufSize = 4*( M + 8 );

  // cbegin points at the delay line of bank 0. Make pointers for banks #1, #2.
  cbuf1 = (int32_t*)( (uintptr_t)cbegin + 1*cbufSize );
  cbuf2 = (int32_t*)( (uintptr_t)cbegin + 2*cbufSize );

  //
  // Break the input signal into 24-samples blocks. For each block, decimate the
  // samples with all three possible offsets, and obtain 3 8-sample chunks.
  // Store each chunk to the delay line of the corresponding filter bank (#0,
  // #1, #2). After that compute 3-banks filter response.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 24 input samples.
    // Q31
    AE_L32X2_IP( t0,  X, +8 );
    AE_L32X2_IP( t1,  X, +8 );
    AE_L32X2_IP( t2,  X, +8 );
    AE_L32X2_IP( t3,  X, +8 );
    AE_L32X2_IP( t4,  X, +8 );
    AE_L32X2_IP( t5,  X, +8 );
    AE_L32X2_IP( t6,  X, +8 );
    AE_L32X2_IP( t7,  X, +8 );
    AE_L32X2_IP( t8,  X, +8 );
    AE_L32X2_IP( t9,  X, +8 );
    AE_L32X2_IP( t10, X, +8 );
    AE_L32X2_IP( t11, X, +8 );

    //
    // To preserve bit-exactness with the reference code, we have to use the 
    // decimation phase 2. Thus the first sample in each triple goes to
    // bank #2. To keep all three filter banks in-phase, the response of banks
    // #0, #1 is implicitly delayed by one sample.
    //

    // Input samples for bank #2. 
    r20 = AE_SEL32_HL(  t0, t1  ); //  0,  3
    r21 = AE_SEL32_HL(  t3, t4  ); //  6,  9
    r22 = AE_SEL32_HL(  t6, t7  ); // 12, 15
    r23 = AE_SEL32_HL(  t9, t10 ); // 18, 21

    // Input samples for bank #0.
    r00 = AE_SEL32_LH(  t0,  t2 ); //  1.  4
    r01 = AE_SEL32_LH(  t3,  t5 ); //  7. 10
    r02 = AE_SEL32_LH(  t6,  t8 ); // 13, 16
    r03 = AE_SEL32_LH(  t9, t11 ); // 19, 22

    // Input samples for bank #1.
    r10 = AE_SEL32_HL(  t1,  t2 ); //  2,  5
    r11 = AE_SEL32_HL(  t4,  t5 ); //  8, 11
    r12 = AE_SEL32_HL(  t7,  t8 ); // 14, 17
    r13 = AE_SEL32_HL( t10, t11 ); // 20, 23

    // Setup circular addressing for bank #0. D_wr always looks in bank #0.
    WUR_AE_CBEGIN0( (uintptr_t)cbegin );
    WUR_AE_CEND0  ( (uintptr_t)cbuf1  );

    // Q31
    AE_S32X2_X( r20, D_wr, 2*cbufSize ); // #2
    AE_S32X2_X( r10, D_wr, 1*cbufSize ); // #1 
    AE_S32X2_XC( r00, D_wr, +8         ); // #0

    // Q31
    AE_S32X2_X( r21, D_wr, 2*cbufSize ); // #2
    AE_S32X2_X( r11, D_wr, 1*cbufSize ); // #1
    AE_S32X2_XC( r01, D_wr, +8         ); // #0

    // Q31
    AE_S32X2_X( r22, D_wr, 2*cbufSize ); // #2
    AE_S32X2_X( r12, D_wr, 1*cbufSize ); // #1
    AE_S32X2_XC( r02, D_wr, +8         ); // #0

    // Q31
    AE_S32X2_X( r23, D_wr, 2*cbufSize ); // #2
    AE_S32X2_X( r13, D_wr, 1*cbufSize ); // #1
    AE_S32X2_XC( r03, D_wr, +8         ); // #0

    // Reset the coefficients pointer. Now it looks at the tap corresponding
    // to the oldest sample in the delay line of bank #0. The structure of the
    // coefficients block (3*M+4 entries):
    //   <M coefs of #0> <M coefs of #1> 0 <M coefs of #2> 0 0 0
    C = (const ae_f16x4*)h;

    //--------------------------------------------------------------------------
    // Bank #0
    //

    // Circular buffer pointer looks at the oldest sample in bank #0:
    // M+8 samples back from the newest one.
    D_rd = (const ae_int32x2 *)D_wr;

    // Load 8 oldest samples from the delay line.
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
    // Inner loop prologue: process the first 4 taps for each of 8 accumulators.
    //

    // Load next 4 samples. Altogether we have 12 samples residing in 6 AE
    // registers.
    // Q31
    AE_L32X2_XC( t4, D_rd, +8 );
    AE_L32X2_XC( t5, D_rd, +8 );

    d4 = ( t4 );
    d5 = ( t5 );

    // Load 4 tap coefficients.
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

    // First 4 taps are done. Move 4 input samples out of the registers file.
    d0 = d2; d1 = d3; d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators on each trip. 
    // Filter's response is implicitly delayed by 1 sample.
    //

    for ( m=0; m<(M>>2)-1; m++ )
    {
      // Load next 4 samples. 
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

      // 4 taps are done. Move 4 input samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    //--------------------------------------------------------------------------
    // Bank #1
    //

    // Setup circular addressing for bank #1.
    WUR_AE_CBEGIN0( (uintptr_t)cbuf1 );
    WUR_AE_CEND0  ( (uintptr_t)cbuf2 );

    // Circular buffer pointer looks at the oldest sample in bank #1:
    // M+8 samples back from the newest one.
    D_rd = (const ae_int32x2 *)( (uintptr_t)D_wr + 1*cbufSize );

    // Load 8 oldest samples from the delay line.
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
    // Inner loop: process 4 taps for 8 accumulators on each trip, add the 
    // result to the output samples of bank #0. Filter's response is implicitly
    // delayed by 1 sample. 
    //

    for ( m=0; m<(M>>2); m++ )
    {
      // Load next 4 samples. 
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

      // 4 taps are done. Move 4 input samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    //--------------------------------------------------------------------------
    // Bank #2
    //

    // Setup circular addressing for bank #2.
    WUR_AE_CBEGIN0( (uintptr_t)cbuf2 );
    WUR_AE_CEND0  ( (uintptr_t)cend  );

    // Circular buffer pointer looks at the oldest sample in bank #2:
    // M+8 samples back from the newest one.
    D_rd = (const ae_int32x2 *)( (uintptr_t)D_wr + 2*cbufSize );

    // Load 8 oldest samples from the delay line.
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
    // Inner loop: process 4 taps for 8 accumulators on each trip, add the 
    // result to the output samples of banks #0, #1. There is no delay 
    // introduced to the response of filter bank #2. To accomplish this zero-
    // delay we actually perform M+4 MACs for each accumulator, with 4 MACs
    // falling on zero coefficients inserted in the impulse response of bank
    // #2 during initialization.
    //

    for ( m=0; m<(M>>2)+1; m++ )
    {
      // Load next 4 samples. 
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

      // 4 taps are done. Move 4 input samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    // Convert and store 8 filter outputs.
    // Q31 <- Q16.47 - 16 w/ rounding and saturation
    AE_S32RA64S_IP( q0, Y, +4 );
    AE_S32RA64S_IP( q1, Y, +4 );
    AE_S32RA64S_IP( q2, Y, +4 );
    AE_S32RA64S_IP( q3, Y, +4 );
    AE_S32RA64S_IP( q4, Y, +4 );
    AE_S32RA64S_IP( q5, Y, +4 );
    AE_S32RA64S_IP( q6, Y, +4 );
    AE_S32RA64S_IP( q7, Y, +4 );
  }

  return ( (int32_t*)D_wr );

} /* firdec3_proc() */

/* Data processing function for a factor 4 decimating FIR filter. */
static int32_t * firdec4_proc( int32_t * restrict y,
                               int32_t *          cbegin,
                               int32_t *          cend,
                               int32_t * restrict cpos,
                         const int32_t * restrict x,
                         const int16_t * restrict h,
                         int D, int N, int M )
{
        ae_int32x2 * restrict D_wr;
  const ae_int32x2 *          D_rd;
  const ae_f16x4   *          C;
  const ae_int32x2 *          X;
        ae_f32     * restrict Y;

  int32_t * cbuf1;
  int32_t * cbuf2;
  int32_t * cbuf3;
  int       cbufSize;

  ae_f16x4   c;
  ae_f32x2   d0, d1, d2, d3, d4, d5;
  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_int32x2 t0, t1, t2, t3, t4, t5, t6, t7;
  ae_int32x2 t8, t9, t10, t11, t12, t13, t14, t15;
  ae_int32x2 r00, r01, r02, r03;
  ae_int32x2 r10, r11, r12, r13;
  ae_int32x2 r20, r21, r22, r23;
  ae_int32x2 r30, r31, r32, r33;

  int n, m;

  ASSERT( y && cbegin && cend && cpos && x && h
          && 
          D == 4 && N > 0 && M > 0 );

  ASSERT( !( N & 7 ) && !( M & 3 ) );

  ASSERT( IS_ALIGN( cbegin ) &&
          IS_ALIGN( cend   ) &&
          IS_ALIGN( cpos   ) &&
          IS_ALIGN( x      ) &&
          IS_ALIGN( h      ) );

  //
  // Setup pointers and 4-bank circular delay line buffers.
  //

  D_wr = (      ae_int32x2 *)cpos;
  X    = (const ae_int32x2 *)x;
  Y    = (      ae_f32     *)y;

  // Bank delay line size, in bytes. M holds the number of taps per bank.
  cbufSize = 4*( M + 8 );

  // cbegin points at the delay line of bank 0. Make pointers for banks #1,#2,#3
  cbuf1 = (int32_t*)( (uintptr_t)cbegin + 1*cbufSize );
  cbuf2 = (int32_t*)( (uintptr_t)cbegin + 2*cbufSize );
  cbuf3 = (int32_t*)( (uintptr_t)cbegin + 3*cbufSize );

  //
  // Break the input signal into 32-samples blocks. For each block, decimate the
  // samples with all four possible offsets, and obtain 4 8-sample chunks.
  // Store each chunk to the delay line of the corresponding filter bank (#0,
  // #1, #2, #3). After that compute 4-banks filter response.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 24 input samples.
    // Q31
    AE_L32X2_IP( t0,  X, +8 );
    AE_L32X2_IP( t1,  X, +8 );
    AE_L32X2_IP( t2,  X, +8 );
    AE_L32X2_IP( t3,  X, +8 );
    AE_L32X2_IP( t4,  X, +8 );
    AE_L32X2_IP( t5,  X, +8 );
    AE_L32X2_IP( t6,  X, +8 );
    AE_L32X2_IP( t7,  X, +8 );
    AE_L32X2_IP( t8,  X, +8 );
    AE_L32X2_IP( t9,  X, +8 );
    AE_L32X2_IP( t10, X, +8 );
    AE_L32X2_IP( t11, X, +8 );
    AE_L32X2_IP( t12, X, +8 );
    AE_L32X2_IP( t13, X, +8 );
    AE_L32X2_IP( t14, X, +8 );
    AE_L32X2_IP( t15, X, +8 );

    //
    // To preserve bit-exactness with the reference code, we have to use the 
    // decimation phase 3. Thus the first sample in each triple goes to
    // bank #3. To keep all four filter banks in-phase, the response of banks
    // #0, #1, #2 is implicitly delayed by one sample.
    //

    // Input samples for bank #3. 
    r30 = AE_SEL32_HH(  t0,  t2 ); //  0,  4
    r31 = AE_SEL32_HH(  t4,  t6 ); //  8, 12
    r32 = AE_SEL32_HH(  t8, t10 ); // 16, 20
    r33 = AE_SEL32_HH( t12, t14 ); // 24, 28

    // Input samples for bank #0. 
    r00 = AE_SEL32_LL(  t0,  t2 ); //  1,  5
    r01 = AE_SEL32_LL(  t4,  t6 ); //  9, 13
    r02 = AE_SEL32_LL(  t8, t10 ); // 17, 21
    r03 = AE_SEL32_LL( t12, t14 ); // 25, 29

    // Input samples for bank #1. 
    r10 = AE_SEL32_HH(  t1,  t3 ); //  2,  6
    r11 = AE_SEL32_HH(  t5,  t7 ); // 10, 14
    r12 = AE_SEL32_HH(  t9, t11 ); // 18, 22
    r13 = AE_SEL32_HH( t13, t15 ); // 26, 30

    // Input samples for bank #2. 
    r20 = AE_SEL32_LL(  t1,  t3 ); //  3,  7
    r21 = AE_SEL32_LL(  t5,  t7 ); // 11, 15
    r22 = AE_SEL32_LL(  t9, t11 ); // 19, 23
    r23 = AE_SEL32_LL( t13, t15 ); // 27, 31

    // Setup circular addressing for bank #0. D_wr always looks in bank #0.
    WUR_AE_CBEGIN0( (uintptr_t)cbegin );
    WUR_AE_CEND0  ( (uintptr_t)cbuf1  );

    // Q31
    AE_S32X2_X( r30, D_wr, 3*cbufSize ); // #3
    AE_S32X2_X( r20, D_wr, 2*cbufSize ); // #2
    AE_S32X2_X( r10, D_wr, 1*cbufSize ); // #1
    AE_S32X2_XC( r00, D_wr, +8         ); // #0

    // Q31
    AE_S32X2_X( r31, D_wr, 3*cbufSize ); // #3
    AE_S32X2_X( r21, D_wr, 2*cbufSize ); // #2
    AE_S32X2_X( r11, D_wr, 1*cbufSize ); // #1
    AE_S32X2_XC( r01, D_wr, +8         ); // #0

    // Q31
    AE_S32X2_X( r32, D_wr, 3*cbufSize ); // #3
    AE_S32X2_X( r22, D_wr, 2*cbufSize ); // #2
    AE_S32X2_X( r12, D_wr, 1*cbufSize ); // #1
    AE_S32X2_XC( r02, D_wr, +8         ); // #0

    // Q31
    AE_S32X2_X( r33, D_wr, 3*cbufSize ); // #3
    AE_S32X2_X( r23, D_wr, 2*cbufSize ); // #2
    AE_S32X2_X( r13, D_wr, 1*cbufSize ); // #1
    AE_S32X2_XC( r03, D_wr, +8         ); // #0

    // Reset the coefficients pointer. Now it looks at the tap corresponding
    // to the oldest sample in the delay line of bank #0. The structure of the
    // coefficients block (3*M+4 entries):
    //   <M coefs of #0> <M coefs of #1> 0 <M coefs of #2> 0 0 0
    C = (const ae_f16x4*)h;

    //--------------------------------------------------------------------------
    // Bank #0
    //

    // Circular buffer pointer looks at the oldest sample in bank #0:
    // M+8 samples back from the newest one.
    D_rd = (const ae_int32x2 *)D_wr;

    // Load 8 oldest samples from the delay line.
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
    // Inner loop prologue: process the first 4 taps for each of 8 accumulators.
    //

    // Load next 4 samples. Altogether we have 12 samples residing in 6 AE
    // registers.
    // Q31
    AE_L32X2_XC( t4, D_rd, +8 );
    AE_L32X2_XC( t5, D_rd, +8 );

    d4 = ( t4 );
    d5 = ( t5 );

    // Load 4 tap coefficients.
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

    // First 4 taps are done. Move 4 input samples out of the registers file.
    d0 = d2; d1 = d3; 
    d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators on each trip. 
    // Filter's response is implicitly delayed by 1 sample.
    //

    for ( m=0; m<(M>>2)-1; m++ )
    {
      // Load next 4 samples. 
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

      // 4 taps are done. Move 4 input samples out of the registers file.
      d0 = d2; d1 = d3;
      d2 = d4; d3 = d5;
    }

    //--------------------------------------------------------------------------
    // Bank #1
    //

    // Setup circular addressing for bank #1.
    WUR_AE_CBEGIN0( (uintptr_t)cbuf1 );
    WUR_AE_CEND0  ( (uintptr_t)cbuf2 );

    // Circular buffer pointer looks at the oldest sample in bank #1:
    // M+8 samples back from the newest one.
    D_rd = (const ae_int32x2 *)( (uintptr_t)D_wr + 1*cbufSize );

    // Load 8 oldest samples from the delay line.
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
    // Inner loop: process 4 taps for 8 accumulators on each trip, add the 
    // result to the output samples of bank #0. Filter's response is implicitly
    // delayed by 1 sample. 
    //

    for ( m=0; m<(M>>2); m++ )
    {
      // Load next 4 samples. 
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

      // 4 taps are done. Move 4 input samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    //--------------------------------------------------------------------------
    // Bank #2
    //

    // Setup circular addressing for bank #2.
    WUR_AE_CBEGIN0( (uintptr_t)cbuf2 );
    WUR_AE_CEND0  ( (uintptr_t)cbuf3 );

    // Circular buffer pointer looks at the oldest sample in bank #2:
    // M+8 samples back from the newest one.
    D_rd = (const ae_int32x2 *)( (uintptr_t)D_wr + 2*cbufSize );

    // Load 8 oldest samples from the delay line.
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
    // Inner loop: process 4 taps for 8 accumulators on each trip, add the 
    // result to the output samples of banks #0, #1. Filter's response is
    // implicitly delayed by 1 sample. 
    //

    for ( m=0; m<(M>>2); m++ )
    {
      // Load next 4 samples. 
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

      // 4 taps are done. Move 4 input samples out of the registers file.
      d0 = d2; d1 = d3; d2 = d4; d3 = d5;
    }

    //--------------------------------------------------------------------------
    // Bank #3
    //

    // Setup circular addressing for bank #3.
    WUR_AE_CBEGIN0( (uintptr_t)cbuf3 );
    WUR_AE_CEND0  ( (uintptr_t)cend  );

    // Circular buffer pointer looks at the oldest sample in bank #3:
    // M+8 samples back from the newest one.
    D_rd = (const ae_int32x2 *)( (uintptr_t)D_wr + 3*cbufSize );

    // Load 8 oldest samples from the delay line.
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
    // Inner loop: process 4 taps for 8 accumulators on each trip, add the 
    // result to the output samples of banks #0, #1, #2. There is no delay 
    // introduced to the response of filter bank #3. To accomplish this zero-
    // delay we actually perform M+4 MACs for each accumulator, with 4 MACs
    // falling on zero coefficients inserted in the impulse response of bank
    // #3 during initialization.
    //

    for ( m=0; m<(M>>2)+1; m++ )
    {
      // Load next 4 samples. 
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

      // 4 taps are done. Move 4 input samples out of the registers file.
      d0 = d2; d1 = d3; 
      d2 = d4; d3 = d5;
    }

    // Convert and store 8 filter outputs.
    // Q31 <- Q16.47 - 16 w/ rounding and saturation
    AE_S32RA64S_IP( q0, Y, +4 );
    AE_S32RA64S_IP( q1, Y, +4 );
    AE_S32RA64S_IP( q2, Y, +4 );
    AE_S32RA64S_IP( q3, Y, +4 );
    AE_S32RA64S_IP( q4, Y, +4 );
    AE_S32RA64S_IP( q5, Y, +4 );
    AE_S32RA64S_IP( q6, Y, +4 );
    AE_S32RA64S_IP( q7, Y, +4 );
  }

  return ( (int32_t*)D_wr );

} /* firdec4_proc() */

/* Generic data processing function for a decimating FIR filter. */
static int32_t * firdecX_proc( int32_t * restrict y,
                               int32_t *          cbegin,
                               int32_t *          cend,
                               int32_t * restrict cpos,
                         const int32_t * restrict x,
                         const int16_t * restrict h,
                         int D, int N, int M )
{
  const ae_int32   *          D_tmp;
  const ae_int32x2 *          D_rd0;
  const ae_int32x2 *          D_rd1;
  const ae_int32x2 *          D_rd2;
  const ae_int32x2 *          D_rd3;
  const ae_int32x2 *          D_rd4;
  const ae_int32x2 *          D_rd5;
  const ae_int32x2 *          D_rd6;
  const ae_int32x2 *          D_rd7;
        ae_int32x2 * restrict D_wr;
  const ae_int32x2 *          X;
        ae_f32     * restrict Y;
  const ae_f16x4   *          C;

  ae_valign D_va1, D_va3, D_va5, D_va7;

  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_f32x2   d0, d1, d2, d3, d4, d5, d6, d7;
  ae_int32x2 t0, t1, t2, t3, t4, t5, t6, t7;
  ae_f16x4   c0;

  int n, m;

  ASSERT( y && cbegin && cend && cpos && x && h && (D > 1) && (N > 0) && (M > 0) );

  ASSERT( !( N & 7 ) && !( M & 1 ) );

  ASSERT( IS_ALIGN( cbegin ) &&
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
  // Break the input signal into 8*D-samples blocks. For each block, store
  // 8*D samples to the delay line buffer, and compute 8 samples of decimated
  // response signal.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    for ( m=0; m<(8>>1)*D; m++ )
    {
      AE_L32X2_IP( t0,    X, +8 );
      AE_S32X2_XC ( t0, D_wr, +8 );
    }

    //
    // Setup 8-way delay line reading pointers, one per an accumulator.
    //

    D_tmp = (const ae_int32   *)D_wr;
    D_rd0 = (const ae_int32x2 *)D_wr;

    AE_L32_XC( t0, D_tmp, 4*D );
    D_rd1 = (const ae_int32x2*)D_tmp;

    AE_L32_XC( t0, D_tmp, 4*D );
    D_rd2 = (const ae_int32x2*)D_tmp;

    AE_L32_XC( t0, D_tmp, 4*D );
    D_rd3 = (const ae_int32x2*)D_tmp;

    AE_L32_XC( t0, D_tmp, 4*D );
    D_rd4 = (const ae_int32x2*)D_tmp;

    AE_L32_XC( t0, D_tmp, 4*D );
    D_rd5 = (const ae_int32x2*)D_tmp;

    AE_L32_XC( t0, D_tmp, 4*D );
    D_rd6 = (const ae_int32x2*)D_tmp;

    AE_L32_XC( t0, D_tmp, 4*D );
    D_rd7 = (const ae_int32x2*)D_tmp;

    AE_LA32X2POS_PC( D_va1, D_rd1 );
    AE_LA32X2POS_PC( D_va3, D_rd3 );
    AE_LA32X2POS_PC( D_va5, D_rd5 );
    AE_LA32X2POS_PC( D_va7, D_rd7 );

    //
    // Inner loop prologue: process the first 4 taps for each of 8 accumulators.
    //

    // Reset the coefficients pointer. Now it looks at the tap corresponding
    // to the oldest sample.
    C = (const ae_f16x4*)h;
    // Load 4 tap coefficients.
    // Q15
    ae_f16x4_loadip( c0, C, +8 );

    // Load 2 samples for each even-numbered accumulator.
    // Q31
    AE_L32X2_XC( t0, D_rd0, +8 );
    AE_L32X2_XC( t2, D_rd2, +8 );
    AE_L32X2_XC( t4, D_rd4, +8 );
    AE_L32X2_XC( t6, D_rd6, +8 );

    d0 = t0;
    d2 = t2;
    d4 = t4;
    d6 = t6;

    // Load 2 samples for each odd-numbered accumulator.
    // Q31
    AE_LA32X2_IC( t1, D_va1, D_rd1 );
    AE_LA32X2_IC( t3, D_va3, D_rd3 );
    AE_LA32X2_IC( t5, D_va5, D_rd5 );
    AE_LA32X2_IC( t7, D_va7, D_rd7 );

    d1 = t1;
    d3 = t3;
    d5 = t5;
    d7 = t7;

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    q0 = AE_MULZAAFD32X16_H3_L2( d0, c0 );
    q1 = AE_MULZAAFD32X16_H3_L2( d1, c0 );
    q2 = AE_MULZAAFD32X16_H3_L2( d2, c0 );
    q3 = AE_MULZAAFD32X16_H3_L2( d3, c0 );
    q4 = AE_MULZAAFD32X16_H3_L2( d4, c0 );
    q5 = AE_MULZAAFD32X16_H3_L2( d5, c0 );
    q6 = AE_MULZAAFD32X16_H3_L2( d6, c0 );
    q7 = AE_MULZAAFD32X16_H3_L2( d7, c0 );

    // Load 2 samples for each even-numbered accumulator.
    // Q31
    AE_L32X2_XC( t0, D_rd0, +8 );
    AE_L32X2_XC( t2, D_rd2, +8 );
    AE_L32X2_XC( t4, D_rd4, +8 );
    AE_L32X2_XC( t6, D_rd6, +8 );

    d0 = t0;
    d2 = t2;
    d4 = t4;
    d6 = t6;

    // Load 2 samples for each odd-numbered accumulator.
    // Q31
    AE_LA32X2_IC(t1, D_va1, D_rd1);
    AE_LA32X2_IC(t3, D_va3, D_rd3);
    AE_LA32X2_IC(t5, D_va5, D_rd5);
    AE_LA32X2_IC(t7, D_va7, D_rd7);

    d1 = t1;
    d3 = t3;
    d5 = t5;
    d7 = t7;

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAAFD32X16_H1_L0( q0, d0, c0 );
    AE_MULAAFD32X16_H1_L0( q1, d1, c0 );
    AE_MULAAFD32X16_H1_L0( q2, d2, c0 );
    AE_MULAAFD32X16_H1_L0( q3, d3, c0 );
    AE_MULAAFD32X16_H1_L0( q4, d4, c0 );
    AE_MULAAFD32X16_H1_L0( q5, d5, c0 );
    AE_MULAAFD32X16_H1_L0( q6, d6, c0 );
    AE_MULAAFD32X16_H1_L0( q7, d7, c0 );

    //
    // Inner loop: process 4 taps for 8 accumulators on each trip. Totally we 
    // perform M+4 MACs for each accumulator, 4 of which fall on zero taps
    // inserted into the impulse response during initialization.
    //

    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip( c0, C, +8 );

      // Load 2 samples for each even-numbered accumulator.
      // Q31
      AE_L32X2_XC(t0, D_rd0, +8);
      AE_L32X2_XC(t2, D_rd2, +8);
      AE_L32X2_XC(t4, D_rd4, +8);
      AE_L32X2_XC(t6, D_rd6, +8);

      d0 = t0;
      d2 = t2;
      d4 = t4;
      d6 = t6;

      // Load 2 samples for each odd-numbered accumulator.
      // Q31
      AE_LA32X2_IC( t1, D_va1, D_rd1 );
      AE_LA32X2_IC( t3, D_va3, D_rd3 );
      AE_LA32X2_IC( t5, D_va5, D_rd5 );
      AE_LA32X2_IC( t7, D_va7, D_rd7 );

      d1 = t1;
      d3 = t3;
      d5 = t5;
      d7 = t7;

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAAFD32X16_H3_L2( q0, d0, c0 );
      AE_MULAAFD32X16_H3_L2( q1, d1, c0 );
      AE_MULAAFD32X16_H3_L2( q2, d2, c0 );
      AE_MULAAFD32X16_H3_L2( q3, d3, c0 );
      AE_MULAAFD32X16_H3_L2( q4, d4, c0 );
      AE_MULAAFD32X16_H3_L2( q5, d5, c0 );
      AE_MULAAFD32X16_H3_L2( q6, d6, c0 );
      AE_MULAAFD32X16_H3_L2( q7, d7, c0 );

      // Load 2 samples for each even-numbered accumulator.
      // Q31
      AE_L32X2_XC( t0, D_rd0, +8 );
      AE_L32X2_XC( t2, D_rd2, +8 );
      AE_L32X2_XC( t4, D_rd4, +8 );
      AE_L32X2_XC( t6, D_rd6, +8 );

      d0 = t0;
      d2 = t2;
      d4 = t4;
      d6 = t6;

      // Load 2 samples for each odd-numbered accumulator.
      // Q31
      AE_LA32X2_IC( t1, D_va1, D_rd1 );
      AE_LA32X2_IC( t3, D_va3, D_rd3 );
      AE_LA32X2_IC( t5, D_va5, D_rd5 );
      AE_LA32X2_IC( t7, D_va7, D_rd7 );

      d1 = t1;
      d3 = t3;
      d5 = t5;
      d7 = t7;

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAAFD32X16_H1_L0( q0, d0, c0 );
      AE_MULAAFD32X16_H1_L0( q1, d1, c0 );
      AE_MULAAFD32X16_H1_L0( q2, d2, c0 );
      AE_MULAAFD32X16_H1_L0( q3, d3, c0 );
      AE_MULAAFD32X16_H1_L0( q4, d4, c0 );
      AE_MULAAFD32X16_H1_L0( q5, d5, c0 );
      AE_MULAAFD32X16_H1_L0( q6, d6, c0 );
      AE_MULAAFD32X16_H1_L0( q7, d7, c0 );
    }

    // Convert and store filter outputs.
    // Q31 <- Q16.47 - 16 w/ rounding and saturation.
    AE_S32RA64S_IP( q0, Y, +4 );
    AE_S32RA64S_IP( q1, Y, +4 );
    AE_S32RA64S_IP( q2, Y, +4 );
    AE_S32RA64S_IP( q3, Y, +4 );
    AE_S32RA64S_IP( q4, Y, +4 );
    AE_S32RA64S_IP( q5, Y, +4 );
    AE_S32RA64S_IP( q6, Y, +4 );
    AE_S32RA64S_IP( q7, Y, +4 );
  }

  return (int32_t*)D_wr;

} /* firdecX_proc() */
