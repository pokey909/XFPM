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
    Interpolating block real FIR filter, 16x16-bit
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
#define MAGIC     0x589FE2D6

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
    ((size_t)(size)+(align)-1)

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

#define sz_i16    sizeof(int16_t)

#if !(defined(AE_MULAAAAQ16))
/* code for no Quad MAC option */
/* Interpolator instance structure. */
typedef int16_t * (proc_fxn_t)( int16_t * restrict y,
                                int16_t *          cbegin,
                                int16_t *          cend,
                                int16_t * restrict cpos,
                          const int16_t * restrict x,
                          const int16_t * restrict h,
                          int D, int N, int M );

static proc_fxn_t firinterpX_proc;
static proc_fxn_t firinterp2_proc;
static proc_fxn_t firinterp3_proc;
static proc_fxn_t firinterp4_proc;

/* Interpolator instance structure. */
typedef struct tag_firinterp16x16_t
{
  uint32_t        magic;     // Instance pointer validation number
  int             D;         // Interpolation factor
  int             M;         // Number of filter coefficients
  proc_fxn_t *    procFxn;   // Filter data processing function
  const int16_t * coef;      // Filter coefficients
  int16_t *       delayLine; // Delay line for complex samples
  int             delayLen;  // Delay line length, in complex samples
  int16_t *       delayPos;  // Delay line slot to be filled next

} firinterp16x16_t, *firinterp16x16_ptr_t;

/* Calculate the memory block size for an interpolator with given 
 * attributes. */
size_t firinterp16x16_alloc( int D, int M )
{
  int delayLen, coefNum;

  ASSERT( (D>1)  && (M>0) );

  ASSERT(!(M&3));

  delayLen = (M + 8);
  coefNum  = (M + 4)*D;

  return ( ALIGNED_SIZE(sizeof(firinterp16x16_t), 4)
           + // Delay line
           ALIGNED_SIZE(delayLen*sz_i16, 8)
           + // Coefficients
           ALIGNED_SIZE(coefNum*sz_i16, 8));

} /* firinterp16x16_alloc() */

/* Initialize the interpolator structure. The delay line is zeroed. */
firinterp16x16_handle_t firinterp16x16_init( void * objmem, 
                                             int D, int M, 
                                             const int16_t * restrict h )
{
  firinterp16x16_ptr_t firinterp;
  void *               ptr;
  int16_t *            delayLine;
  int                  delayLen;
  int16_t *            coefBuf;
  int16_t *            coefBank;
  int                  coefNum;
  proc_fxn_t *         procFxn;

  int d, m;

  ASSERT( objmem && (D>1) && (M>0) && h );

  ASSERT(  !(M&3) && IS_ALIGN(h) );

  //
  // Select the processing function, delay line length and coefficients
  // block layout.
  //

  delayLen = (M + 8);
  coefNum  = (M + 4)*D;
  procFxn  = (D == 2 ? &firinterp2_proc :
              D == 3 ? &firinterp3_proc :
              D == 4 ? &firinterp4_proc :
                       &firinterpX_proc );

  //
  // Partition the memory block.
  //
  ptr       = objmem;
  firinterp = (firinterp16x16_ptr_t)(ALIGNED_ADDR( ptr, 4 ));
  ptr       = firinterp + 1;
  delayLine = (int16_t*)(ALIGNED_ADDR( ptr, 8 ));
  ptr       = delayLine + delayLen;
  coefBuf   = (int16_t*)(ALIGNED_ADDR( ptr, 8 ));
  ptr       = coefBuf + coefNum;

  ASSERT( (int8_t*)ptr - (int8_t*)objmem <= (int)firinterp16x16_alloc( D, M ) );

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

} /* firinterp16x16_init() */

/* Put a chunk of input signal into the delay line and compute the filter
 * response. */
void firinterp16x16_process( firinterp16x16_handle_t _firinterp, 
                             int16_t * restrict      y,
                       const int16_t * restrict      x, int N )
{
  firinterp16x16_ptr_t firinterp = (firinterp16x16_ptr_t)_firinterp;

  ASSERT( firinterp && ((firinterp->magic) == MAGIC) && y && x );

  ASSERT(IS_ALIGN(x));

  //
  // Call filter's data processing function. It will store the block of input 
  // samples to the delay line, and compute the filter response. Returns the
  // updated next position pointer into the delay line buffer.
  //

  ASSERT(firinterp->procFxn);

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

} /* firinterp16x16_process() */

/* Data processing function for a factor 2 interpolating FIR filter. */
static int16_t * firinterp2_proc( int16_t * restrict y,
                                  int16_t *          cbegin,
                                  int16_t *          cend,
                                  int16_t * restrict cpos,
                            const int16_t * restrict x,
                            const int16_t * restrict h,
                            int D, int N, int M )
{
        ae_int16x4 * restrict D_wr;
  const ae_int16x4 *          D_rd;
  const ae_int16x4 *          X;
        ae_int16x4 * restrict Y;
  const ae_f16x4   *          C;

  ae_f32x2   d0, d1, d2, d3, d4, d5;
  ae_int16x4 t0, t1, t2;
  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_f64     qq0, qq1, qq2, qq3, qq4, qq5, qq6, qq7;
  ae_int64   u0, u1, u2, u3, u4, u5, u6, u7;
  ae_int32x2 temp0, temp1;
  ae_f16x4   c;
  ae_valign ay;

  int n, m;

  ASSERT(y && cbegin && cend && cpos && x && h && (D == 2) && (N > 0) && (M > 0));

  ASSERT(!(N&7) && !(M&3));

  ASSERT( IS_ALIGN(cbegin) && IS_ALIGN(cend) &&
          IS_ALIGN(cpos  ) && IS_ALIGN(x   ) && IS_ALIGN(h) );

  //
  // Setup pointers and circular delay line buffer.
  //

  D_wr = (      ae_int16x4 *)cpos;
  X    = (const ae_int16x4 *)x;
  Y    = (      ae_int16x4 *)y;
  ay = AE_ZALIGN64();

  WUR_AE_CBEGIN0((uintptr_t)cbegin);
  WUR_AE_CEND0  ((uintptr_t)cend  );

  //
  // Break input signal into 8-samples blocks. For each block, store 8 samples
  // to the delay line buffer and compute 16 samples of interpolating FIR filter
  // response.
  //

  for (n=0; n<(N>>3); n++)
  {
    // Load 8 input samples.
    // Q15
    AE_L16X4_IP(t0, X, +8);
    AE_L16X4_IP(t1, X, +8);

    // Store input samples into the circular delay line.
    // Q15
    AE_S16X4_XC(t0, D_wr, +8);
    AE_S16X4_XC(t1, D_wr, +8);

    //--------------------------------------------------------------------------
    // Bank 0
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q15
    AE_L16X4_XC(t0, D_rd, +8);
    AE_L16X4_XC(t1, D_rd, +8);
    // Q31<-Q15
    d0 = AE_CVT32X2F16_32(t0);
    d1 = AE_CVT32X2F16_10(t0);
    d2 = AE_CVT32X2F16_32(t1);
    d3 = AE_CVT32X2F16_10(t1);

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q15
    AE_L16X4_XC(t2, D_rd, +8);
    // Q31<-Q15
    d4 = AE_CVT32X2F16_32(t2);
    d5 = AE_CVT32X2F16_10(t2);

    // Reset the coefficients pointer. Now it looks at the oldest sample
    // coefficient, bank #0.
    C = (const ae_f16x4 *)h;

    // Load first 4 tap coefficients of bank #0.
    // Q15
    ae_f16x4_loadip(c, C, +8);

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH(q0, q1, d0, d1, c);
    AE_MULFD32X16X2_FIR_HH(q2, q3, d1, d2, c);
    AE_MULFD32X16X2_FIR_HH(q4, q5, d2, d3, c);
    AE_MULFD32X16X2_FIR_HH(q6, q7, d3, d4, c);

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL(q0, q1, d1, d2, c);
    AE_MULAFD32X16X2_FIR_HL(q2, q3, d2, d3, c);
    AE_MULAFD32X16X2_FIR_HL(q4, q5, d3, d4, c);
    AE_MULAFD32X16X2_FIR_HL(q6, q7, d4, d5, c);

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3; 
    d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //

    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q15
      AE_L16X4_XC(t2, D_rd, +8);
      // Q31<-Q15
      d4 = AE_CVT32X2F16_32(t2);
      d5 = AE_CVT32X2F16_10(t2);

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip(c, C, +8);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH(q0, q1, d0, d1, c);
      AE_MULAFD32X16X2_FIR_HH(q2, q3, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HH(q4, q5, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HH(q6, q7, d3, d4, c);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL(q0, q1, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HL(q2, q3, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HL(q4, q5, d3, d4, c);
      AE_MULAFD32X16X2_FIR_HL(q6, q7, d4, d5, c);

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3; 
      d2 = d4; d3 = d5;
    }

    //
    // Scale outputs by the interpolation factor 2.
    //

    // Q16.47 <- Q16.47 + Q16.47 w/ saturation
    u0 = ( q0 );
    u0 = AE_ADD64S( u0, u0 );
    qq0 = ( u0  );

    u1 = ( q1 );
    u1 = AE_ADD64S( u1, u1 );
    qq1 = ( u1  );

    u2 = ( q2 );
    u2 = AE_ADD64S( u2, u2 );
    qq2 = ( u2  );

    u3 = ( q3 );
    u3 = AE_ADD64S( u3, u3 );
    qq3 = ( u3  );

    u4 = ( q4 );
    u4 = AE_ADD64S( u4, u4 );
    qq4 = ( u4  );

    u5 = ( q5 );
    u5 = AE_ADD64S( u5, u5 );
    qq5 = ( u5  );

    u6 = ( q6 );
    u6 = AE_ADD64S( u6, u6 );
    qq6 = ( u6  );

    u7 = ( q7 );
    u7 = AE_ADD64S( u7, u7 );
    qq7 = ( u7  );


    //--------------------------------------------------------------------------
    // Bank 1
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q15
    AE_L16X4_XC(t0, D_rd, +8);
    AE_L16X4_XC(t1, D_rd, +8);
    // Q31<-Q15
    d0 = AE_CVT32X2F16_32(t0);
    d1 = AE_CVT32X2F16_10(t0);
    d2 = AE_CVT32X2F16_32(t1);
    d3 = AE_CVT32X2F16_10(t1);

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q15
    AE_L16X4_XC(t2, D_rd, +8);
    //Q31<-Q15
    d4 = AE_CVT32X2F16_32(t2);
    d5 = AE_CVT32X2F16_10(t2);

    // Load first 4 tap coefficients of bank #1.
    // Q15
    ae_f16x4_loadip(c, C, +8);

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH(q0, q1, d0, d1, c);
    AE_MULFD32X16X2_FIR_HH(q2, q3, d1, d2, c);
    AE_MULFD32X16X2_FIR_HH(q4, q5, d2, d3, c);
    AE_MULFD32X16X2_FIR_HH(q6, q7, d3, d4, c);

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL(q0, q1, d1, d2, c);
    AE_MULAFD32X16X2_FIR_HL(q2, q3, d2, d3, c);
    AE_MULAFD32X16X2_FIR_HL(q4, q5, d3, d4, c);
    AE_MULAFD32X16X2_FIR_HL(q6, q7, d4, d5, c);

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3;
    d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //

    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q15
      AE_L16X4_XC(t2, D_rd, +8);
      // Q31<-Q15
      d4 = AE_CVT32X2F16_32(t2);
      d5 = AE_CVT32X2F16_10(t2);

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip(c, C, +8);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH(q0, q1, d0, d1, c);
      AE_MULAFD32X16X2_FIR_HH(q2, q3, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HH(q4, q5, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HH(q6, q7, d3, d4, c);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL(q0, q1, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HL(q2, q3, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HL(q4, q5, d3, d4, c);
      AE_MULAFD32X16X2_FIR_HL(q6, q7, d4, d5, c);

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3;
      d2 = d4; d3 = d5;
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

    // Q15 <- Q16.47 - 16 w/ rounding and saturation
    temp0 = AE_TRUNCA32X2F64S(qq0, q0, 16);
    temp1 = AE_TRUNCA32X2F64S(qq1, q1, 16);
    AE_SA16X4_IP(AE_ROUND16X4F32SASYM(temp0, temp1), ay, Y); 
    temp0 = AE_TRUNCA32X2F64S(qq2, q2, 16);
    temp1 = AE_TRUNCA32X2F64S(qq3, q3, 16);
    AE_SA16X4_IP(AE_ROUND16X4F32SASYM(temp0, temp1), ay, Y); 
    temp0 = AE_TRUNCA32X2F64S(qq4, q4, 16);
    temp1 = AE_TRUNCA32X2F64S(qq5, q5, 16);
    AE_SA16X4_IP(AE_ROUND16X4F32SASYM(temp0, temp1), ay, Y); 
    temp0 = AE_TRUNCA32X2F64S(qq6, q6, 16);
    temp1 = AE_TRUNCA32X2F64S(qq7, q7, 16);
    AE_SA16X4_IP(AE_ROUND16X4F32SASYM(temp0, temp1), ay, Y); 
  }
  AE_SA64POS_FP(ay, Y);
  return (int16_t*)D_wr;

} // firinterp2_proc()

/* Data processing function for a factor 3 interpolating FIR filter. */
static int16_t * firinterp3_proc( int16_t * restrict y,
                                  int16_t *          cbegin,
                                  int16_t *          cend,
                                  int16_t * restrict cpos,
                            const int16_t * restrict x,
                            const int16_t * restrict h,
                            int D, int N, int M )
{
        ae_int16x4 * restrict D_wr;
  const ae_int16x4 *          D_rd;
  const ae_int16x4 *          X;
        ae_int16x4 * restrict Y;
  const ae_f16x4   *          C;

  ae_f32x2   d0, d1, d2, d3, d4, d5;
  ae_int16x4 t0, t1, t2;
  ae_f64     q0,   q1,   q2,   q3,   q4,   q5,   q6,   q7;
  ae_f64     qq0,  qq1,  qq2,  qq3,  qq4,  qq5,  qq6,  qq7;
  ae_f64     qqq0, qqq1, qqq2, qqq3, qqq4, qqq5, qqq6, qqq7;
  ae_int64   u0, u1, u2, u3, u4, u5, u6, u7;
  ae_int64   v0, v1, v2, v3, v4, v5, v6, v7;
  ae_f16x4   c;
  ae_int32x2 temp0, temp1;
  ae_valign ay;

  int n, m;

  ASSERT( y && cbegin && cend && cpos && x && h && (D == 3) && (N > 0) && (M > 0) );

  ASSERT(!(N&7) && !(M&3));

  ASSERT( IS_ALIGN(cbegin) && IS_ALIGN(cend) &&
          IS_ALIGN(cpos  ) && IS_ALIGN(x   ) && IS_ALIGN(h) );

  //
  // Setup pointers and circular delay line buffer.
  //

  D_wr = (      ae_int16x4 *)cpos;
  X    = (const ae_int16x4 *)x;
  Y    = (      ae_int16x4 *)y;

  WUR_AE_CBEGIN0( (uintptr_t)cbegin );
  WUR_AE_CEND0  ( (uintptr_t)cend   );
  ay = AE_ZALIGN64();

  //
  // Break input signal into 8-samples blocks. For each block, store 8 samples
  // to the delay line buffer and compute 16 samples of interpolating FIR filter
  // response.
  //

  for (n=0; n<(N>>3); n++)
  {
    // Load 8 input samples.
    // Q15
    AE_L16X4_IP(t0, X, +8);
    AE_L16X4_IP(t1, X, +8);

    // Store input samples into the circular delay line.
    // Q15
    AE_S16X4_XC(t0, D_wr, +8);
    AE_S16X4_XC(t1, D_wr, +8);

    //--------------------------------------------------------------------------
    // Bank 0
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q15
    AE_L16X4_XC(t0, D_rd, +8);
    AE_L16X4_XC(t1, D_rd, +8);
    // Q31<-Q15
    d0 = AE_CVT32X2F16_32(t0);
    d1 = AE_CVT32X2F16_10(t0);
    d2 = AE_CVT32X2F16_32(t1);
    d3 = AE_CVT32X2F16_10(t1);

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q15
    AE_L16X4_XC( t2, D_rd, +8 );
    // Q31<-Q15
    d4 = AE_CVT32X2F16_32(t2);
    d5 = AE_CVT32X2F16_10(t2);

    // Reset the coefficients pointer. Now it looks at the oldest sample
    // coefficient, bank #0.
    C = (const ae_f16x4 *)h;

    // Load first 4 tap coefficients of bank #0.
    // Q15
    ae_f16x4_loadip(c, C, +8);

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH(q0, q1, d0, d1, c);
    AE_MULFD32X16X2_FIR_HH(q2, q3, d1, d2, c);
    AE_MULFD32X16X2_FIR_HH(q4, q5, d2, d3, c);
    AE_MULFD32X16X2_FIR_HH(q6, q7, d3, d4, c);

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL(q0, q1, d1, d2, c);
    AE_MULAFD32X16X2_FIR_HL(q2, q3, d2, d3, c);
    AE_MULAFD32X16X2_FIR_HL(q4, q5, d3, d4, c);
    AE_MULAFD32X16X2_FIR_HL(q6, q7, d4, d5, c);

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3;
    d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //

    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q15
      AE_L16X4_XC(t2, D_rd, +8);
      // Q31<-Q15
      d4 = AE_CVT32X2F16_32(t2);
      d5 = AE_CVT32X2F16_10(t2);

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip(c, C, +8);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH(q0, q1, d0, d1, c);
      AE_MULAFD32X16X2_FIR_HH(q2, q3, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HH(q4, q5, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HH(q6, q7, d3, d4, c);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL(q0, q1, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HL(q2, q3, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HL(q4, q5, d3, d4, c);
      AE_MULAFD32X16X2_FIR_HL(q6, q7, d4, d5, c);

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3; 
      d2 = d4; d3 = d5;
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

    //--------------------------------------------------------------------------
    // Bank 1
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q15
    AE_L16X4_XC(t0, D_rd, +8);
    AE_L16X4_XC(t1, D_rd, +8);
    // Q31<-Q15
    d0 = AE_CVT32X2F16_32(t0);
    d1 = AE_CVT32X2F16_10(t0);
    d2 = AE_CVT32X2F16_32(t1);
    d3 = AE_CVT32X2F16_10(t1);

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q15
    AE_L16X4_XC(t2, D_rd, +8);
    // Q31<-Q15
    d4 = AE_CVT32X2F16_32(t2);
    d5 = AE_CVT32X2F16_10(t2);

    // Load first 4 tap coefficients of bank #1.
    // Q15
    ae_f16x4_loadip(c, C, +8);

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH(qq0, qq1, d0, d1, c);
    AE_MULFD32X16X2_FIR_HH(qq2, qq3, d1, d2, c);
    AE_MULFD32X16X2_FIR_HH(qq4, qq5, d2, d3, c);
    AE_MULFD32X16X2_FIR_HH(qq6, qq7, d3, d4, c);

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL(qq0, qq1, d1, d2, c);
    AE_MULAFD32X16X2_FIR_HL(qq2, qq3, d2, d3, c);
    AE_MULAFD32X16X2_FIR_HL(qq4, qq5, d3, d4, c);
    AE_MULAFD32X16X2_FIR_HL(qq6, qq7, d4, d5, c);

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3;
    d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //

    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q15
      AE_L16X4_XC(t2, D_rd, +8);
      // Q31<-Q15
      d4 = AE_CVT32X2F16_32(t2);
      d5 = AE_CVT32X2F16_10(t2);

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip(c, C, +8);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH(qq0, qq1, d0, d1, c);
      AE_MULAFD32X16X2_FIR_HH(qq2, qq3, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HH(qq4, qq5, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HH(qq6, qq7, d3, d4, c);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL(qq0, qq1, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HL(qq2, qq3, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HL(qq4, qq5, d3, d4, c);
      AE_MULAFD32X16X2_FIR_HL(qq6, qq7, d4, d5, c);

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3; 
      d2 = d4; d3 = d5;
    }

    //
    // Scale outputs by the interpolation factor 3.
    //

    // Q16.47 <- Q16.47 + Q16.47 + Q16.47 w/ saturation
    u0 = ( qq0 );
    v0 = AE_ADD64S( u0, u0 );
    v0 = AE_ADD64S( u0, v0 );
    qq0 = ( v0 );

    u1 = ( qq1 );
    v1 = AE_ADD64S( u1, u1 );
    v1 = AE_ADD64S( u1, v1 );
    qq1 = ( v1 );

    u2 = ( qq2 );
    v2 = AE_ADD64S( u2, u2 );
    v2 = AE_ADD64S( u2, v2 );
    qq2 = ( v2 );

    u3 = ( qq3 );
    v3 = AE_ADD64S( u3, u3 );
    v3 = AE_ADD64S( u3, v3 );
    qq3 = ( v3 );

    u4 = ( qq4 );
    v4 = AE_ADD64S( u4, u4 );
    v4 = AE_ADD64S( u4, v4 );
    qq4 = ( v4 );

    u5 = ( qq5 );
    v5 = AE_ADD64S( u5, u5 );
    v5 = AE_ADD64S( u5, v5 );
    qq5 = ( v5 );

    u6 = ( qq6 );
    v6 = AE_ADD64S( u6, u6 );
    v6 = AE_ADD64S( u6, v6 );
    qq6 = ( v6 );

    u7 = ( qq7 );
    v7 = AE_ADD64S( u7, u7 );
    v7 = AE_ADD64S( u7, v7 );
    qq7 = ( v7 );

    //--------------------------------------------------------------------------
    // Bank 2
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q15
    AE_L16X4_XC(t0, D_rd, +8);
    AE_L16X4_XC(t1, D_rd, +8);
    // Q31<-Q15
    d0 = AE_CVT32X2F16_32(t0);
    d1 = AE_CVT32X2F16_10(t0);
    d2 = AE_CVT32X2F16_32(t1);
    d3 = AE_CVT32X2F16_10(t1);

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q15
    AE_L16X4_XC(t2, D_rd, +8);
    // Q31<-Q15
    d4 = AE_CVT32X2F16_32(t2);
    d5 = AE_CVT32X2F16_10(t2);

    // Load first 4 tap coefficients of bank #2.
    // Q15
    ae_f16x4_loadip(c, C, +8);

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH(qqq0, qqq1, d0, d1, c);
    AE_MULFD32X16X2_FIR_HH(qqq2, qqq3, d1, d2, c);
    AE_MULFD32X16X2_FIR_HH(qqq4, qqq5, d2, d3, c);
    AE_MULFD32X16X2_FIR_HH(qqq6, qqq7, d3, d4, c);

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL(qqq0, qqq1, d1, d2, c);
    AE_MULAFD32X16X2_FIR_HL(qqq2, qqq3, d2, d3, c);
    AE_MULAFD32X16X2_FIR_HL(qqq4, qqq5, d3, d4, c);
    AE_MULAFD32X16X2_FIR_HL(qqq6, qqq7, d4, d5, c);

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3;
    d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //

    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q15
      AE_L16X4_XC(t2, D_rd, +8);
      // Q31<-Q15
      d4 = AE_CVT32X2F16_32(t2);
      d5 = AE_CVT32X2F16_10(t2);

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip(c, C, +8);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH(qqq0, qqq1, d0, d1, c);
      AE_MULAFD32X16X2_FIR_HH(qqq2, qqq3, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HH(qqq4, qqq5, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HH(qqq6, qqq7, d3, d4, c);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL(qqq0, qqq1, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HL(qqq2, qqq3, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HL(qqq4, qqq5, d3, d4, c);
      AE_MULAFD32X16X2_FIR_HL(qqq6, qqq7, d4, d5, c);

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3;
      d2 = d4; d3 = d5;
    }

    //
    // Scale outputs by the interpolation factor 3.
    //

    // Q16.47 <- Q16.47 + Q16.47 + Q16.47 w/ saturation
    u0 = ( qqq0 );
    v0 = AE_ADD64S( u0, u0 );
    v0 = AE_ADD64S( u0, v0 );
    qqq0 = ( v0 );

    u1 = ( qqq1 );
    v1 = AE_ADD64S( u1, u1 );
    v1 = AE_ADD64S( u1, v1 );
    qqq1 = ( v1 );

    u2 = ( qqq2 );
    v2 = AE_ADD64S( u2, u2 );
    v2 = AE_ADD64S( u2, v2 );
    qqq2 = ( v2 );

    u3 = ( qqq3 );
    v3 = AE_ADD64S( u3, u3 );
    v3 = AE_ADD64S( u3, v3 );
    qqq3 = ( v3 );

    u4 = ( qqq4 );
    v4 = AE_ADD64S( u4, u4 );
    v4 = AE_ADD64S( u4, v4 );
    qqq4 = ( v4 );

    u5 = ( qqq5 );
    v5 = AE_ADD64S( u5, u5 );
    v5 = AE_ADD64S( u5, v5 );
    qqq5 = ( v5 );

    u6 = ( qqq6 );
    v6 = AE_ADD64S( u6, u6 );
    v6 = AE_ADD64S( u6, v6 );
    qqq6 = ( v6 );

    u7 = ( qqq7 );
    v7 = AE_ADD64S( u7, u7 );
    v7 = AE_ADD64S( u7, v7 );
    qqq7 = ( v7 );

    //
    // Store 24 outputs to the output array with 3-samples stride.
    //

    // Q15 <- Q16.47 - 16 w/ rounding and saturation
    temp0 = AE_TRUNCA32X2F64S(q0, qq0, 16);
    temp1 = AE_TRUNCA32X2F64S(qqq0, q1, 16);
    AE_SA16X4_IP(AE_ROUND16X4F32SASYM(temp0, temp1), ay, Y); 
    temp0 = AE_TRUNCA32X2F64S(qq1, qqq1, 16);
    temp1 = AE_TRUNCA32X2F64S(q2, qq2, 16);
    AE_SA16X4_IP(AE_ROUND16X4F32SASYM(temp0, temp1), ay, Y); 
    temp0 = AE_TRUNCA32X2F64S(qqq2, q3, 16);
    temp1 = AE_TRUNCA32X2F64S(qq3, qqq3, 16);
    AE_SA16X4_IP(AE_ROUND16X4F32SASYM(temp0, temp1), ay, Y); 
    temp0 = AE_TRUNCA32X2F64S(q4, qq4, 16);
    temp1 = AE_TRUNCA32X2F64S(qqq4, q5, 16);
    AE_SA16X4_IP(AE_ROUND16X4F32SASYM(temp0, temp1), ay, Y); 
    temp0 = AE_TRUNCA32X2F64S(qq5, qqq5, 16);
    temp1 = AE_TRUNCA32X2F64S(q6, qq6, 16);
    AE_SA16X4_IP(AE_ROUND16X4F32SASYM(temp0, temp1), ay, Y); 
    temp0 = AE_TRUNCA32X2F64S(qqq6, q7, 16);
    temp1 = AE_TRUNCA32X2F64S(qq7, qqq7, 16);
    AE_SA16X4_IP(AE_ROUND16X4F32SASYM(temp0, temp1), ay, Y);
  }
  AE_SA64POS_FP(ay, Y);
  return ( (int16_t*)D_wr );
} // firinterp3_proc()

/* Data processing function for a factor 4 interpolating FIR filter. */
static int16_t * firinterp4_proc( int16_t * restrict y,
                                  int16_t *          cbegin,
                                  int16_t *          cend,
                                  int16_t * restrict cpos,
                            const int16_t * restrict x,
                            const int16_t * restrict h,
                            int D, int N, int M )
{
        ae_int16x4 * restrict D_wr;
  const ae_int16x4 *          D_rd;
  const ae_int16x4 *          X;
        ae_int16   * restrict Y;
  const ae_f16x4   *          C;

  ae_f32x2   d0, d1, d2, d3, d4, d5;
  ae_int16x4 t0, t1, t2;
  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_int64   u0, u1, u2, u3, u4, u5, u6, u7;
  ae_f16x4   c;
  ae_int32x2 temp;

  int n, m;

  ASSERT( y && cbegin && cend && cpos && x && h && (D == 4) && (N > 0) && (M > 0) );

  ASSERT( !(N&7) && !(M&3) );

  ASSERT( IS_ALIGN(cbegin) && IS_ALIGN(cend) && 
          IS_ALIGN(cpos  ) && IS_ALIGN(x   ) && IS_ALIGN(h) );

  //
  // Setup pointers and circular delay line buffer.
  //

  D_wr = (      ae_int16x4 *)cpos;
  X    = (const ae_int16x4 *)x;
  Y    = (      ae_int16   *)y;

  WUR_AE_CBEGIN0((uintptr_t)cbegin);
  WUR_AE_CEND0  ((uintptr_t)cend  );

  //
  // Break input signal into 8-samples blocks. For each block, store 8 samples
  // to the delay line buffer and compute 16 samples of interpolating FIR filter
  // response.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 8 input samples.
    // Q15
    AE_L16X4_IP( t0, X, +8 );
    AE_L16X4_IP( t1, X, +8 );

    // Store input samples into the circular delay line.
    // Q15
    AE_S16X4_XC( t0, D_wr, +8 );
    AE_S16X4_XC( t1, D_wr, +8 );

    //--------------------------------------------------------------------------
    // Bank 0
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q15
    AE_L16X4_XC(t0, D_rd, +8);
    AE_L16X4_XC(t1, D_rd, +8);
    //Q31<-Q15
    d0 = AE_CVT32X2F16_32(t0);
    d1 = AE_CVT32X2F16_10(t0);
    d2 = AE_CVT32X2F16_32(t1);
    d3 = AE_CVT32X2F16_10(t1);

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q15
    AE_L16X4_XC(t2, D_rd, +8);
    //Q31<-Q15
    d4 = AE_CVT32X2F16_32(t2);
    d5 = AE_CVT32X2F16_10(t2);

    // Reset the coefficients pointer. Now it looks at the oldest sample
    // coefficient, bank #0.
    C = (const ae_f16x4 *)h;

    // Load first 4 tap coefficients of bank #0.
    // Q15
    ae_f16x4_loadip( c, C, +8 );

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH(q0, q1, d0, d1, c );
    AE_MULFD32X16X2_FIR_HH(q2, q3, d1, d2, c );
    AE_MULFD32X16X2_FIR_HH(q4, q5, d2, d3, c );
    AE_MULFD32X16X2_FIR_HH(q6, q7, d3, d4, c );

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL(q0, q1, d1, d2, c);
    AE_MULAFD32X16X2_FIR_HL(q2, q3, d2, d3, c);
    AE_MULAFD32X16X2_FIR_HL(q4, q5, d3, d4, c);
    AE_MULAFD32X16X2_FIR_HL(q6, q7, d4, d5, c);

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3;
    d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //

    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q15
      AE_L16X4_XC(t2, D_rd, +8);
      // Q31<-Q15
      d4 = AE_CVT32X2F16_32(t2);
      d5 = AE_CVT32X2F16_10(t2);

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip(c, C, +8);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH(q0, q1, d0, d1, c);
      AE_MULAFD32X16X2_FIR_HH(q2, q3, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HH(q4, q5, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HH(q6, q7, d3, d4, c);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL(q0, q1, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HL(q2, q3, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HL(q4, q5, d3, d4, c);
      AE_MULAFD32X16X2_FIR_HL(q6, q7, d4, d5, c);

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3; 
      d2 = d4; d3 = d5;
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

    // Q15 <- Q16.47 - 16 w/ rounding and saturation
    temp = AE_TRUNCA32X2F64S(0, q0, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q1, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q2, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q3, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q4, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q5, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q6, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q7, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, -27*2);

    //--------------------------------------------------------------------------
    // Bank 1
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q15
    AE_L16X4_XC(t0, D_rd, +8);
    AE_L16X4_XC(t1, D_rd, +8);
    // Q31<-Q15
    d0 = AE_CVT32X2F16_32(t0);
    d1 = AE_CVT32X2F16_10(t0);
    d2 = AE_CVT32X2F16_32(t1);
    d3 = AE_CVT32X2F16_10(t1);

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q15
    AE_L16X4_XC( t2, D_rd, +8 );
    // Q31<-Q15
    d4 = AE_CVT32X2F16_32(t2);
    d5 = AE_CVT32X2F16_10(t2);

    // Load first 4 tap coefficients of bank #1.
    // Q15
    ae_f16x4_loadip(c, C, +8);

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH(q0, q1, d0, d1, c);
    AE_MULFD32X16X2_FIR_HH(q2, q3, d1, d2, c);
    AE_MULFD32X16X2_FIR_HH(q4, q5, d2, d3, c);
    AE_MULFD32X16X2_FIR_HH(q6, q7, d3, d4, c);

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL(q0, q1, d1, d2, c);
    AE_MULAFD32X16X2_FIR_HL(q2, q3, d2, d3, c);
    AE_MULAFD32X16X2_FIR_HL(q4, q5, d3, d4, c);
    AE_MULAFD32X16X2_FIR_HL(q6, q7, d4, d5, c);

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3; 
    d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //

    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q15
      AE_L16X4_XC( t2, D_rd, +8 );
      // Q31<-Q15
      d4 = AE_CVT32X2F16_32(t2);
      d5 = AE_CVT32X2F16_10(t2);

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip(c, C, +8);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH(q0, q1, d0, d1, c);
      AE_MULAFD32X16X2_FIR_HH(q2, q3, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HH(q4, q5, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HH(q6, q7, d3, d4, c);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL(q0, q1, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HL(q2, q3, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HL(q4, q5, d3, d4, c);
      AE_MULAFD32X16X2_FIR_HL(q6, q7, d4, d5, c);

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3; 
      d2 = d4; d3 = d5;
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

    // Q15 <- Q16.47 - 16 w/ rounding and saturation
    temp = AE_TRUNCA32X2F64S(0, q0, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q1, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q2, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q3, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q4, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q5, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q6, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q7, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, -27*2);

    //--------------------------------------------------------------------------
    // Bank 2
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q15
    AE_L16X4_XC(t0, D_rd, +8);
    AE_L16X4_XC(t1, D_rd, +8);
    // Q31<-Q15
    d0 = AE_CVT32X2F16_32(t0);
    d1 = AE_CVT32X2F16_10(t0);
    d2 = AE_CVT32X2F16_32(t1);
    d3 = AE_CVT32X2F16_10(t1);

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q15
    AE_L16X4_XC(t2, D_rd, +8);
    // Q31<-Q15
    d4 = AE_CVT32X2F16_32(t2);
    d5 = AE_CVT32X2F16_10(t2);

    // Load first 4 tap coefficients of bank #2.
    // Q15
    ae_f16x4_loadip(c, C, +8);

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH(q0, q1, d0, d1, c);
    AE_MULFD32X16X2_FIR_HH(q2, q3, d1, d2, c);
    AE_MULFD32X16X2_FIR_HH(q4, q5, d2, d3, c);
    AE_MULFD32X16X2_FIR_HH(q6, q7, d3, d4, c);

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL(q0, q1, d1, d2, c);
    AE_MULAFD32X16X2_FIR_HL(q2, q3, d2, d3, c);
    AE_MULAFD32X16X2_FIR_HL(q4, q5, d3, d4, c);
    AE_MULAFD32X16X2_FIR_HL(q6, q7, d4, d5, c);

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3;
    d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //

    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q15
      AE_L16X4_XC( t2, D_rd, +8 );
      // Q31<-Q15
      d4 = AE_CVT32X2F16_32(t2);
      d5 = AE_CVT32X2F16_10(t2);

      // Load 4 tap coefficients.
      // Q15
      ae_f16x4_loadip(c, C, +8);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH(q0, q1, d0, d1, c);
      AE_MULAFD32X16X2_FIR_HH(q2, q3, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HH(q4, q5, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HH(q6, q7, d3, d4, c);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL(q0, q1, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HL(q2, q3, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HL(q4, q5, d3, d4, c);
      AE_MULAFD32X16X2_FIR_HL(q6, q7, d4, d5, c);

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3;
      d2 = d4; d3 = d5;
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

    // Q15 <- Q16.47 - 16 w/ rounding and saturation
    temp = AE_TRUNCA32X2F64S(0, q0, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q1, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q2, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q3, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q4, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q5, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q6, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q7, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, -27*2);

    //--------------------------------------------------------------------------
    // Bank 3
    //

    // Start reading the delay line from the oldest sample, M+8 samples back
    // from the newest one.
    D_rd = D_wr;

    // Load 8 oldest samples.
    // Q15
    AE_L16X4_XC(t0, D_rd, +8);
    AE_L16X4_XC(t1, D_rd, +8);
    // Q31<-Q15
    d0 = AE_CVT32X2F16_32(t0);
    d1 = AE_CVT32X2F16_10(t0);
    d2 = AE_CVT32X2F16_32(t1);
    d3 = AE_CVT32X2F16_10(t1);

    //
    // Inner loop prologue: process first 4 taps for 8 accumulators.
    //

    // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
    // registers to perform 4 taps processing for 8 accumulators.
    // Q15
    AE_L16X4_XC( t2, D_rd, +8 );
    // Q31<-Q15
    d4 = AE_CVT32X2F16_32(t2);
    d5 = AE_CVT32X2F16_10(t2);

    // Load first 4 tap coefficients of bank #3.
    // Q15
    ae_f16x4_loadip(c, C, +8);

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH(q0, q1, d0, d1, c);
    AE_MULFD32X16X2_FIR_HH(q2, q3, d1, d2, c);
    AE_MULFD32X16X2_FIR_HH(q4, q5, d2, d3, c);
    AE_MULFD32X16X2_FIR_HH(q6, q7, d3, d4, c);

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL(q0, q1, d1, d2, c);
    AE_MULAFD32X16X2_FIR_HL(q2, q3, d2, d3, c);
    AE_MULAFD32X16X2_FIR_HL(q4, q5, d3, d4, c);
    AE_MULAFD32X16X2_FIR_HL(q6, q7, d4, d5, c);

    // 4 taps are done. Move 4 samples out of the registers file.
    d0 = d2; d1 = d3; 
    d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
    // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
    // coefficients appended to the impulse response to avoid a 1-sample
    // response delay.
    //

    for ( m=0; m<(M>>2); m++ )
    {
      // Load 4 samples.
      // Q15
      AE_L16X4_XC(t2, D_rd, +8);
      // Q31<-Q15
      d4 = AE_CVT32X2F16_32(t2);
      d5 = AE_CVT32X2F16_10(t2);

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
      d0 = d2; d1 = d3;
      d2 = d4; d3 = d5;
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

    // Q15 <- Q16.47 - 16 w/ rounding and saturation
    temp = AE_TRUNCA32X2F64S(0, q0, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q1, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q2, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q3, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q4, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q5, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q6, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 4*2);
    temp = AE_TRUNCA32X2F64S(0, q7, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, 1*2);
  }
  return (int16_t*)D_wr;
} // firinterp4_proc()

/* Data processing function for a generic interpolating FIR filter with an
 * arbitrary interpolation factor. */
static int16_t * firinterpX_proc( int16_t * restrict y,
                                  int16_t *          cbegin,
                                  int16_t *          cend,
                                  int16_t * restrict cpos,
                            const int16_t * restrict x,
                            const int16_t * restrict h,
                            int D, int N, int M )
{
        ae_int16x4 * restrict D_wr;
  const ae_int16x4 *          D_rd;
  const ae_int16x4 *          X;
        ae_int16   * restrict Y;
  const ae_f16x4   *          C;

  ae_f32x2   d0, d1, d2, d3, d4, d5;
  ae_int16x4 t0, t1, t2;
  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_int64   u0, u1, u2, u3, u4, u5, u6, u7;
  ae_f16x4   c;
  ae_int32x2 temp;

  int        g_exp;
  ae_p24x2s g_frac;

  int d, n, m;

  ASSERT( y && cbegin && cend && cpos && x && h && (D > 1) && (N > 0) && (M > 0) );

  ASSERT( !( N & 7 ) && !( M & 3 ) );

  ASSERT( IS_ALIGN(cbegin) && IS_ALIGN(cend) &&
          IS_ALIGN(cpos  ) && IS_ALIGN(x   ) && IS_ALIGN(h) );

  //
  // Calculate the scaling exponent and fraction.
  //

  // Q(15+8-g_exp) <- Q0 + 15 + 8 - g_exp
  g_exp  = 31 - AE_NSAZ32_L(D);
  {
    int q_frac_int =  D << ( 15 + 8 - g_exp );
    g_frac = *(ae_p24s *) &q_frac_int;
  }

  WUR_AE_SAR( g_exp );

  //
  // Setup pointers and circular delay line buffer.
  //

  D_wr = (      ae_int16x4 *)cpos;
  X    = (const ae_int16x4 *)x;
  Y    = (      ae_int16   *)y;

  WUR_AE_CBEGIN0((uintptr_t)cbegin);
  WUR_AE_CEND0  ((uintptr_t)cend  );

  //
  // Break input signal into 8-samples blocks. For each block, store 8 samples
  // to the delay line buffer and compute 8 outputs per filter bank,
  // in all 8*D samples
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 8 input samples.
    // Q15
    AE_L16X4_IP(t0, X, +8);
    AE_L16X4_IP(t1, X, +8);

    // Store input samples into the circular delay line.
    // Q15
    AE_S16X4_XC(t0, D_wr, +8);
    AE_S16X4_XC(t1, D_wr, +8);

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
      // Q15
      AE_L16X4_XC(t0, D_rd, +8);
      AE_L16X4_XC(t1, D_rd, +8);
      //Q31<-Q15
      d0 = AE_CVT32X2F16_32(t0);
      d1 = AE_CVT32X2F16_10(t0);
      d2 = AE_CVT32X2F16_32(t1);
      d3 = AE_CVT32X2F16_10(t1);

      //
      // Inner loop prologue: process first 4 taps for 8 accumulators.
      //

      // Load 4 more samples. Altogether we need 12 samples residing in 6 AE
      // registers to perform 4 taps processing for 8 accumulators.
      // Q15
      AE_L16X4_XC(t2, D_rd, +8);
      //Q31<-Q15
      d4 = AE_CVT32X2F16_32(t2);
      d5 = AE_CVT32X2F16_10(t2);

      // Load first 4 tap coefficients of bank #d.
      // Q15
      ae_f16x4_loadip(c, C, +8);

      // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULFD32X16X2_FIR_HH(q0, q1, d0, d1, c);
      AE_MULFD32X16X2_FIR_HH(q2, q3, d1, d2, c);
      AE_MULFD32X16X2_FIR_HH(q4, q5, d2, d3, c);
      AE_MULFD32X16X2_FIR_HH(q6, q7, d3, d4, c);

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL(q0, q1, d1, d2, c);
      AE_MULAFD32X16X2_FIR_HL(q2, q3, d2, d3, c);
      AE_MULAFD32X16X2_FIR_HL(q4, q5, d3, d4, c);
      AE_MULAFD32X16X2_FIR_HL(q6, q7, d4, d5, c);

      // 4 taps are done. Move 4 samples out of the registers file.
      d0 = d2; d1 = d3;
      d2 = d4; d3 = d5;

      //
      // Inner loop kernel: process 4 taps for 8 accumulators. Actually we 
      // perform M+4 MACs for each accumulator. 4 extra MACs fall on zero
      // coefficients appended to the impulse response to avoid a 1-sample
      // response delay.
      //

      for ( m=0; m<(M>>2); m++ )
      {
        // Load 4 samples.
        // Q15
        AE_L16X4_XC(t2, D_rd, +8);
        // Q31<-Q15
        d4 = AE_CVT32X2F16_32(t2);
        d5 = AE_CVT32X2F16_10(t2);

        // Load 4 tap coefficients.
        // Q15
        ae_f16x4_loadip(c, C, +8);

        // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
        AE_MULAFD32X16X2_FIR_HH(q0, q1, d0, d1, c);
        AE_MULAFD32X16X2_FIR_HH(q2, q3, d1, d2, c);
        AE_MULAFD32X16X2_FIR_HH(q4, q5, d2, d3, c);
        AE_MULAFD32X16X2_FIR_HH(q6, q7, d3, d4, c);

        // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
        AE_MULAFD32X16X2_FIR_HL(q0, q1, d1, d2, c);
        AE_MULAFD32X16X2_FIR_HL(q2, q3, d2, d3, c);
        AE_MULAFD32X16X2_FIR_HL(q4, q5, d3, d4, c);
        AE_MULAFD32X16X2_FIR_HL(q6, q7, d4, d5, c);

        // 4 taps are done. Move 4 samples out of the registers file.
        d0 = d2; d1 = d3; 
        d2 = d4; d3 = d5;
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
      temp = AE_TRUNCA32X2F64S(0, q0, 16);
      AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, D*2);
      temp = AE_TRUNCA32X2F64S(0, q1, 16);
      AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, D*2);
      temp = AE_TRUNCA32X2F64S(0, q2, 16);
      AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, D*2);
      temp = AE_TRUNCA32X2F64S(0, q3, 16);
      AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, D*2);
      temp = AE_TRUNCA32X2F64S(0, q4, 16);
      AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, D*2);
      temp = AE_TRUNCA32X2F64S(0, q5, 16);
      AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, D*2);
      temp = AE_TRUNCA32X2F64S(0, q6, 16);
      AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, D*2);
      temp = AE_TRUNCA32X2F64S(0, q7, 16);
      AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp), Y, (1-7*D)*2);
    }
    Y = (ae_int16*)( (uintptr_t)Y + 7*D*2 );
  }
  return (int16_t*)D_wr;
} // firinterpX_proc()

#else
/* code for Quad MAC option */

/* Data processing function of a particular interpolating filter */
typedef int (proc_fxn_t)( int16_t * restrict y,
                          int16_t * delayLine,
                          int       delayLen,
                    const int16_t * restrict x,
                    const int16_t * restrict h,
                          int wrIx, int D, int N, int M );

proc_fxn_t firinterp16x16_D2_proc;
proc_fxn_t firinterp16x16_D3_proc;
proc_fxn_t firinterp16x16_D4_proc;
proc_fxn_t firinterp16x16_DX_proc;

/* Interpolator instance structure. */
typedef struct tag_firinterpf_t
{
    uint32_t        magic;     // Instance pointer validation number
    int             D;         // Interpolation factor
    int             M;         // Number of filter coefficients
    proc_fxn_t    * procFxn;   // Filter data processing function
    const int16_t * coef;      // Filter coefficients
    int16_t *       delayLine; // Delay line for complex samples
    int             delayLen;  // Delay line length, in complex samples
    int             wrIx;      // Index of the oldest sample

} firinterp16x16_t, *firinterp16x16_ptr_t;
/* Calculate the memory block size for an interpolator with given
* attributes. */
size_t firinterp16x16_alloc(int D, int M)
{
    int delayLen, coefNum;
    NASSERT( (D > 1) && (M > 0));

    NASSERT(!(M & 3));
    if (D > 4)
        delayLen = (M + 8);
    else
        delayLen = (M + 4);

    coefNum = (M + 4)*D;

    return (ALIGNED_SIZE(sizeof(firinterp16x16_t), 4)
        + // Delay line
        // Allocated memory for delay line has increased by 2 samples
        // to avoid memory bank conflicts and get the best performance
        ALIGNED_SIZE((delayLen + 2)*sz_i16, 8)
        + // Coefficients
        ALIGNED_SIZE(coefNum*sz_i16, 8));
} // firinterp16x16_alloc()

/* Initialize the interpolator structure. The delay line is zeroed. */
firinterp16x16_handle_t firinterp16x16_init( void * objmem,
                                             int D, int M,
                                       const int16_t * restrict h)
{
    firinterp16x16_ptr_t firinterp;
    void            *    ptr;
    int16_t         *    delLine;
    int                  delLen;
    int16_t         *    coef;
    int16_t         *    coefB;
    int d, m;
    proc_fxn_t      *    procFxn;
    NASSERT(objmem && D > 1 && M > 0 && h);

    NASSERT(!(M & 3) && IS_ALIGN(h));

    if (D > 4)
        delLen = (M + 8);
    else
        delLen = (M + 4);

    procFxn = ( D == 2 ? &firinterp16x16_D2_proc :
                D == 3 ? &firinterp16x16_D3_proc :
                D == 4 ? &firinterp16x16_D4_proc :
                         &firinterp16x16_DX_proc);

    //
    // Partition the memory block.
    //
    ptr         = objmem;
    firinterp   = (firinterp16x16_ptr_t)ALIGNED_ADDR(ptr, 4);
    ptr         = firinterp + 1;
    delLine     = (int16_t*)ALIGNED_ADDR(ptr, 8);
    // Allocated memory for delay line has increased by 2 samples
    // to avoid memory bank conflicts and get the best performance
    ptr         = delLine + delLen + 2;
    coef        = (int16_t*)ALIGNED_ADDR(ptr, 8);
    ptr         = coef + M*D;

    NASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)firinterp16x16_alloc(D, M));


    //
    // Break the impulse response into D coefficients banks and copy them in
    // reverted order.
    //

    for (d = 0; d<D; d++)
    {
        coefB = coef + d*(M + 4);

        // To avoid a 1-sample delay, insert a zero coefficient that will match the
        // oldest sample in the delay line. To keep the filter length a multiple of
        // 4, append 3 zeros after the last coefficient.
        coefB[0] = coefB[M + 1] = coefB[M + 2] = coefB[M + 3] = 0;

        // Copy bank's coefficients in reverted order.
        for (m = 1; m<M + 1; m++)
        {
            coefB[m] = h[D*(M - m) + d];
        }
    }

    //
    // Zero the delay line.
    //

    for (m = 0; m<delLen; m++)
    {
        delLine[m] = 0;
    }

    //
    // Initialize the interpolator instance.
    //

    firinterp->magic      = MAGIC;
    firinterp->D          = D;
    firinterp->M          = M;
    firinterp->coef       = coef;
    firinterp->procFxn    = procFxn;
    firinterp->delayLine  = delLine;
    firinterp->delayLen   = delLen;
    firinterp->wrIx       = 0;

    return (firinterp);
} // firinterp16x16_init()


/* Put a chunk of input signal into the delay line and compute the filter
* response. */
void firinterp16x16_process(firinterp16x16_handle_t  _firinterp,
                             int16_t * restrict      y,
                       const int16_t * restrict      x, int N)
{
    firinterp16x16_ptr_t firinterp = (firinterp16x16_ptr_t)_firinterp;
    NASSERT(firinterp && firinterp->magic == MAGIC && y && x);
    if (N <= 0) return;

    NASSERT(N % 8 == 0);
    NASSERT(IS_ALIGN(x));
    //
    // Call filter's data processing function. It will store the block of input 
    // samples to the delay line, and compute the filter response. Returns the
    // updated next position pointer into the delay line buffer.
    //

    NASSERT(firinterp->procFxn);

    firinterp->wrIx = (*firinterp->procFxn)(
        y,
        firinterp->delayLine,
        firinterp->delayLen,
        x,
        firinterp->coef,
        firinterp->wrIx,
        firinterp->D,
        N,
        firinterp->M);
} // firinterp16x16_process()
#endif

