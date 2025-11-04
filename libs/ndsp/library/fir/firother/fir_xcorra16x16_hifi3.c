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
    Real data circular cross-correlation, 16x16-bit, no requirements on vectors 
    length and alignment.
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Circular Correlation
  Estimates the circular cross-correlation between vectors x (of length N) 
  and y (of length M)  resulting in vector r of length N. It is a similar 
  to correlation but x is read in opposite direction.
  These functions implement the circular correlation algorithm with no 
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
              FIR_XCORRA16X16_SCRATCH_SIZE( N, M )
              FIR_XCORRA24X24_SCRATCH_SIZE( N, M ) or
              FIR_XCORRA32X16_SCRATCH_SIZE( N, M ) or
              FIR_XCORRA32X32_SCRATCH_SIZE( N, M ) or
              FIR_XCORRAF_SCRATCH_SIZE( N, M ) or
              CXFIR_XCORRAF_SCRATCH_SIZE( N, M ) bytes

  x[N]      input data Q15, Q31 or floating point
  y[M]      input data Q15, Q31 or floating point
  N         length of x
  M         length of y
  Output:
  r[N]      output data, Q15, Q31 or floating point

  Restrictions:
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

#include "raw_corr16x16.h"

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

#if (defined(AE_MULZAAAAQ16))
void fir_xcorra16x16( void    * restrict s,
                      int16_t * restrict r,
                const int16_t * restrict x,
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

    void     * s_ptr;
    int16_t  * x_buf;
    int16_t  * y_buf;
    const ae_int16x4 *          S;
    ae_int16x4       * restrict D;
    ae_int16         * restrict D_s;
    ae_valign S_va, D_va;
    ae_int16x4 y0, x0;
    ae_int16x4 z_f16x4;

    int n, m;

    NASSERT(s);
    NASSERT(r);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(s, 8);
    NASSERT(N>0 && M>0);
    NASSERT(N >= M - 1);

    //----------------------------------------------------------------------------
    // Partition the scratch memory area.
    s_ptr = s;
    x_buf = (int16_t*)ALIGNED_ADDR(s_ptr, 8);
    s_ptr = x_buf + N + M - 1 + 3;
    y_buf = (int16_t*)ALIGNED_ADDR(s_ptr, 8);
    s_ptr = y_buf + M + 3;

    ASSERT((int8_t *)s_ptr - (int8_t *)s <= (int)FIR_XCORRA16X16_SCRATCH_SIZE(N, M));
    //----------------------------------------------------------------------------
    // Copy x[N] data into the aligned buffer in a way that simplifies the
    // correlation calculation:
    //  x_buf[N+M-1+3] = { x[0]..x[N-1] x[0]..x[M-2] X X X }
    // Three X locations are reserved to allow for block-4 processing.
    // Copy x[N].
    S = (const ae_int16x4 *)x;
    D = (ae_int16x4 *)x_buf;
    S_va = AE_LA64_PP(S);
    for (n = 0; n<((N + 3) >> 2); n++)
    {
        AE_LA16X4_IP(x0, S_va, S);
        AE_S16X4_IP(x0, D, +8);
    }
    // Copy first M-1 entries of x{N].
    S = (const ae_int16x4 *)x;
    D = (ae_int16x4 *)((int16_t *)x_buf + N);
    S_va = AE_LA64_PP(S);
    D_va = AE_ZALIGN64();
    for (m = 0; m<((M + 3 - 1)>> 2); m++)
    {
        AE_LA16X4_IP(x0, S_va, S);
        AE_SA16X4_IP(x0, D_va, D);
    }
    AE_SA64POS_FP(D_va, D);
    // Append three zeros to allow for block-4 processing.
    D_s = (ae_int16 *)((int16_t *)x_buf + (N + M - 1));
    z_f16x4 = AE_ZERO16();
    for (m = 0; m<3; m++)
    {
        AE_S16_0_IP(z_f16x4, D_s, +2);
    }

    //----------------------------------------------------------------------------
    // Copy y[M] data into the aligned buffer and append 3 zeros:
    //  y_buf[M+3] = { y[0]..y[M-1] 0 0 0 }
    S = (const ae_int16x4 *)y;
    D = (ae_int16x4 *)y_buf;
    S_va = AE_LA64_PP(S);
    for (m = 0; m<((M + 3) >> 2); m++)
    {
        AE_LA16X4_IP(y0, S_va, S);
        AE_S16X4_IP(y0, D, +8);
    }
    // Append three zeros to allow for block-4 processing.
    D_s = (ae_int16 *)((int16_t *)y_buf + M);
    z_f16x4 = AE_ZERO16();
    for (m = 0; m<3; m++)
    {
        AE_S16_0_IP(z_f16x4, D_s, +2);
    }
    raw_corr16x16(r, x_buf, y_buf, N, (M + 3)&~3);
} /* fir_xcorra16x16() */

#else
void fir_xcorra16x16( void    * restrict s,
                      int16_t * restrict r,
                const int16_t * restrict x,
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

  void    * s_ptr;
  int16_t * x_buf;
  int16_t * y_buf;

  const ae_int16x4 *          S;
        ae_int16x4 * restrict D;
        ae_int16   * restrict D_s;
  const ae_f16x4   *          X;
  const ae_f16x4   *          SH;
  const ae_f16x4   *          Y;
        ae_f16     * restrict R;
        ae_f16x4   * restrict Rr;

  ae_valign S_va, D_va;

  ae_f64     q0, q1, q2, q3, q4, q5, q6, q7;
  ae_f16x4   x0, x1, x2, x4;
  ae_f32x2   d0, d1, d2, d3, d4, d5;
  ae_f16x4   y0;
  ae_int16x4 z_f16x4;
  ae_int32x2 temp0, temp1;
  ae_valign ar;

  int _M;
  int n, m;

  NASSERT(s);
  NASSERT(r);
  NASSERT(x);
  NASSERT(y);
  NASSERT_ALIGN(s,8);
  NASSERT((N>0) && (M>0));
  NASSERT(N>=(M-1));

  //----------------------------------------------------------------------------
  // Partition the scratch memory area.

  s_ptr = s;
  x_buf = (int16_t*)ALIGNED_ADDR( s_ptr, 8 );
  s_ptr = x_buf + N + M-1 + 3;
  y_buf = (int16_t*)ALIGNED_ADDR( s_ptr, 8 );
  s_ptr = y_buf + M + 3;

  ASSERT( (int8_t *)s_ptr - (int8_t *)s <= (int)FIR_XCORRA32X16_SCRATCH_SIZE( N, M ) );

  ar = AE_ZALIGN64();
  //----------------------------------------------------------------------------
  // Copy x[N] data into the aligned buffer in a way that simplifies the
  // correlation calculation:
  //  x_buf[N+M-1+3] = { x[0]..x[N-1] x[0]..x[M-2] X X X }
  // Three X locations are reserved to allow for block-4 processing.

  //
  // Copy x[N].
  //

  S = (const ae_int16x4 *)x;
  D = (ae_int16x4 *)x_buf;
  S_va = AE_LA64_PP(S);
  {
    ae_int16x4 x0, y0;
    for (n = 0; n<((N + 3) >> 2); n++)
    {
        AE_LA16X4_IP(x0, S_va, S);
        AE_S16X4_IP(x0, D, +8);
    }
    // Copy first M-1 entries of x{N].
    S = (const ae_int16x4 *)x;
    D = (ae_int16x4 *)((int16_t *)x_buf + N);
    S_va = AE_LA64_PP(S);
    D_va = AE_ZALIGN64();
    for (m = 0; m<((M + 3 - 1)>> 2); m++)
    {
        AE_LA16X4_IP(x0, S_va, S);
        AE_SA16X4_IP(x0, D_va, D);
    }
    AE_SA64POS_FP(D_va, D);
    // Append three zeros to allow for block-4 processing.
    D_s = (ae_int16 *)((int16_t *)x_buf + (N + M - 1));
    z_f16x4 = AE_ZERO16();
    for (m = 0; m<3; m++)
    {
        AE_S16_0_IP(z_f16x4, D_s, +2);
    }
    
    //----------------------------------------------------------------------------
    // Copy y[M] data into the aligned buffer and append 3 zeros:
    //  y_buf[M+3] = { y[0]..y[M-1] 0 0 0 }
    S = (const ae_int16x4 *)y;
    D = (ae_int16x4 *)y_buf;
    S_va = AE_LA64_PP(S);
    for (m = 0; m<((M + 3) >> 2); m++)
    {
        AE_LA16X4_IP(y0, S_va, S);
        AE_S16X4_IP (y0, D, +8);
    }
  }
  // Append three zeros to allow for block-4 processing.
  D_s = (ae_int16 *)((int16_t *)y_buf + M);
  z_f16x4 = AE_ZERO16();
  for (m = 0; m<3; m++)
  {
      AE_S16_0_IP(z_f16x4, D_s, +2);
  }

  //----------------------------------------------------------------------------
  // Compute (N&~7) correlation results.

  X = (const ae_f16x4 *)x_buf;
  Rr = (     ae_f16x4 *)r;

  _M = M + ( -M & 3 );

  //
  // Process vector x data in 8-entries blocks.
  //

  for ( n=0; n<(N>>3); n++ )
  {
    // Load 8 left-most unprocessed x[] entries, the first is x[8*n].
    // Q31
    ae_f16x4_loadip( x0, X, +8 );
    ae_f16x4_loadip( x1, X, +8 );

    d0 = AE_CVT32X2F16_32(x0);
    d1 = AE_CVT32X2F16_10(x0);
    d2 = AE_CVT32X2F16_32(x1);
    d3 = AE_CVT32X2F16_10(x1);

    // Use the shuttle pointer when computing the correlation. Preserve the X
    // pointer for the next iteration.
    SH = X;

    //
    // Inner loop prologue: process first 4 y[] entries for 8 accumulators.
    //

    // Load 4 more x[] entries. Altogether we have 12 x[] entries stored in
    // 6 AE registers.
    // Q31
    ae_f16x4_loadip( x4, SH, +8 );

    d4 = AE_CVT32X2F16_32(x4);
    d5 = AE_CVT32X2F16_10(x4);

    Y = (const ae_f16x4 *)y_buf;

    // Load y[0]..y[3].
    // Q15
    ae_f16x4_loadip( y0, Y, +8 );

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH( q0, q1, d0, d1, y0 );
    AE_MULFD32X16X2_FIR_HH( q2, q3, d1, d2, y0 );
    AE_MULFD32X16X2_FIR_HH( q4, q5, d2, d3, y0 );
    AE_MULFD32X16X2_FIR_HH( q6, q7, d3, d4, y0 );

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, y0 );
    AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, y0 );
    AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, y0 );
    AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, y0 );

    // 4 x[] entries are done, shift them out of the registers.
    d0 = d2; d1 = d3;
    d2 = d4; d3 = d5;

    //
    // Inner loop kernel: process 4 y[] entries for 8 accumulators. 12-entries 
    // register delay line is updated similarly to the loop prologue.
    //

    for ( m=0; m<(_M>>2)-1; m++ )
    {
      // Q31
      ae_f16x4_loadip( x4, SH, +8 );

      d4 = AE_CVT32X2F16_32(x4);
      d5 = AE_CVT32X2F16_10(x4);

      // Load y[4*(m+1)]..y[4*(m+1)+3].
      // Q15
      ae_f16x4_loadip( y0, Y, +8 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH( q0, q1, d0, d1, y0 );
      AE_MULAFD32X16X2_FIR_HH( q2, q3, d1, d2, y0 );
      AE_MULAFD32X16X2_FIR_HH( q4, q5, d2, d3, y0 );
      AE_MULAFD32X16X2_FIR_HH( q6, q7, d3, d4, y0 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, y0 );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, y0 );
      AE_MULAFD32X16X2_FIR_HL( q4, q5, d3, d4, y0 );
      AE_MULAFD32X16X2_FIR_HL( q6, q7, d4, d5, y0 );

      // 4 x[] entries are done, shift them out of the registers.
      d0 = d2; d1 = d3;
      d2 = d4; d3 = d5;
    }

    //
    // Convert and save 8 correlation results.
    //

    // Q31 <- Q16.47 - 16 w/ rounding and saturation.
    temp0 = AE_TRUNCA32X2F64S(q0, q1, 16);
    temp1 = AE_TRUNCA32X2F64S(q2, q3, 16);
    AE_SA16X4_IP(AE_ROUND16X4F32SASYM(temp0, temp1), ar, Rr); 
    temp0 = AE_TRUNCA32X2F64S(q4, q5, 16);               
    temp1 = AE_TRUNCA32X2F64S(q6, q7, 16);   
    AE_SA16X4_IP(AE_ROUND16X4F32SASYM(temp0, temp1), ar, Rr); 
  }

  //----------------------------------------------------------------------------
  // (N&~7) correlation results are done by now. If N - (N&~7) >= 4, then
  // compute a block of 4 results.

  if ( N & 4 )
  {
    // Load 4 left-most unprocessed x[] entries.
    // Q15
    ae_f16x4_loadip( x0, X, +8 );
    // Q31<-Q15
    d0 = AE_CVT32X2F16_32(x0);
    d1 = AE_CVT32X2F16_10(x0);

    // Use the shuttle pointer when computing the correlation. Preserve the X
    // pointer for the last block of (N&3) results.
    SH = X;

    //
    // Inner loop prologue: process first 4 y[] entries for 4 accumulators.
    //

    // Load 4 more x[] entries. Altogether we have 8 x[] entries residing in
    // 4 AE registers.
    // Q15
    ae_f16x4_loadip( x2, SH, +8 );
    // Q31<-Q15
    d2 = AE_CVT32X2F16_32(x2);
    d3 = AE_CVT32X2F16_10(x2);

    Y = (const ae_f16x4 *)y_buf;

    // Load y[0]..y[3].
    // Q15
    ae_f16x4_loadip( y0, Y, +8 );

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH( q0, q1, d0, d1, y0 );
    AE_MULFD32X16X2_FIR_HH( q2, q3, d1, d2, y0 );

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, y0 );
    AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, y0 );

    // 4 x[] entries are done, shift them out of the registers.
    d0 = d2; d1 = d3;

    //
    // Inner loop kernel: process 4 y[] entries for 4 accumulators. 8-entries 
    // register delay line is updated similarly to the loop prologue.
    //

    for ( m=0; m<(_M>>2)-1; m++ )
    {
      // Q15
      ae_f16x4_loadip( x2, SH, +8 );
      // Q31<-Q15
      d2 = AE_CVT32X2F16_32(x2);
      d3 = AE_CVT32X2F16_10(x2);

      // Load y[4*(m+1)]..y[4*(m+1)+3].
      // Q15
      ae_f16x4_loadip( y0, Y, +8 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH( q0, q1, d0, d1, y0 );
      AE_MULAFD32X16X2_FIR_HH( q2, q3, d1, d2, y0 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, y0 );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, y0 );

      // 4 x[] entries are done, shift them out of the registers.
      d0 = d2; d1 = d3;
    }

    //
    // Convert and save 4 correlation results.
    //

    // Q31 <- Q16.47 - 16 w/ rounding and saturation.
    temp0 = AE_TRUNCA32X2F64S(q0, q1, 16);
    temp1 = AE_TRUNCA32X2F64S(q2, q3, 16);
    AE_SA16X4_IP(AE_ROUND16X4F32SASYM(temp0, temp1), ar, Rr); 
  }

  //----------------------------------------------------------------------------
  // (N&~3) correlation results are done by now. If N is not a multiple of 4,
  // compute the last (N&3) results.

  if ( N & 3 )
  {
    // Load 4 left-most unprocessed x[] entries.
    // Q15
    ae_f16x4_loadip( x0, X, +8 );
    // Q31<-Q15
    d0 = AE_CVT32X2F16_32(x0);
    d1 = AE_CVT32X2F16_10(x0);

    //
    // Inner loop prologue: process first 4 y[] entries for 4 accumulators.
    //

    // Load 4 more x[] entries. Altogether we have 8 x[] entries residing in
    // 4 AE registers.
    // Q15
    ae_f16x4_loadip( x2, X, +8 );
    // Q31<-Q15
    d2 = AE_CVT32X2F16_32(x2);
    d3 = AE_CVT32X2F16_10(x2);

    Y = (const ae_f16x4 *)y_buf;

    // Load y[0]..y[3].
    // Q15
    ae_f16x4_loadip( y0, Y, +8 );

    // Q16.47 <- [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULFD32X16X2_FIR_HH( q0, q1, d0, d1, y0 );
    AE_MULFD32X16X2_FIR_HH( q2, q3, d1, d2, y0 );

    // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
    AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, y0 );
    AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, y0 );

    // 4 x[] entries are done, shift them out of the registers.
    d0 = d2; 
    d1 = d3;

    //
    // Inner loop kernel: process 4 y[] entries for 4 accumulators. 8-entries 
    // register delay line is updated similarly to the loop prologue.
    //

    for ( m=0; m<(_M>>2)-1; m++ )
    {
      // Q15
      ae_f16x4_loadip( x2, X, +8 );
      // Q31<-Q15
      d2 = AE_CVT32X2F16_32(x2);
      d3 = AE_CVT32X2F16_10(x2);

      // Load y[4*(m+1)]..y[4*(m+1)+3].
      // Q15
      ae_f16x4_loadip( y0, Y, +8 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HH( q0, q1, d0, d1, y0 );
      AE_MULAFD32X16X2_FIR_HH( q2, q3, d1, d2, y0 );

      // Q16.47 <- Q16.47 + [ Q31*Q15 + 1 ] + [ Q31*Q15 + 1 ]
      AE_MULAFD32X16X2_FIR_HL( q0, q1, d1, d2, y0 );
      AE_MULAFD32X16X2_FIR_HL( q2, q3, d2, d3, y0 );

      // 4 x[] entries are done, shift them out of the registers.
      d0 = d2;
      d1 = d3;
    }

    //
    // Convert and save (N&3) correlation results.
    //

    R = (ae_f16 *)( r + N-1 );

    switch ( N & 3 )
    {
      // Q31 <- Q16.47 - 16 w/ rounding and saturation.
      case 3: 
      {
        temp0 = AE_TRUNCA32X2F64S(0, q2, 16);
        AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp0), R, -2);
      }
      case 2: 
      {
        temp0 = AE_TRUNCA32X2F64S(0, q1, 16);
        AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp0), R, -2);
      }
    }

    // Q31 <- Q16.47 - 16 w/ rounding and saturation.
    temp0 = AE_TRUNCA32X2F64S(0, q0, 16);
    AE_S16_0_XP(AE_ROUND16X4F32SASYM(0, temp0), R, -2);
  }
  AE_SA64POS_FP(ar, Rr);
} /* fir_xcorra16x16() */

#endif
