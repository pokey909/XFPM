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
  NatureDSP Signal Processing Library. IIR part
    Lattice Block Real IIR, 32x32-bit
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Lattice Block Real IIR
  Computes a real cascaded lattice autoregressive IIR filter using reflection 
  coefficients stored in vector k. The real data input are stored in vector x.
  The filter output result is stored in vector r. Input scaling is done before 
  the first cascade for normalization and overflow protection.

  Precision: 
  16x16   16-bit data, 16-bit coefficients
  24x24   24-bit data, 24-bit coefficients
  32x16   32-bit data, 16-bit coefficients
  32x32   32-bit data, 32-bit coefficients
  f       single precision floating point

  Input:
  N       length of input sample block
  M       number of reflection coefficients
  scale   input scale factor g, Q15, Q31 or floating point
  k[M]    reflection coefficients, Q15, Q31 or floating point
  x[N]    input samples, Q15, Q31 or floating point
  Output:
  r[N]    output data, Q15, Q31 or floating point

  Restriction:
  x,r,k   should not overlap

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  M - from the range 1...8
-------------------------------------------------------------------------*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_iir.h"
/* Common utility and macros declarations. */
#include "common.h"

/* Instance pointer validation number. */
#define MAGIC     0x1c6951b2

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
      ( (size_t)(size) + (align) - 1 )

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
      (void*)( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

#define sz_i32 sizeof(int32_t)

/* Puts the high 32-bit element to the low part of a register *
 * and extends the sign to high 32 bits                       */
#define AE_SEXT_64_32X2_H(inout)                    \
        {                                           \
            ae_f64 tmp = AE_MOVF64_FROMF32X2(inout);\
            tmp = AE_SRAI64(tmp, 32);               \
            inout = AE_MOVF32X2_FROMF64(tmp);       \
        }

/* Lattice filter data processing function. */
typedef void (proc_fxn_t)( int32_t * r,      // r[N]     [out   ] Q15
                     const int32_t * x,      // x[N]     [in    ] Q15
                           int32_t * delLine,// dline[M] [in/out] Q14
                     const int32_t * coef,   // coef[M]  [in    ] Q15
                           int32_t scale,    // scale    [in    ] Q15
                           int N, int M );

/* Custom data processing functions for particular lattice orders. */
static proc_fxn_t latr1_proc;
static proc_fxn_t latr2_proc;
static proc_fxn_t latr3_proc;
static proc_fxn_t latr4_proc;
static proc_fxn_t latr5_proc;
static proc_fxn_t latr6_proc;
static proc_fxn_t latr7_proc;
static proc_fxn_t latr8_proc;
static proc_fxn_t latrX_proc;
/* Custom data processing functions for particular lattice orders. */
static proc_fxn_t latr1_proc;
static proc_fxn_t latr2_proc;
static proc_fxn_t latr3_proc;
static proc_fxn_t latr4_proc;
static proc_fxn_t latr5_proc;
static proc_fxn_t latr6_proc;
static proc_fxn_t latr7_proc;
static proc_fxn_t latr8_proc;
static proc_fxn_t latrX_proc;
/* Table of processing functions */
static proc_fxn_t * const proc_fxn_tbl[9] =
{&latr1_proc, &latr2_proc, &latr3_proc, &latr4_proc, 
 &latr5_proc, &latr6_proc, &latr7_proc, &latr8_proc,
 &latrX_proc};


/* Filter instance structure. */
typedef struct tag_latr32x32_t
{
  uint32_t        magic;     // Instance pointer validation number
  int             M;         // Lattice filter order
  proc_fxn_t    * procFxn;   // Custom data processing function
  int32_t         scale;     // Input signal prescaling factor, Q31
  const int32_t * coef;      // M reflection coefficients, Q31
  int32_t       * delayLine; // M delay elements, Q30

} latr32x32_t, *latr32x32_ptr_t;

// Determine the memory area size for a filter instance.
size_t latr32x32_alloc( int M )
{
  ASSERT( M > 0 );

  return ( ALIGNED_SIZE( sizeof( latr32x32_t ), 4 )
           + // M delay elements
           ALIGNED_SIZE( M*sz_i32, 2*sz_i32 )
           + // M reflection coefficients
           ALIGNED_SIZE( M*sz_i32, 2*sz_i32 ) );

} /* latr32x32_alloc() */

// Given a memory area of sufficient size, initialize a filter instance and 
// return a handle to it.
latr32x32_handle_t latr32x32_init( void *             objmem, 
                                   int                M,
                             const int32_t * restrict k,
                                   int32_t            scale )
{
  latr32x32_ptr_t latr;
  int32_t *       delayLine;
  int32_t *       coef;
  void *          ptr;

  int m;

  ASSERT( objmem && M > 0 && k );

  //
  // Partition the memory block
  //

  ptr = objmem;

  latr      = (latr32x32_ptr_t)ALIGNED_ADDR( ptr, 4 );
  ptr       = latr + 1;
  delayLine = (int32_t*)ALIGNED_ADDR( ptr, 2*sz_i32 );
  ptr       = delayLine + M;
  coef      = (int32_t*)ALIGNED_ADDR( ptr, 2*sz_i32 );
  ptr       = coef + M;

  ASSERT( (int8_t*)ptr - (int8_t*)objmem <= (int)latr32x32_alloc( M ) );

  //
  // Copy reflection coefficients, zero the delay line.
  //

  for ( m=0; m<M; m++ )
  {
    coef[m] = k[m];

    delayLine[m] = 0;
  }

  //
  // Initialize the filter instance.
  //

  latr->magic     = MAGIC;
  latr->M         = M;
  latr->scale     = scale;
  latr->coef      = coef;
  latr->delayLine = delayLine;

  //
  // Set the correct processing function.
  //

  M = M>8 ? 8 : M-1;
  latr->procFxn = proc_fxn_tbl[M];

  return (latr);

} /* latr32x32_init() */

// Process data. The filter instance is identified with a handle returned by
// the initializing function.
void latr32x32_process( latr32x32_handle_t _latr, 
                        int32_t * restrict     r,
                  const int32_t *              x, int N )
{
  latr32x32_ptr_t latr = (latr32x32_ptr_t)_latr;

        int32_t * delLine;
  const int32_t * coef;
  int32_t scale;
  int M;

  ASSERT( latr && latr->magic == MAGIC && r && x );

  delLine = latr->delayLine;
  coef    = latr->coef;
  scale   = latr->scale;
  M       = latr->M;

  if (N <= 0) return;
  latr->procFxn(r, x, delLine, coef, scale, N, M);
} /* latr32x32_process() */
/* Order 1 lattice filter data processing function. */
static void latr1_proc( int32_t * restrict r,      // r[N]     [out   ] Q31
                  const int32_t * restrict x,      // x[N]     [in    ] Q31
                        int32_t *          delLine,// dline[M] [in/out] Q30
                  const int32_t *          coef,   // coef[M]  [in    ] Q31
                        int32_t            scale,  // scale    [in    ] Q31
                        int N, int M )
{
  const ae_f32x2 * restrict in;
        ae_f32x2 * restrict out;
  ae_f64  ACC;
  ae_int32x2 X;
  ae_f32x2 xin, rout;
  ae_f32x2 dl, sccf;
  int n;

  ASSERT(M == 1);
  // Set the input and output pointers
  in  = (const ae_f32x2 *)x;
  out = (      ae_f32x2 *)r;
  // Load the scale and reflection coefficients, the delay element
  dl = AE_L32_I((const ae_int32 *)delLine, 0);
  // sc : Q30 <- Q31 - 1 w/ rounding
  // cf : Q31
  sccf = AE_MOVDA32X2(((scale+1)>>1), *coef);
  rout = AE_ZERO32();

  __Pragma("loop_count min=1")
  for ( n=0; n<N; n++ )
  {
    // Load the input sample
    AE_L32_IP(X, castxcc(const ae_int32, in), sz_i32);
    xin = X;
    // Scale the input sample.
    // Q17.46 <- ( Q31*Q30 + 1 ) - 16 /w rounding
    ACC = AE_MULF32R_LH(xin, sccf);
    // Q17.46 <- Q17.46 - [( Q30*Q31 + 1 ) - 16]
    AE_MULSF32R_LL(ACC, dl, sccf);
    // Update the delay element with the resulting sample
    // Q30 <- Q17.46 - 16 w/ rounding and saturation
    dl = AE_ROUND32X2F48SASYM(ACC, ACC);
    // Make and store the output sample.
    // Q31 <- Q17.46 - 15 w/ rounding and saturation
    AE_PKSR32(rout, ACC, 1);
    AE_S32_L_IP(rout, castxcc( ae_int32, out), sz_i32);
  }
  // Save the delay line
  AE_S32_L_I(dl, (ae_int32 *)delLine, 0);
} /* latr1_proc() */

/* Order 2 lattice filter data processing function. */
static void latr2_proc( int32_t * restrict r,      // r[N]     [out   ] Q31
                  const int32_t * restrict x,      // x[N]     [in    ] Q31
                        int32_t *          delLine,// dline[M] [in/out] Q30
                  const int32_t *          coef,   // coef[M]  [in    ] Q31
                        int32_t            scale,  // scale    [in    ] Q31
                        int N, int M )
{
  const ae_f32x2   * restrict in;
        ae_f32x2   * restrict out;
        ae_int32   * restrict pdl;
  const ae_int32x2 * restrict pcoef;
  ae_f64  ACC, D1, D0;
  ae_int32x2 X;
  ae_f32x2 xin, r0, rout;
  ae_f32x2 dl1, dl0, cf01, sc;
  int n;

  ASSERT(M == 2);
  // Set the input and output pointers
  in    = (const ae_f32x2   *)x;
  out   = (      ae_f32x2   *)r;
  pdl   = (      ae_int32   *)delLine;
  pcoef = (const ae_int32x2 *)coef;
  // Load the scale and reflection coefficients, the delay elements
  dl0  = AE_L32_I(pdl  , 0);
  dl1  = AE_L32_I(pdl  , sz_i32);
  cf01 = AE_L32X2_I(pcoef, 0);
  sc = AE_MOVDA32X2(scale, scale);
  // D0 <- dl0
  // Q17.46 <- Q1.30
  D0 = AE_MOVF64_FROMF32X2(dl0);
  D0 = AE_SLAI64(D0, 32);
  D0 = AE_SRAI64(D0, 16);
  rout = AE_ZERO32();

  //
  // Pass the input samples block through the AR lattice. n-th response sample
  // and lattice state update are defined as follows:
  //
  //   r0 = sc*x[n] - cf1*dl1 - cf0*dl0
  //
  //   dl1 = dl0 + cf0*r0;
  //   dl0 = r0;
  //
  //   r[n] = r0;
  //
  // The inner loop is fully unrolled.
  //

  __Pragma("loop_count min=1")
  for ( n=0; n<N; n++ )
  {
    // Load the input sample
    AE_L32_IP(X, castxcc(const ae_int32, in), sz_i32);
    xin = X;

    // Scale the input sample.
    // Q17.46 <- ( Q31*Q31 + 1 ) - 17
    ACC = AE_MULF32R_LL(xin, sc);
    ACC = AE_SRAI64(ACC, 1);
    // Compute the output sample
    // Q17.46 <- Q17.46 - ( Q30*Q31 + 1 )
    AE_MULSF32R_LH(ACC, dl0, cf01);
    AE_MULSF32R_LL(ACC, dl1, cf01);
    // Q30 <- Q17.46 - 16 w/ rounding and saturation
    r0 = AE_ROUND32X2F48SASYM(ACC, ACC);
    // Update delay elements
    // dl1 = dl0 + r[n]*cf0
    // dl0 = r[n]
    D1 = D0;
    AE_MULAF32R_LH(D1, r0, cf01);
    // Q30 <- Q17.46 - 16 w/ rounding and saturation
    AE_PKSR32(dl1, D1, 0);
    D0 = ACC;
    dl0 = r0;

    // Make and store the output sample.
    // Q31 <- Q17.46 - 15 w/ rounding and saturation
    AE_PKSR32(rout, ACC, 1);
    AE_S32_L_IP(rout, castxcc(ae_int32, out), sz_i32);
  }
  // Save the updated delay line elements
  AE_S32_L_I(dl0, pdl, 0);
  AE_S32_L_I(dl1, pdl, sz_i32);

} /* latr2_proc() */

/* Order 3 lattice filter data processing function. */
static void latr3_proc( int32_t * restrict r,      // r[N]     [out   ] Q31
                  const int32_t * restrict x,      // x[N]     [in    ] Q31
                        int32_t *          delLine,// dline[M] [in/out] Q30
                  const int32_t *          coef,   // coef[M]  [in    ] Q31
                        int32_t            scale,  // scale    [in    ] Q31
                        int N, int M )
{
#if (XCHAL_HAVE_HIFI3Z)
    const ae_int32* restrict px=(const ae_int32*)x;
          ae_int32* restrict pr=(      ae_int32*)r;
    ae_int32x2 s,c01,c12,c23,xn,rn;
    ae_f32x2 d23,d01,d12,d_0,t0;
    ae_f64 q1,q0;
    int n;

    NASSERT(M==3);
    NASSERT_ALIGN8(coef);
    NASSERT_ALIGN8(delLine);
    rn=0;
    c01=AE_L32X2_I((const ae_int32x2*)coef,0);
    c23=AE_L32X2_I((const ae_int32x2*)coef,8);
    c12=AE_SEL32_LH(c01,c23);

    d01=AE_L32X2_I((const ae_int32x2*)delLine,0);
    d23=AE_L32X2_I((const ae_int32x2*)delLine,8);
    d12=AE_SEL32_LH(d01,d23);
    d_0=AE_SEL32_LH(  0,d01);
    s=AE_MOVDA32(scale>>1);
    __Pragma("loop_count min=1")
    for ( n=0; n<N; n++ )
    {
        AE_L32_IP(xn,px,4);
        q0=AE_MULF32S_LL(xn,s);
        AE_MULSSFD32S_HH_LL(q0,d12,c12);
        q1=q0;
        AE_MULSF32S_LH(q1,d_0,c01);
        t0=AE_ROUND32X2F64SASYM(q1,q0);
        d12=AE_SEL32_LH(d_0,d12);
        AE_MULAFP32X2RAS(d12,t0,c01);
        d_0=AE_ROUND32X2F64SASYM(q1,q1);
        rn=AE_TRUNCI32F64S_L(rn,q1,1);  // Make the output sample.
        AE_S32_L_IP(rn,pr,4);
    }
    d23=AE_SEL32_LH(d12,0);
    d01=AE_SEL32_LH(d_0,d12);
    AE_S32X2_I(d01,(ae_int32x2*)delLine,0);
    AE_S32X2_I(d23,(ae_int32x2*)delLine,8);
#else
  const ae_f32x2   * restrict in;
        ae_f32x2   * restrict out;
        ae_int32   * restrict pdl;
  const ae_int32x2 * restrict pcoef;
  ae_f64  ACC, R1, R0, D2, D1;
  ae_int32x2 X;
  ae_f32x2 xin, r1, r0, rout;
  ae_f32x2 dl2, dl1, dl0,
           cf01, cf2sc;
  int32_t cf2, ae_sar;
  int n;

  ASSERT(M == 3);

  // Set the input and output pointers
  in    = (const ae_f32x2   *)x;
  out   = (      ae_f32x2   *)r;
  pdl   = (      ae_int32   *)delLine;
  pcoef = (const ae_int32x2 *)coef;
  // Load the scale and reflection coefficients, the delay elements
  dl0 = AE_L32_I(pdl, 0*sz_i32);
  dl1 = AE_L32_I(pdl, 1*sz_i32);
  dl2 = AE_L32_I(pdl, 2*sz_i32);
  cf01 = AE_L32X2_I(pcoef, 0*sz_i32);
  cf2 = ((const int32_t *)pcoef)[2];
  cf2sc = AE_MOVDA32X2(cf2, scale);
  AE_SEXT_64_32X2_H(dl0);
  AE_SEXT_64_32X2_H(dl1);
  // Save and reset AE_SAR register
  ae_sar = RUR_AE_SAR();
  WUR_AE_SAR(0);
  
  //
  // Pass the input samples block through the AR lattice. n-th response sample
  // and lattice state update are defined as follows:
  //
  //   r1 = sc*x[n] - cf2*dl2 - cf1*dl1
  //   r0 = r1 - cf0*dl0;
  //
  //   dl2 = dl1 + cf1*r1;
  //   dl1 = dl0 + cf0*r0;
  //   dl0 = r0;
  //
  //   r[n] = r0;
  //
  // The inner loop is fully unrolled.
  //

  __Pragma("loop_count min=1")
  for ( n=0; n<N; n++ )
  {
    // Load the input sample
    AE_L32_IP(X, castxcc(const ae_int32, in), sz_i32);
    xin = X;
    // Scale the input sample.
    // Q17.46 <- ( Q31*Q31 + 1 ) - 17
    ACC = AE_MULF32R_LL(xin, cf2sc);
    ACC = AE_SRAI64(ACC, 1);
    // Compute output samples of each section
    // Q17.46 <- Q17.46 - ( Q30*Q31 + 1 )
    AE_MULSF32R_LH(ACC, dl2, cf2sc);
    AE_MULSF32R_LL(ACC, dl1, cf01);
    R1 = ACC;
    AE_MULSF32R_LH(ACC, dl0, cf01);
    R0 = ACC;
    // Q30 <- Q17.46 - 16 w/ rounding and saturation
    AE_PKSR32(r1, R1, 0);
    AE_PKSR32(r0, R0, 0);
    
    // Update delay elements
    // dl2 = dl1 + cf1*r1;
    // dl1 = dl0 + cf0*r0;
    // dl0 = r0;
    // D2 <- dl1; D1 <- dl0
    // Q17.46 <- Q30
    D2 = AE_SRA64_32(dl1, 0);
    D1 = AE_SRA64_32(dl0, 0);
    // Q17.46 <- Q17.46 + [(Q30*Q31 + 1) - 16] /w rounding
    AE_MULAF32R_LL(D2, r1, cf01);
    AE_MULAF32R_LH(D1, r0, cf01);
    // Q30 <- Q17.46 - 16 w/ rounding and saturation
    AE_PKSR32(dl2, D2, 0);
    AE_PKSR32(dl1, D1, 0);
    dl0 = r0;

    // Make and store the output sample.
    // Q31 <- Q17.46 - 15 w/ rounding and saturation
    AE_PKSR32(rout, ACC, 1);
    AE_S32_L_IP(rout, castxcc(ae_int32, out), sz_i32);
  }
  // Save the updated delay elements
  AE_S32_L_I(dl0, pdl, 0*sz_i32);
  AE_S32_L_I(dl1, pdl, 1*sz_i32);
  AE_S32_L_I(dl2, pdl, 2*sz_i32);
  // Restore AE_SAR
  WUR_AE_SAR(ae_sar);
#endif
} /* latr3_proc() */

/* Order 4 lattice filter data processing function. */
static void latr4_proc( int32_t * restrict r,      // r[N]     [out   ] Q31
                  const int32_t * restrict x,      // x[N]     [in    ] Q31
                        int32_t *          delLine,// dline[M] [in/out] Q30
                  const int32_t *          coef,   // coef[M]  [in    ] Q31
                        int32_t            scale,  // scale    [in    ] Q31
                        int N, int M )
{
    const ae_int32* restrict px=(const ae_int32*)x;
          ae_int32* restrict pr=(      ae_int32*)r;
    ae_int32x2 s,c01,c12,c23,c_0,xn,rn;
    ae_f32x2 d23,d01,t0;
    ae_f64 q0,q1;
    int n;

    NASSERT(M==4);
    NASSERT_ALIGN8(coef);
    NASSERT_ALIGN8(delLine);
    rn=0;
    c01=AE_L32X2_I((const ae_int32x2*)coef,0);
    c23=AE_L32X2_I((const ae_int32x2*)coef,8);
    c_0=AE_SEL32_HH(AE_MOVDA32(0x7fffffff),c01);
    c12=AE_SEL32_LH(c01,c23);
    d01=AE_L32X2_I((const ae_int32x2*)delLine,0);
    d23=AE_L32X2_I((const ae_int32x2*)delLine,8);
    s=AE_MOVDA32(scale>>1);

    __Pragma("loop_count min=1")
    __Pragma("ymemory(px)")
    for ( n=0; n<N; n++ )
    {
        AE_L32_IP(xn,px,4);
        q0=AE_MULF32S_LL(xn,s);
#if (XCHAL_HAVE_HIFI3Z)
        AE_MULSSFD32S_HH_LL(q0,d23,c23);
#else
        AE_MULSF32S_HH(q0,d23,c23);
        AE_MULSF32S_LL(q0,d23,c23);
#endif
        q1=q0;
        AE_MULSF32S_LL(q1,d01,c01);
        t0=AE_ROUND32X2F64SASYM(q1,q0);
        d23=AE_SEL32_LH(d01,d23);
        AE_MULAFP32X2RS(d23,t0,c12);

        AE_MULSF32S_HL(q1,d01,c_0);
        d01=AE_SEL32_HH(0,d01);
        t0=AE_ROUND32X2F64SASYM(q1,q1);
        AE_MULAFP32X2RS(d01,t0,c_0);

        rn=AE_TRUNCI32F64S_L(rn,q1,1);  // Make the output sample.
        AE_S32_L_IP(rn,pr,4);
    }
    AE_S32X2_I(d01,(ae_int32x2*)delLine,0);
    AE_S32X2_I(d23,(ae_int32x2*)delLine,8);
} /* latr4_proc() */

/* Order 5 lattice filter data processing function. */
static void latr5_proc( int32_t * restrict r,      // r[N]     [out   ] Q31
                  const int32_t * restrict x,      // x[N]     [in    ] Q31
                        int32_t *          delLine,// dline[M] [in/out] Q30
                  const int32_t *          coef,   // coef[M]  [in    ] Q31
                        int32_t            scale,  // scale    [in    ] Q31
                        int N, int M )
{
    const ae_int32* restrict px=(const ae_int32*)x;
          ae_int32* restrict pr=(      ae_int32*)r;
    ae_int32x2 s,c01,c12,c23,c45,c_0,c3,xn,rn;
    ae_f32x2 d45,d23,d01,t0;
    ae_f64 q0,q1;

    int  n;
    NASSERT(M==5);
    NASSERT_ALIGN8(coef);
    NASSERT_ALIGN8(delLine);
    rn=0;
    c01=AE_L32X2_I((const ae_int32x2*)coef,0);
    c23=AE_L32X2_I((const ae_int32x2*)coef,8);
    c45=AE_L32X2_I((const ae_int32x2*)coef,16);
    c_0=AE_SEL32_HH(AE_MOVDA32(0x7fffffff),c01);
    c12=AE_SEL32_LH(c01,c23);
    c3 =AE_SEL32_LL(c23,c23);
    d01=AE_L32X2_I((const ae_int32x2*)delLine,0);
    d23=AE_L32X2_I((const ae_int32x2*)delLine,8);
    d45=AE_L32X2_I((const ae_int32x2*)delLine,16);
    s=AE_MOVDA32(scale>>1);
    __Pragma("loop_count min=1")
    for ( n=0; n<N; n++ )
    {
        AE_L32_IP(xn,px,4);
        q0=AE_MULF32S_LL(xn,s);
        AE_MULSF32S_HH(q0,d45,c45);
        q1=q0;
        AE_MULSF32S_LL(q0,d23,c23);
        t0=AE_ROUND32X2F64SASYM(q0,q0);
        d45=AE_SEL32_LL(d23,d23);
        AE_MULAFP32X2RS(d45,t0,c3);

#if (XCHAL_HAVE_HIFI3Z)
        AE_MULSSFD32S_HH_LL(q1,d23,c23);
        q0=q1;
#else
        AE_MULSF32S_HH(q0,d23,c23);
        q1=q0;
#endif
        AE_MULSF32S_LL(q1,d01,c01);
        t0=AE_ROUND32X2F64SASYM(q1,q0);
        d23=AE_SEL32_LH(d01,d23);
        AE_MULAFP32X2RS(d23,t0,c12);

        AE_MULSF32S_HL(q1,d01,c_0);
        d01=AE_SEL32_HH(0,d01);
        t0=AE_ROUND32X2F64SASYM(q1,q1);
        AE_MULAFP32X2RS(d01,t0,c_0);

        rn=AE_TRUNCI32F64S_L(rn,q1,1);  // Make the output sample.
        AE_S32_L_IP(rn,pr,4);
    }
    AE_S32X2_I(d01,(ae_int32x2*)delLine,0);
    AE_S32X2_I(d23,(ae_int32x2*)delLine,8);
    AE_S32X2_I(d45,(ae_int32x2*)delLine,16);
} /* latr5_proc() */

/* Order 6 lattice filter data processing function. */
static void latr6_proc( int32_t * restrict r,      // r[N]     [out   ] Q31
                  const int32_t * restrict x,      // x[N]     [in    ] Q31
                        int32_t *          delLine,// dline[M] [in/out] Q30
                  const int32_t *          coef,   // coef[M]  [in    ] Q31
                        int32_t            scale,  // scale    [in    ] Q31
                        int N, int M )
{
    const ae_int32* restrict px=(const ae_int32*)x;
          ae_int32* restrict pr=(      ae_int32*)r;
    ae_int32x2 s,c01,c12,c23,c45,c34,c_0,xn,rn;
    ae_f32x2 d45,d23,d01,t0;
    ae_f64 q0,q1;

    int n;
    NASSERT(M==6);
    NASSERT_ALIGN8(coef);
    NASSERT_ALIGN8(delLine);
    rn=0;
    c01=AE_L32X2_I((const ae_int32x2*)coef,0);
    c23=AE_L32X2_I((const ae_int32x2*)coef,8);
    c45=AE_L32X2_I((const ae_int32x2*)coef,16);
    c_0=AE_SEL32_HH(AE_MOVDA32(0x7fffffff),c01);
    c34=AE_SEL32_LH(c23,c45);
    c12=AE_SEL32_LH(c01,c23);
    d01=AE_L32X2_I((const ae_int32x2*)delLine,0);
    d23=AE_L32X2_I((const ae_int32x2*)delLine,8);
    d45=AE_L32X2_I((const ae_int32x2*)delLine,16);
    s=AE_MOVDA32(scale>>1);

    __Pragma("loop_count min=1")
    for ( n=0; n<N; n++ )
    {
        AE_L32_IP(xn,px,4);
        q0=AE_MULF32S_LL(xn,s);
#if (XCHAL_HAVE_HIFI3Z)
        AE_MULSSFD32S_HH_LL(q0,d45,c45);
#else
        AE_MULSF32S_HH(q0,d45,c45);
        AE_MULSF32S_LL(q0,d45,c45);
#endif
        q1=q0;
        AE_MULSF32S_LL(q1,d23,c23);
        t0=AE_ROUND32X2F64SASYM(q1,q0);
        d45=AE_SEL32_LH(d23,d45);
        AE_MULAFP32X2RS(d45,t0,c34);

        AE_MULSF32S_HH(q1,d23,c23);
        q0=q1;
        AE_MULSF32S_LL(q0,d01,c01);
        t0=AE_ROUND32X2F64SASYM(q0,q1);
        d23=AE_SEL32_LH(d01,d23);
        AE_MULAFP32X2RS(d23,t0,c12);

        AE_MULSF32S_HL(q0,d01,c_0);
        d01=AE_SEL32_HH(0,d01);
        t0=AE_ROUND32X2F64SASYM(q0,q0);
        AE_MULAFP32X2RS(d01,t0,c_0);

        rn=AE_TRUNCI32F64S_L(rn,q0,1);  // Make the output sample.
        AE_S32_L_IP(rn,pr,4);
    }
    AE_S32X2_I(d01,(ae_int32x2*)delLine,0);
    AE_S32X2_I(d23,(ae_int32x2*)delLine,8);
    AE_S32X2_I(d45,(ae_int32x2*)delLine,16);
} /* latr6_proc() */

/* Order 7 lattice filter data processing function. */
static void latr7_proc( int32_t * restrict r,      // r[N]     [out   ] Q31
                  const int32_t * restrict x,      // x[N]     [in    ] Q31
                        int32_t *          delLine,// dline[M] [in/out] Q30
                  const int32_t *          coef,   // coef[M]  [in    ] Q31
                        int32_t            scale,  // scale    [in    ] Q31
                        int N, int M )
{
  const ae_f32x2 * restrict in;
        ae_f32x2 * restrict out;
        ae_f32x2 * restrict pdl;
  const ae_f32x2 * restrict pcoef;
  ae_f64  ACC, R5, R4, R3, R2, R1, R0;
  ae_int32x2 X;
  ae_f32x2 xin, r54, r32, r10, rout;
  ae_f32x2 dl65, dl43, dl21, dl0,
           sccf6, cf54, cf32, cf10;
  int32_t  cf6;
  int n;

  ASSERT(M == 7);
  
  // Set the input and output pointers
  in    = (const ae_f32x2 *)x;
  out   = (      ae_f32x2 *)r;
  pdl   = (      ae_f32x2 *)delLine;
  pcoef = (const ae_f32x2 *)coef;
  // Load the scale and reflection coefficients, the delay elements
  dl65 = ae_f32x2_loadi(pdl, 0*sz_i32);
  dl43 = ae_f32x2_loadi(pdl, 2*sz_i32);
  dl21 = ae_f32x2_loadi(pdl, 4*sz_i32);
  dl0  = AE_L32_I((ae_int32 *)pdl, 6*sz_i32);
  cf10 = ae_f32x2_loadi(pcoef, 0*sz_i32);
  cf32 = ae_f32x2_loadi(pcoef, 2*sz_i32);
  cf54 = ae_f32x2_loadi(pcoef, 4*sz_i32);
  cf10 = AE_SEL32_LH(cf10, cf10);
  cf32 = AE_SEL32_LH(cf32, cf32);
  cf54 = AE_SEL32_LH(cf54, cf54);
  cf6 = ((int32_t *)pcoef)[6];
  sccf6 = AE_MOVDA32X2(scale, cf6);
  
  //
  // Pass the input samples block through the AR lattice. n-th response sample
  // and lattice state update are defined as follows:
  //
  //   r5 = sc*x[n] - cf6*dl6 - cf5*dl5
  //   r4 = r5 - cf4*dl4;
  //   r3 = r4 - cf3*dl3;
  //   r2 = r3 - cf2*dl2;
  //   r1 = r2 - cf1*dl1;
  //   r0 = r1 - cf0*dl0;
  //
  //   dl6 = dl5 + cf5*r5;
  //   dl5 = dl4 + cf4*r4;
  //   dl4 = dl3 + cf3*r3;
  //   dl3 = dl2 + cf2*r2;
  //   dl2 = dl1 + cf1*r1;
  //   dl1 = dl0 + cf0*r0;
  //   dl0 = r0;
  //
  //   r[n] = r0;
  //
  // The inner loop is fully unrolled.
  //

  __Pragma("loop_count min=1")
  for ( n=0; n<N; n++ )
  {
    // Load the input sample
    AE_L32_IP(X, castxcc(const ae_int32, in), sz_i32);
    xin = X;
    // Scale the input sample.
    // Q17.46 <- ( Q31*Q31 + 1 ) - 17
    ACC = AE_MULF32R_LH(xin, sccf6);
    ACC = AE_SRAI64(ACC, 1);
    // Compute output samples of each section
    // Q17.46 <- Q17.46 - [( Q30*Q31 + 1 ) - 16]
    // Update delay elements
    // Q30 <- Q30 + [( Q30*Q31 + 1 ) - 32]

    // r5 = xin - cf6*dl6 - cf5*dl5;
    // r4 = r5  - cf4*dl4;
    AE_MULSF32R_HL(ACC, dl65, sccf6);
    AE_MULSF32R_LH(ACC, dl65, cf54);
    R5 = ACC;
    AE_MULSF32R_HL(ACC, dl43, cf54);
    R4 = ACC;
    r54 = AE_ROUND32X2F48SASYM(R5, R4);
    // dl6 = dl5 + cf5*r5;
    // dl5 = dl4 + cf4*r4;
    dl65 = AE_SEL32_LH(dl65, dl43);
    AE_MULAFP32X2RAS(dl65, r54, cf54);

    // r3 = r4 - cf3*dl3;
    // r2 = r3 - cf2*dl2;
    AE_MULSF32R_LH(ACC, dl43, cf32);
    R3 = ACC;
    AE_MULSF32R_HL(ACC, dl21, cf32);
    R2 = ACC;
    r32 = AE_ROUND32X2F48SASYM(R3, R2);
    // dl4 = dl3 + cf3*r3;
    // dl3 = dl2 + cf2*r2;
    dl43 = AE_SEL32_LH(dl43, dl21);
    AE_MULAFP32X2RAS(dl43, r32, cf32);

    // r1 = r2 - cf1*dl1;
    // r0 = r1 - cf0*dl0;
    AE_MULSF32R_LH(ACC, dl21, cf10);
    R1 = ACC;
    AE_MULSF32R_LL(ACC, dl0, cf10);
    R0 = ACC;
    r10 = AE_ROUND32X2F48SASYM(R1, R0);
    // dl2 = dl1 + cf1*r1;
    // dl1 = dl0 + cf0*r0;
    dl21 = AE_SEL32_LL(dl21, dl0 );
    AE_MULAFP32X2RAS(dl21, r10, cf10);
    // dl0 = r0;
    dl0 = r10;
    
    // Make and store the output sample.
    // Q31 <- Q17.46 - 15 w/ rounding and saturation
    AE_PKSR32(rout, ACC, 1);
    AE_S32_L_IP(rout, castxcc(ae_int32, out), sz_i32);
  }
  // Save the updated delay elements
  ae_f32x2_storei(dl65, pdl, 0*sz_i32);
  ae_f32x2_storei(dl43, pdl, 2*sz_i32);
  ae_f32x2_storei(dl21, pdl, 4*sz_i32);
  AE_S32_L_I(dl0, (ae_int32 *)pdl, 6*sz_i32);

} /* latr7_proc() */

/* Order 8 lattice filter data processing function. */
static void latr8_proc( int32_t * restrict r,      // r[N]     [out   ] Q31
                  const int32_t * restrict x,      // x[N]     [in    ] Q31
                        int32_t *          delLine,// dline[M] [in/out] Q30
                  const int32_t *          coef,   // coef[M]  [in    ] Q31
                        int32_t            scale,  // scale    [in    ] Q31
                        int N, int M )
{
  const ae_f32x2 * restrict in;
        ae_f32x2 * restrict out;
        ae_f32x2 * restrict pdl;
  const ae_f32x2 * restrict pcoef;
  ae_f64  ACC, R6, R5, R4, R3, R2, R1, R0;
  ae_int32x2 X;
  ae_f32x2 xin, r65, r43, r21, r0out;
  ae_f32x2 dl76, dl54, dl32, dl10,
           sccf7, cf65, cf43, cf21, cf0_;
  int n;

  ASSERT(M == 8);
  
  // Set the input and output pointers
  in    = (const ae_f32x2 *)x;
  out   = (      ae_f32x2 *)r;
  pdl   = (      ae_f32x2 *)delLine;
  pcoef = (const ae_f32x2 *)coef;
  // Load the scale and reflection coefficients, the delay elements
  dl76 = ae_f32x2_loadi(pdl, 0*sz_i32);
  dl54 = ae_f32x2_loadi(pdl, 2*sz_i32);
  dl32 = ae_f32x2_loadi(pdl, 4*sz_i32);
  dl10 = ae_f32x2_loadi(pdl, 6*sz_i32);
  cf0_ = ae_f32x2_loadi(pcoef, 0*sz_i32);// cf0, cf1
  cf21 = ae_f32x2_loadi(pcoef, 2*sz_i32);// cf2, cf3
  cf43 = ae_f32x2_loadi(pcoef, 4*sz_i32);// cf4, cf5
  cf65 = ae_f32x2_loadi(pcoef, 6*sz_i32);// cf6, cf7
  sccf7 = AE_MOVDA32X2(scale, scale);    // scale
  // Rearrange coefficients to the correct order
  sccf7 = AE_SEL32_HL(sccf7, cf65);// scale, cf7
  cf65  = AE_SEL32_HL(cf65, cf43); // cf6  , cf5
  cf43  = AE_SEL32_HL(cf43, cf21); // cf4  , cf3
  cf21  = AE_SEL32_HL(cf21, cf0_); // cf2  , cf1
  cf0_  = AE_SEL32_HL(cf0_, 0);    // cf0  , ___
  
  //
  // Pass the input samples block through the AR lattice. n-th response sample
  // and lattice state update are defined as follows:
  //
  //   r6 = sc*x[n] - cf7*dl7 - cf6*dl6
  //   r5 = r6 - cf5*dl5;
  //   r4 = r5 - cf4*dl4;
  //   r3 = r4 - cf3*dl3;
  //   r2 = r3 - cf2*dl2;
  //   r1 = r2 - cf1*dl1;
  //   r0 = r1 - cf0*dl0;
  //
  //   dl7 = dl6 + cf6*r6;
  //   dl6 = dl5 + cf5*r5;
  //   dl5 = dl4 + cf4*r4;
  //   dl4 = dl3 + cf3*r3;
  //   dl3 = dl2 + cf2*r2;
  //   dl2 = dl1 + cf1*r1;
  //   dl1 = dl0 + cf0*r0;
  //   dl0 = r0;
  //
  //   r[n] = r0;
  //
  // The inner loop is fully unrolled.
  //

  __Pragma("loop_count min=1")
  for ( n=0; n<N; n++ )
  {
    // Load the input sample
    AE_L32_IP(X, castxcc(const ae_int32, in), sz_i32);
    xin = X;
    // Scale the input sample.
    // Q17.46 <- ( Q31*Q31 + 1 ) - 17
    ACC = AE_MULF32R_LH(xin, sccf7);
    ACC = AE_SRAI64(ACC, 1);

    // Compute output samples of each section
    // Q17.46 <- Q17.46 - [( Q30*Q31 + 1 ) - 16]
    // Update delay elements
    // Q30 <- Q30 + [( Q30*Q31 + 1 ) - 32]

    // r6 = xin - cf7*dl7 - cf6*dl6;
    // r5 = r6  - cf5*dl5;
    AE_MULSF32R_HL(ACC, dl76, sccf7);
    AE_MULSF32R_LH(ACC, dl76, cf65);
    R6 = ACC;
    AE_MULSF32R_HL(ACC, dl54, cf65);
    R5 = ACC;
    r65   = AE_ROUND32X2F48SASYM(R6, R5);
    // dl7 = dl6 + cf6*r6;
    // dl6 = dl5 + cf5*r5;
    dl76 = AE_SEL32_LH(dl76, dl54);
    AE_MULAFP32X2RAS(dl76, r65, cf65);

    // r4 = r5 - cf4*dl4;
    // r3 = r4 - cf3*dl3;
    AE_MULSF32R_LH(ACC, dl54, cf43);
    R4 = ACC;
    AE_MULSF32R_HL(ACC, dl32, cf43);
    R3 = ACC;
    r43   = AE_ROUND32X2F48SASYM(R4, R3);
    // dl5 = dl4 + cf4*r4;
    // dl4 = dl3 + cf3*r3;
    dl54 = AE_SEL32_LH(dl54, dl32);
    AE_MULAFP32X2RAS(dl54, r43, cf43);

    // r2 = r3 - cf2*dl2;
    // r1 = r2 - cf1*dl1;
    AE_MULSF32R_LH(ACC, dl32, cf21);
    R2 = ACC;
    AE_MULSF32R_HL(ACC, dl10, cf21);
    R1 = ACC;
    r21   = AE_ROUND32X2F48SASYM(R2, R1);
    // dl3 = dl2 + cf2*r2;
    // dl2 = dl1 + cf1*r1;
    dl32 = AE_SEL32_LH(dl32, dl10);
    AE_MULAFP32X2RAS(dl32, r21, cf21);

    // r0 = r1 - cf0*dl0;
    AE_MULSF32R_LH(ACC, dl10, cf0_);
    R0 = ACC;
    // Make the output sample.
    // Q31 <- Q17.46 + 1 - 16 w/ rounding and saturation
    ACC   = AE_SLAI64S(ACC,1);
    r0out = AE_ROUND32X2F48SASYM(R0, ACC);
    // dl1 = dl0 + cf0*r0;
    // dl0 = r0;
    dl10 = AE_SEL32_LH(dl10, r0out);
    AE_MULAFP32X2RAS(dl10, r0out, cf0_);

    // Store the output sample.
    AE_S32_L_IP(r0out, castxcc(ae_int32, out), sz_i32);
  }
  // Save the updated delay elements
  ae_f32x2_storei(dl76, pdl, 0*sz_i32);
  ae_f32x2_storei(dl54, pdl, 2*sz_i32);
  ae_f32x2_storei(dl32, pdl, 4*sz_i32);
  ae_f32x2_storei(dl10, pdl, 6*sz_i32);

} /* latr8_proc() */

/* Data processing function for a lattice filter of arbitrary order. */
static void latrX_proc( int32_t * restrict r,      // r[N]     [out   ] Q31
                  const int32_t * restrict x,      // x[N]     [in    ] Q31
                        int32_t *          delLine,// dline[M] [in/out] Q30
                  const int32_t *          coef,   // coef[M]  [in    ] Q31
                        int32_t            scale,  // scale    [in    ] Q31
                        int N, int M )
{
  const ae_f32x2 * restrict in;
        ae_f32x2 * restrict out;
  const ae_f32x2 * restrict pdl_ld;
        ae_f32x2 * restrict pdl_st;
  const ae_f32x2 * restrict pcoef;
  ae_f64  ACC;
  ae_int32x2 X, DL, CF;
  ae_f32x2 xin, rm, rout;
  ae_f32x2 dl1, dl0, cf1, cf0, sc;
  int n, m;

  ASSERT(M > 8);

  sc = AE_MOVDA32X2(scale, scale);
  // Set the input and output pointers
  in  = (const ae_f32x2 *)x;
  out = (      ae_f32x2 *)r;

  __Pragma("loop_count min=1")
  for ( n=0; n<N; n++ )
  {
    // Set pointers to the delay elements and coefficients
    pdl_st = (      ae_f32x2 *)(delLine + M-1);
    pdl_ld = (const ae_f32x2 *)pdl_st;
    pcoef  = (      ae_f32x2 *)(coef + M-1);
    // Load the input sample x[n].
    AE_L32_IP(X, castxcc(const ae_int32, in), sz_i32);
    xin = X;
    // Scale the input sample.
    // Q17.46 <- ( Q31*Q31 + 1 ) - 17
    ACC = AE_MULF32R_LL(xin, sc);
    ACC = AE_SRAI64(ACC, 1);
    
    // Load (M-1)-th delay element and coefficient
    // dl[M-1] : Q30
    AE_L32_IP(DL, castxcc(const ae_int32, pdl_ld), -(int)sz_i32);
    dl1 = DL;
    // cf[M-1] : Q31
    AE_L32_IP(CF, castxcc(const ae_int32, pcoef), -(int)sz_i32);
    cf1 = CF;
    // acc = acc - dl[M-1]*cf[M-1]
    // Q17.46 <- Q17.46 - [( Q30*Q31 + 1 ) - 16] /w rounding
    AE_MULSF32R_LL(ACC, dl1, cf1);
    
    __Pragma("loop_count min=8")
    for ( m=M-2; m>=0; m-- )
    {
      // Load m-th delay element and coefficient
      // dl[m] : Q30
      AE_L32_IP(DL, castxcc(const ae_int32, pdl_ld), -(int)sz_i32);
      dl0 = DL;
      // cf[m] : Q31
      AE_L32_IP(CF, castxcc(const ae_int32, pcoef), -(int)sz_i32);
      cf0 = CF;
      // acc = acc - dl[m]*cf[m]
      // Q17.46 <- Q17.46 - [( Q30*Q31 + 1 ) - 16] /w rounding
      AE_MULSF32R_LL(ACC, dl0, cf0);
      // Get the output sample of m-th section
      // Q30 <- Q17.46 - 16 w/ rounding and saturation
      rm = AE_ROUND32X2F48SASYM(ACC, ACC);
      // Compute dl[m+1] delay element
      // dl[m+1] = dl[m] + r[m]*cf[m]
      // Q30 <- Q30 + [(Q30*Q31 + 1) - 32] /w rounding and saturation
      dl1 = dl0;
      AE_MULAFP32X2RAS(dl1, rm, cf0);
      // Update the (m+1)-th delay line element.
      AE_S32_L_IP(dl1, castxcc( ae_int32, pdl_st), -(int)sz_i32);
    }

    // Update the first delay line element with the resulting sample
    // Q30
    AE_S32_L_IP(rm, castxcc(ae_int32, pdl_st), -(int)sz_i32);
    
    // Make and store the output sample.
    // Q31 <- Q17.46 - 15 w/ rounding and saturation
    AE_PKSR32(rout, ACC, 1);
    AE_S32_L_IP(rout, castxcc(ae_int32, out), sz_i32);
  }

} /* latrX_proc() */
