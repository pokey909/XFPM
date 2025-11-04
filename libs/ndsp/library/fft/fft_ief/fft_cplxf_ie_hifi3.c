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
	NatureDSP Signal Processing Library. FFT part
    FFT on Complex Data with Optimized Memory Usage
    C code optimized for HiFi3
	Integrit, 2006-2018
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fft.h"
#include "common_fpu.h"

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(int,fft_cplxf_ie,(complex_float * y, complex_float * x, const complex_float* twd, int twdstep, int N ))
#elif (HAVE_VFPU)
/*-------------------------------------------------------------------------
  FFT on Complex Data with Optimized Memory Usage
  These functions make FFT on complex data with optimized memory usage.
  Scaling: 
      +-------------------+----------------------------------------+
      |      Function     |           Scaling options              |
      +-------------------+----------------------------------------+
      |  fft_cplx16x16_ie |  2 - 16-bit dynamic scaling            | 
      |  fft_cplx24x24_ie |  3 - fixed scaling before each stage   | 
      |  fft_cplx32x16_ie |  3 - fixed scaling before each stage   | 
      |  fft_cplx32x32_ie |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      +-------------------+----------------------------------------+
  NOTES:
  1. Bit-reversing reordering is done here.
  2. FFT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED after 
     the call.
  3. FFT of size N may be supplied with constant data
     (twiddle factors) of a larger-sized FFT = N*twdstep.

  Precision: 
  16x16_ie      16-bit input/outputs, 16-bit twiddles
  24x24_ie      24-bit input/outputs, 24-bit twiddles
  32x16_ie      32-bit input/outputs, 16-bit twiddles
  32x32_ie      32-bit input/outputs, 32-bit twiddles
  f_ie          floating point
 
  Input:
  x[N]                  complex input signal. Real and imaginary data 
                        are interleaved and real data goes first
  twd[N*twdstep*3/4]    twiddle factor table of a complex-valued FFT of 
                        size N*twdstep
  N                     FFT size
  twdstep               twiddle step 
  scalingOpt            scaling option (see table above), not applicable
                        to the floating point function 
  Output:
  y[N]                  output spectrum. Real and imaginary data are 
                        interleaved and real data goes first

  Returned value: total number of right shifts occurred during scaling 
                  procedure. Floating point function always return 0.

  Restrictions:
  x,y   should not overlap
  x,y   aligned on 8-bytes boundary
-------------------------------------------------------------------------*/

#define SZ_CF32 (sizeof(complex_float))
#define LOG2_SZ_CF32 3/* log2(sizeof(complex_float)) */

int fft_cplxf_ie    (complex_float * y, complex_float * x, const complex_float* twd, int twdstep, int N )
{
  const xtfloatx2 *restrict p_twd;
  const xtfloatx2 *restrict p0_ld;
  const xtfloatx2 *restrict p1_ld;
  const xtfloatx2 *restrict p2_ld;
  const xtfloatx2 *restrict p3_ld;
        xtfloatx2 *restrict p0_st;
        xtfloatx2 *restrict p1_st;
        xtfloatx2 *restrict p2_st;
        xtfloatx2 *restrict p3_st;
  xtfloatx2 tw1, tw2, tw3;
  xtfloatx2 a0, a1, a2, a3;
  xtfloatx2 b0, b1, b2, b3;
  xtfloatx2 c1;
  int N4, logN, stride;
  int m, n;
  unsigned int idx, bitrevstride;

  NASSERT( x );
  NASSERT( y );
  NASSERT( twd );

  NASSERT( x != y );

  NASSERT_ALIGN( x, 8 );
  NASSERT_ALIGN( y, 8 );
  NASSERT_ALIGN( twd, 8 );

  NASSERT( twdstep >= 1 );
  NASSERT( N>=8 && 0 == (N&(N-1)) );

  N4 = N>>2;
  logN = 32 - NSA(N4);
  /* Helping constant, used for multiplication by j */
  c1 = (xtfloatx2)1.0f;

  /* Set the pointer to the twiddle table            *
   * and set bounds of the table for circular loads. */
  p_twd = (const xtfloatx2 *)(twd);
  WUR_AE_CBEGIN0((uintptr_t)(twd));
  WUR_AE_CEND0  ((uintptr_t)(twd)+AE_MOVAD32_L(AE_MULP32X2(3*twdstep,N))*2);
  //WUR_AE_CEND0  ((uintptr_t)(twd+3*twdstep*(N4)));
  /*----------------------------------------------------------------------------*
   * Perform the first stage. We use DIF, all permutations are deferred until   *
   * the last stage.                                                            */
  {
    stride = N4;
    
    p0_st = (xtfloatx2 *)(x);
    p1_st = (xtfloatx2 *)((uintptr_t)p0_st + 8 * stride);
    p2_st = (xtfloatx2 *)((uintptr_t)p1_st + 8 * stride);
    p3_st = (xtfloatx2 *)((uintptr_t)p2_st + 8 * stride);
    p0_ld = p0_st;
    p1_ld = p1_st;
    p2_ld = p2_st;
    p3_ld = p3_st;
    __Pragma("loop_count min=2")
    for ( n=0; n<stride; n++ )
    {
      tw2 = XT_LSX2I(p_twd,   SZ_CF32);
      XT_LSX2IP(tw1, p_twd, 2*SZ_CF32);
      XT_LSX2XC(tw3, p_twd, (twdstep*3-2)*SZ_CF32);

      XT_LSX2IP(a0, p0_ld, SZ_CF32);
      XT_LSX2IP(a1, p1_ld, SZ_CF32);
      XT_LSX2IP(a2, p2_ld, SZ_CF32);
      XT_LSX2IP(a3, p3_ld, SZ_CF32);

      b0 = a0 + a2;
      b1 = a1 + a3;
      b2 = a0 - a2;
      b3 = a1 - a3;

      a0 = b0 + b1;
      a2 = b0 - b1;
      /* a1 <- b2-j*b3 */
      /* a3 <- b2+j*b3 */
      a1 = a3 = b2;
      XT_MADDMUX_S(a1, c1, b3, 3);
      XT_MADDMUX_S(a3, c1, b3, 1);

      b0 = a0;
      b2 = XT_MULC_S(a2, tw2);
      b1 = XT_MULC_S(a1, tw1);
      b3 = XT_MULC_S(a3, tw3);

      /* Two middle quartiles are swapped on all but the last stage to use the bit reversal
       * permutation instead of the digit reverse. */
      XT_SSX2IP(b0, p0_st, SZ_CF32);
      XT_SSX2IP(b2, p1_st, SZ_CF32);
      XT_SSX2IP(b1, p2_st, SZ_CF32);
      XT_SSX2IP(b3, p3_st, SZ_CF32);
    }
  }
  /*----------------------------------------------------------------------------
   Perform second through the next to last stages.*/

  for ( stride>>=2; stride>1; stride>>=2 )
  {
    twdstep <<= 2;

    p0_st = (xtfloatx2 *)(x);

    for ( m=0; m<N4; m+=stride )
    {
      p1_st = (xtfloatx2 *)((uintptr_t)p0_st + 8 * stride);
      p2_st = (xtfloatx2 *)((uintptr_t)p1_st + 8 * stride);
      p3_st = (xtfloatx2 *)((uintptr_t)p2_st + 8 * stride);
      p0_ld = p0_st;
      p1_ld = p1_st;
      p2_ld = p2_st;
      p3_ld = p3_st;
      
      __Pragma("loop_count min=2")
      for ( n=0; n<stride; n++ )
      {
        tw2 = XT_LSX2I(p_twd,   SZ_CF32);
        XT_LSX2IP(tw1, p_twd, 2*SZ_CF32);
        XT_LSX2XC(tw3, p_twd, (twdstep*3-2)*SZ_CF32);
        
        XT_LSX2IP(a0, p0_ld, SZ_CF32);
        XT_LSX2IP(a1, p1_ld, SZ_CF32);
        XT_LSX2IP(a2, p2_ld, SZ_CF32);
        XT_LSX2IP(a3, p3_ld, SZ_CF32);

        b0 = a0 + a2;
        b1 = a1 + a3;
        b2 = a0 - a2;
        b3 = a1 - a3;
        
        a0 = b0 + b1;
        a2 = b0 - b1;
        /* a1 <- b2-j*b3 */
        /* a3 <- b2+j*b3 */
        a1 = a3 = b2;
        XT_MADDMUX_S(a1, c1, b3, 3);
        XT_MADDMUX_S(a3, c1, b3, 1);
        
        b0 = a0;
        b2 = XT_MULC_S(a2, tw2);
        b1 = XT_MULC_S(a1, tw1);
        b3 = XT_MULC_S(a3, tw3);

        /* Two middle quartiles are swapped on all but the last stage to use the bit reversal
         * permutation instead of the digit reverse. */
        XT_SSX2IP(b0, p0_st, SZ_CF32);
        XT_SSX2IP(b2, p1_st, SZ_CF32);
        XT_SSX2IP(b1, p2_st, SZ_CF32);
        XT_SSX2IP(b3, p3_st, SZ_CF32);
      }
      p0_st = p3_st;
    }
  }

  /*----------------------------------------------------------------------------
   Last stage (radix-4 or radix-2 for odd powers of two) with bit reversal
   permutation.*/
  idx = 0;
  bitrevstride = 0x80000000U >> (logN-3+LOG2_SZ_CF32);
  if ( stride != 0 )
  { 
    p0_ld = (const xtfloatx2 *)(x);
    p1_ld = p0_ld+1;
    p2_ld = p1_ld+1;
    p3_ld = p2_ld+1;
    p0_st = (xtfloatx2 *)(y);
    p1_st = p0_st+N4;
    p2_st = p1_st+N4;
    p3_st = p2_st+N4;
    
    __Pragma("loop_count min=2")
    for ( n=0; n<N4; n++ )
    {
      XT_LSX2IP(a0, p0_ld, 4*SZ_CF32);
      XT_LSX2IP(a1, p1_ld, 4*SZ_CF32);
      XT_LSX2IP(a2, p2_ld, 4*SZ_CF32);
      XT_LSX2IP(a3, p3_ld, 4*SZ_CF32);

      b0 = a0 + a2;
      b1 = a1 + a3;
      b2 = a0 - a2;
      b3 = a1 - a3;
      
      a0 = b0 + b1;
      a2 = b0 - b1;
      /* a1 <- b2-j*b3 */
      /* a3 <- b2+j*b3 */
      a1 = a3 = b2;
      XT_MADDMUX_S(a1, c1, b3, 3);
      XT_MADDMUX_S(a3, c1, b3, 1);
      
      XT_SSX2X(a0, p0_st, idx);
      XT_SSX2X(a1, p1_st, idx);
      XT_SSX2X(a2, p2_st, idx);
      XT_SSX2X(a3, p3_st, idx);

      idx = AE_ADDBRBA32(idx, bitrevstride);
    }
  }
  else
  {
    bitrevstride >>= 1;

    p0_ld = (const xtfloatx2 *)(x);
    p1_ld = p0_ld+1;
    p0_st = (xtfloatx2 *)(y);
    p1_st = (xtfloatx2 *)((uintptr_t)p0_st + 8 * (N4 << 1));

    __Pragma("loop_count min=1")
    for ( n=0; n<(N4<<1); n++ )
    {
      XT_LSX2IP(a0, p0_ld, 2*SZ_CF32);
      XT_LSX2IP(a1, p1_ld, 2*SZ_CF32);

      b0 = a0 + a1;
      b1 = a0 - a1;
      
      XT_SSX2X(b0, p0_st, idx);
      XT_SSX2X(b1, p1_st, idx);

      idx = AE_ADDBRBA32(idx, bitrevstride);
    }
  }
  return 0;
} /* fft_cplxf_ie() */

#else
// for scalar FPU
#define SZ_CF32 (sizeof(complex_float))
#define LOG2_SZ_CF32 3/* log2(sizeof(complex_float)) */

int fft_cplxf_ie    (complex_float * y, complex_float * x, const complex_float* twd, int twdstep, int N )
{
    int logN, stride;
    int m, n;
    unsigned int bitrevstride,idx;
    const xtfloat * restrict pTwd;
    xtfloat* restrict pY0;
    xtfloat* restrict pY1;
    xtfloat* restrict pY2;
    xtfloat* restrict pY3;
    xtfloat* restrict pX;

    NASSERT( x );
    NASSERT( y );
    NASSERT( twd );

    NASSERT( x != y );

    NASSERT_ALIGN( x, 8 );
    NASSERT_ALIGN( y, 8 );
    NASSERT_ALIGN( twd, 8 );

    NASSERT( twdstep >= 1 );
    NASSERT( N>=8 && 0 == (N&(N-1)) );

    logN = 30 - NSA(N);

    /*----------------------------------------------------------------------------*
    * Perform the first through the next to last stages. We use DIF,             *
    * all permutations are deferred until the last stage.                        */
    for ( stride=N/4; stride>1; stride>>=2 )
    {
        for ( m=0; m*(4*stride)<N; m++ )
        {
            pY0=(xtfloat*)(x+(4*m+0)*stride);
            pY1=(xtfloat*)(x+(4*m+1)*stride);
            pY2=(xtfloat*)(x+(4*m+2)*stride);
            pY3=(xtfloat*)(x+(4*m+3)*stride);
            pTwd=(const xtfloat *)twd;
            for ( n=0; n<stride; n++ )
            {
                xtfloat a0_re,a0_im;
                xtfloat a1_re,a1_im;
                xtfloat a2_re,a2_im;
                xtfloat a3_re,a3_im;
                xtfloat b0_re,b0_im;
                xtfloat b1_re,b1_im;
                xtfloat b2_re,b2_im;
                xtfloat b3_re,b3_im;
                xtfloat tw1_re,tw1_im;
                xtfloat tw2_re,tw2_im;
                xtfloat tw3_re,tw3_im;
                                                         tw1_im = XT_LSI(pTwd,1*sizeof(xtfloat)); 
                tw2_re = XT_LSI(pTwd,2*sizeof(xtfloat)); tw2_im = XT_LSI(pTwd,3*sizeof(xtfloat)); 
                tw3_re = XT_LSI(pTwd,4*sizeof(xtfloat)); tw3_im = XT_LSI(pTwd,5*sizeof(xtfloat)); 
                XT_LSXP(tw1_re,pTwd,twdstep*3*sizeof(complex_float)); 
                a0_re=XT_LSI(pY0,0); a0_im=XT_LSI(pY0,sizeof(xtfloat));
                a1_re=XT_LSI(pY1,0); a1_im=XT_LSI(pY1,sizeof(xtfloat));
                a2_re=XT_LSI(pY2,0); a2_im=XT_LSI(pY2,sizeof(xtfloat));
                a3_re=XT_LSI(pY3,0); a3_im=XT_LSI(pY3,sizeof(xtfloat));

                b0_re = XT_ADD_S(a0_re , a2_re); b0_im = XT_ADD_S(a0_im , a2_im);
                b1_re = XT_ADD_S(a1_re , a3_re); b1_im = XT_ADD_S(a1_im , a3_im);
                b2_re = XT_SUB_S(a0_re , a2_re); b2_im = XT_SUB_S(a0_im , a2_im);
                b3_re = XT_SUB_S(a1_re , a3_re); b3_im = XT_SUB_S(a1_im , a3_im);
 
                a0_re = XT_ADD_S(b0_re , b1_re); a0_im = XT_ADD_S(b0_im , b1_im);
                a1_re = XT_ADD_S(b2_re , b3_im); a1_im = XT_SUB_S(b2_im , b3_re);
                a2_re = XT_SUB_S(b0_re , b1_re); a2_im = XT_SUB_S(b0_im , b1_im);
                a3_re = XT_SUB_S(b2_re , b3_im); a3_im = XT_ADD_S(b2_im , b3_re);

                b0_re = a0_re;
                b0_im = a0_im;
                b1_re = XT_MUL_S(a1_re,tw1_re); XT_MSUB_S(b1_re,a1_im,tw1_im);
                b1_im = XT_MUL_S(a1_im,tw1_re); XT_MADD_S(b1_im,a1_re,tw1_im);
                b2_re = XT_MUL_S(a2_re,tw2_re); XT_MSUB_S(b2_re,a2_im,tw2_im);
                b2_im = XT_MUL_S(a2_im,tw2_re); XT_MADD_S(b2_im,a2_re,tw2_im);
                b3_re = XT_MUL_S(a3_re,tw3_re); XT_MSUB_S(b3_re,a3_im,tw3_im);
                b3_im = XT_MUL_S(a3_im,tw3_re); XT_MADD_S(b3_im,a3_re,tw3_im);

                /* Two middle quartiles are swapped on all but the last stage to use the bit reversal
                    * permutation instead of the digit reverse. */
                XT_SSIP(b0_re,pY0,sizeof(xtfloat)); XT_SSIP(b0_im,pY0,sizeof(xtfloat));
                XT_SSIP(b2_re,pY1,sizeof(xtfloat)); XT_SSIP(b2_im,pY1,sizeof(xtfloat));
                XT_SSIP(b1_re,pY2,sizeof(xtfloat)); XT_SSIP(b1_im,pY2,sizeof(xtfloat));
                XT_SSIP(b3_re,pY3,sizeof(xtfloat)); XT_SSIP(b3_im,pY3,sizeof(xtfloat));
            }
        }
        twdstep<<=2;
    }

    /*----------------------------------------------------------------------------
    Last stage (radix-8 or radix-4 for even powers of two) with bit reversal
    permutation.*/
    __Pragma("no_reorder");
    bitrevstride = 0x80000000U >> (logN-3+LOG2_SZ_CF32);
    idx=0;
    pX =(xtfloat*)(x);
    if ( stride == 0 )
    {
        /* radix-2 butterfly */
        bitrevstride >>= 1;
        pY0=(xtfloat*)(y+0*N/2);
        pY1=(xtfloat*)(y+1*N/2);
        __Pragma("loop_count min=1")
        for ( n=0; n<N/2; n++ )
        {
            xtfloat a0_re,a0_im;
            xtfloat a1_re,a1_im;
            xtfloat b0_re,b0_im;
            xtfloat b1_re,b1_im;
            XT_LSIP(a0_re,pX,sizeof(xtfloat)); XT_LSIP(a0_im,pX,sizeof(xtfloat));
            XT_LSIP(a1_re,pX,sizeof(xtfloat)); XT_LSIP(a1_im,pX,sizeof(xtfloat));

            b0_re = XT_ADD_S(a0_re , a1_re); b0_im = XT_ADD_S(a0_im , a1_im);
            b1_re = XT_SUB_S(a0_re , a1_re); b1_im = XT_SUB_S(a0_im , a1_im);
      
            XT_SSX(b0_re,pY0,idx); XT_SSX(b0_im,pY0,idx+sizeof(float32_t));
            XT_SSX(b1_re,pY1,idx); XT_SSX(b1_im,pY1,idx+sizeof(float32_t));
            idx = AE_ADDBRBA32(idx, bitrevstride);
        }
  }
  else
  {
    /* radix-4 butterfly */
    pY0=(xtfloat*)(y+0*N/4);
    pY1=(xtfloat*)(y+1*N/4);
    pY2=(xtfloat*)(y+2*N/4);
    pY3=(xtfloat*)(y+3*N/4);
    for ( n=0; n<N/4; n++ )
    {
        xtfloat a0_re,a0_im;
        xtfloat a1_re,a1_im;
        xtfloat a2_re,a2_im;
        xtfloat a3_re,a3_im;
        xtfloat b0_re,b0_im;
        xtfloat b1_re,b1_im;
        xtfloat b2_re,b2_im;
        xtfloat b3_re,b3_im;

        XT_LSIP(a0_re,pX,sizeof(xtfloat)); XT_LSIP(a0_im,pX,sizeof(xtfloat));
        XT_LSIP(a1_re,pX,sizeof(xtfloat)); XT_LSIP(a1_im,pX,sizeof(xtfloat));
        XT_LSIP(a2_re,pX,sizeof(xtfloat)); XT_LSIP(a2_im,pX,sizeof(xtfloat));
        XT_LSIP(a3_re,pX,sizeof(xtfloat)); XT_LSIP(a3_im,pX,sizeof(xtfloat));

        b0_re = XT_ADD_S(a0_re , a2_re); b0_im = XT_ADD_S(a0_im , a2_im);
        b1_re = XT_ADD_S(a1_re , a3_re); b1_im = XT_ADD_S(a1_im , a3_im);
        b2_re = XT_SUB_S(a0_re , a2_re); b2_im = XT_SUB_S(a0_im , a2_im);
        b3_re = XT_SUB_S(a1_re , a3_re); b3_im = XT_SUB_S(a1_im , a3_im);

        a0_re = XT_ADD_S(b0_re , b1_re); a0_im = XT_ADD_S(b0_im , b1_im);
        a1_re = XT_ADD_S(b2_re , b3_im); a1_im = XT_SUB_S(b2_im , b3_re);
        a2_re = XT_SUB_S(b0_re , b1_re); a2_im = XT_SUB_S(b0_im , b1_im);
        a3_re = XT_SUB_S(b2_re , b3_im); a3_im = XT_ADD_S(b2_im , b3_re);
        XT_SSX(a0_re,pY0,idx); XT_SSX(a0_im,pY0,idx+sizeof(float32_t));
        XT_SSX(a1_re,pY1,idx); XT_SSX(a1_im,pY1,idx+sizeof(float32_t));
        XT_SSX(a2_re,pY2,idx); XT_SSX(a2_im,pY2,idx+sizeof(float32_t));
        XT_SSX(a3_re,pY3,idx); XT_SSX(a3_im,pY3,idx+sizeof(float32_t));
        idx = AE_ADDBRBA32(idx, bitrevstride);
    }
  }
  return 0;
}

#endif
