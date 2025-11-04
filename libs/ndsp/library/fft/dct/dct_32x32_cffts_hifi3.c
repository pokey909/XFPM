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
    Discrete Cosine Transform, Type II
    C code optimized for HiFi3
   Integrit, 2006-2017
*/

/* Signal Processing Library API. */
#include "NatureDSP_Signal_fft.h"
#include "common.h"
#include "dct2_twd.h"
#include "fft_twiddles32x32.h"

static const int32_t ALIGN(8) __fft8_tw1[] =
{
  (int32_t)0x7FFFFFFF, (int32_t)0x00000000,
  (int32_t)0x2D413CCD, (int32_t)0xD2BEC333,
  (int32_t)0x00000000, (int32_t)0xC0000000,
  (int32_t)0xD2BEC333, (int32_t)0xD2BEC333,
};

#ifndef AE_MULFC32RAS
inline_ ae_int32x2 func_AE_MULFC32RAS(ae_int32x2 *a0, ae_int32x2 *a1)
{
    ae_f32x2 vA, vB;
    ae_f64 re, im;
    ae_f32x2 prod;

    vA = *a0; vB = *a1;
    re = AE_MULF32S_HH(vA, vB); AE_MULSF32S_LL(re, vA, vB);
    im = AE_MULF32S_LH(vA, vB); AE_MULAF32S_HL(im, vA, vB);

    prod = AE_ROUND32X2F64SASYM(re, im);
    return prod;
}
#define AE_MULFC32RAS(a,b) func_AE_MULFC32RAS(&(a),&(b))
#endif

#define DFT4X1(x0, x1, x2, x3)\
{   \
\
    ae_int32x2 s0, s1, d0, d1;                      \
    s0 = AE_ADD32S(x0, x2);                         \
    d0 = AE_SUB32S(x0, x2);                         \
    s1 = AE_ADD32S(x1, x3);                         \
    d1 = AE_SUB32S(x1, x3);                         \
    d1 = AE_SEL32_LH(d1, d1);                       \
    x0 = AE_ADD32S(s0, s1);                         \
    x2 = AE_SUB32S(s0, s1);                         \
    x1 = AE_ADDSUB32S(d0, d1);                      \
    x3 = AE_SUBADD32S(d0, d1);                      \
}

/*
   scaled fft with reordering
   N=16 - size of FFT
   NOTE: y is input and output, x - temporary
*/
void fft16_32x32(int32_t *y, int32_t *x, const int32_t *ptwd)
{
    #define N 16
    #define stride 4
    int k;
    ae_int32x2  * restrict p_x0,
                * restrict p_x1,
                * restrict p_x2,
                * restrict p_x3,
                * restrict p_y0;
    const ae_int32x2 * restrict p_twd;
    ae_int32x2  vA0, vA1, vA2, vA3, vB0, vB1, vB2, vB3;

    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);

    {
        p_x0  = (ae_int32x2 *)(x);
        p_x1  = (ae_int32x2 *)(x+2*N-2);
        p_y0  = (ae_int32x2 *)(y);

        /* reordering */
        __Pragma("loop_count min=2");
    	for (k=0; k<(N>>2); k++) 
        { 
            AE_L32X2_IP(vA0, p_y0, sizeof(ae_int32x2));
            AE_L32X2_IP(vA1, p_y0, sizeof(ae_int32x2));
            AE_L32X2_IP(vA2, p_y0, sizeof(ae_int32x2));
            AE_L32X2_IP(vA3, p_y0, sizeof(ae_int32x2));
            vB0 = AE_SEL32_HH(vA0, vA1);
            vB3 = AE_SEL32_LL(vA1, vA0);
            vB1 = AE_SEL32_HH(vA2, vA3);
            vB2 = AE_SEL32_LL(vA3, vA2);
            AE_S32X2_IP(vB0, p_x0,       sizeof(ae_int32x2));
            AE_S32X2_XP(vB3, p_x1, -(int)sizeof(ae_int32x2));
            AE_S32X2_IP(vB1, p_x0,       sizeof(ae_int32x2));
            AE_S32X2_XP(vB2, p_x1, -(int)sizeof(ae_int32x2));
        }
    }
#if 0
    fft_cplx32x32(y, x, cfft32_16, 3);/* fft of half-size with reordering */
#else
    __Pragma("no_reorder");
    p_x0 = (ae_int32x2 *)(x + 0*stride*2);
    p_x1 = (ae_int32x2 *)(x + 1*stride*2);
    p_x2 = (ae_int32x2 *)(x + 2*stride*2);
    p_x3 = (ae_int32x2 *)(x + 3*stride*2);
    p_y0 = (ae_int32x2 *)(y);

    p_twd = (const ae_int32x2 *)ptwd;
    __Pragma("loop_count min=3");
    /* hifi3z: 10 cycles per stage, unroll=1*/
    for (k = 0; k < stride; k++)
    {
      ae_int32x2 x0, x1, x2, x3;
      ae_int32x2 tw1, tw2, tw3;

      AE_L32X2_IP(x0, p_x0, sizeof(ae_int32x2));
      AE_L32X2_IP(x1, p_x1, sizeof(ae_int32x2));
      AE_L32X2_IP(x2, p_x2, sizeof(ae_int32x2));
      AE_L32X2_IP(x3, p_x3, sizeof(ae_int32x2));

      AE_L32X2_IP(tw1, p_twd, sizeof(ae_int32x2));
      AE_L32X2_IP(tw2, p_twd, sizeof(ae_int32x2));
      AE_L32X2_IP(tw3, p_twd, sizeof(ae_int32x2));

      x0 = AE_SRAI32R(x0, 3);
      x1 = AE_SRAI32R(x1, 3);
      x2 = AE_SRAI32R(x2, 3);
      x3 = AE_SRAI32R(x3, 3);

      DFT4X1(x0, x1, x2, x3);

      x1 = AE_MULFC32RAS(x1, tw1);
      x2 = AE_MULFC32RAS(x2, tw2);
      x3 = AE_MULFC32RAS(x3, tw3);


      AE_S32X2_IP(x0, p_y0, sizeof(ae_int32x2));
      AE_S32X2_IP(x1, p_y0, sizeof(ae_int32x2));
      AE_S32X2_IP(x2, p_y0, sizeof(ae_int32x2));
      AE_S32X2_IP(x3, p_y0, sizeof(ae_int32x2));
    }

    /* Last stage */
    __Pragma("no_reorder");
    {
      int j;
      p_x0 = (ae_int32x2 *)(0 * sizeof(ae_int32x2) + (uintptr_t)y);
      p_y0 = (ae_int32x2 *)(0 * sizeof(ae_int32x2) + (uintptr_t)y);
      __Pragma("loop_count min=3");
      for (j = 0; j < stride; j++)
      {
        ae_int32x2 x0, x1, x2, x3;

        x1 = AE_L32X2_I(p_x0, sizeof(ae_int32x2)*stride);
        x3 = AE_L32X2_X(p_x0, sizeof(ae_int32x2)*stride*3);
        x2 = AE_L32X2_X(p_x0, sizeof(ae_int32x2)*stride*2);
        AE_L32X2_IP(x0, p_x0, sizeof(ae_int32x2));

        x0 = AE_SRAI32R(x0, 2);
        x1 = AE_SRAI32R(x1, 2);
        x2 = AE_SRAI32R(x2, 2);
        x3 = AE_SRAI32R(x3, 2);

        DFT4X1(x0, x1, x2, x3);

        AE_S32X2_IP(x0, p_y0, stride*sizeof(ae_int32x2));
        AE_S32X2_IP(x1, p_y0, stride*sizeof(ae_int32x2));
        AE_S32X2_IP(x2, p_y0, stride*sizeof(ae_int32x2));
        AE_S32X2_XP(x3, p_y0, (-3*stride+1)*(int)sizeof(ae_int32x2));
      }
    }
#endif
    #undef stride
    #undef N
}




/*
   scaled fft with reordering
   N=32 - size of FFT
   NOTE: y is input and output, x - temporary
*/
void fft32_32x32(int32_t *y, int32_t *x, const int32_t *ptwd)
{
    #define N 32
    #define stride 8
    int k;
    ae_int32x2  * restrict p_x0,
                * restrict p_x1,
                * restrict p_x2,
                * restrict p_x3,
                * restrict p_y0;
    const ae_int32x2 * restrict p_twd;
    ae_int32x2  vA0, vA1, vA2, vA3, vB0, vB1, vB2, vB3;
    {
        p_x0  = (ae_int32x2 *)(x);
        p_x1  = (ae_int32x2 *)(x+2*N-2);
        p_y0  = (ae_int32x2 *)(y);

        /* reordering */
        __Pragma("loop_count min=2");
    	for (k=0; k<(N>>2); k++) 
        { 
            AE_L32X2_IP(vA0, p_y0, sizeof(ae_int32x2));
            AE_L32X2_IP(vA1, p_y0, sizeof(ae_int32x2));
            AE_L32X2_IP(vA2, p_y0, sizeof(ae_int32x2));
            AE_L32X2_IP(vA3, p_y0, sizeof(ae_int32x2));
            vB0 = AE_SEL32_HH(vA0, vA1);
            vB3 = AE_SEL32_LL(vA1, vA0);
            vB1 = AE_SEL32_HH(vA2, vA3);
            vB2 = AE_SEL32_LL(vA3, vA2);
            AE_S32X2_IP(vB0, p_x0,       sizeof(ae_int32x2));
            AE_S32X2_XP(vB3, p_x1, -(int)sizeof(ae_int32x2));
            AE_S32X2_IP(vB1, p_x0,       sizeof(ae_int32x2));
            AE_S32X2_XP(vB2, p_x1, -(int)sizeof(ae_int32x2));
        }
    }
#if 0
    fft_cplx32x32(y, x, cfft32_32, 3);/* fft of half-size with reordering */
#else
    __Pragma("no_reorder");
    p_x0 = (ae_int32x2 *)(x + 0*stride*2);
    p_x1 = (ae_int32x2 *)(x + 1*stride*2);
    p_x2 = (ae_int32x2 *)(x + 2*stride*2);
    p_x3 = (ae_int32x2 *)(x + 3*stride*2);
    p_y0 = (ae_int32x2 *)(y);

    p_twd = (const ae_int32x2 *)ptwd;
    __Pragma("loop_count min=3");
    /* hifi3z: 10 cycles per stage, unroll=1*/
    for (k = 0; k < stride; k++)
    {
      ae_int32x2 x0, x1, x2, x3;
      ae_int32x2 tw1, tw2, tw3;

      AE_L32X2_IP(x0, p_x0, sizeof(ae_int32x2));
      AE_L32X2_IP(x1, p_x1, sizeof(ae_int32x2));
      AE_L32X2_IP(x2, p_x2, sizeof(ae_int32x2));
      AE_L32X2_IP(x3, p_x3, sizeof(ae_int32x2));

      AE_L32X2_IP(tw1, p_twd, sizeof(ae_int32x2));
      AE_L32X2_IP(tw2, p_twd, sizeof(ae_int32x2));
      AE_L32X2_IP(tw3, p_twd, sizeof(ae_int32x2));

      x0 = AE_SRAI32R(x0, 3);
      x1 = AE_SRAI32R(x1, 3);
      x2 = AE_SRAI32R(x2, 3);
      x3 = AE_SRAI32R(x3, 3);

      DFT4X1(x0, x1, x2, x3);

      x1 = AE_MULFC32RAS(x1, tw1);
      x2 = AE_MULFC32RAS(x2, tw2);
      x3 = AE_MULFC32RAS(x3, tw3);


      AE_S32X2_IP(x0, p_y0, sizeof(ae_int32x2));
      AE_S32X2_IP(x1, p_y0, sizeof(ae_int32x2));
      AE_S32X2_IP(x2, p_y0, sizeof(ae_int32x2));
      AE_S32X2_IP(x3, p_y0, sizeof(ae_int32x2));
    }

    /* Last stage */
    __Pragma("no_reorder");
    {
      int j, stride_;
      stride_ = stride >> 1;
      p_x0 = (ae_int32x2 *)(7 * stride_ * sizeof(ae_int32x2) + (uintptr_t)y);
      p_y0 = (ae_int32x2 *)(0 * sizeof(ae_int32x2) + (uintptr_t)y);
      p_twd = (const ae_int32x2 *)__fft8_tw1;
      __Pragma("loop_count min=3");
      for (j = 0; j < stride_; j++)
      {
        ae_int32x2 tw1, tw2, tw3;
        ae_int32x2 x0, x1, x2, x3;
        ae_int32x2 x4, x5, x6, x7;

        tw1 = AE_L32X2_I(p_twd, 1 * sizeof(ae_int32x2));
        tw2 = AE_L32X2_I(p_twd, 2 * sizeof(ae_int32x2));
        tw3 = AE_L32X2_I(p_twd, 3 * sizeof(ae_int32x2));

        AE_L32X2_XP(x7, p_x0, -1*(sizeof(ae_int32x2)* stride_) );
        AE_L32X2_XP(x6, p_x0, -1*(sizeof(ae_int32x2)* stride_) );
        AE_L32X2_XP(x5, p_x0, -1*(sizeof(ae_int32x2)* stride_) );
        AE_L32X2_XP(x4, p_x0, -1*(sizeof(ae_int32x2)* stride_) );
        AE_L32X2_XP(x3, p_x0, -1*(sizeof(ae_int32x2)* stride_) );
        AE_L32X2_XP(x2, p_x0, -1*(sizeof(ae_int32x2)* stride_) );
        AE_L32X2_XP(x1, p_x0, -1*(sizeof(ae_int32x2)* stride_) );
        AE_L32X2_XP(x0, p_x0,     sizeof(ae_int32x2)*(stride_*7 + 1 ));

        x0 = AE_SRAI32R(x0, 2);
        x1 = AE_SRAI32R(x1, 2);
        x2 = AE_SRAI32R(x2, 2);
        x3 = AE_SRAI32R(x3, 2);
        x4 = AE_SRAI32R(x4, 2);
        x5 = AE_SRAI32R(x5, 2);
        x6 = AE_SRAI32R(x6, 2);
        x7 = AE_SRAI32R(x7, 2);

        DFT4X1(x0, x2, x4, x6); 
        DFT4X1(x1, x3, x5, x7);

        x3 = AE_MULFC32RAS(x3, tw1);
        x5 = AE_MULFC32RAS(x5, tw2);
        x7 = AE_MULFC32RAS(x7, tw3);

        {
          ae_int32x2 s0, s1, s2, s3;
          ae_int32x2 d0, d1, d2, d3;

          x0 = AE_SRAI32(x0, 1);
          x1 = AE_SRAI32(x1, 1);
          x2 = AE_SRAI32(x2, 1);
          x4 = AE_SRAI32(x4, 1);
          x6 = AE_SRAI32(x6, 1);

          s0 = AE_ADD32S(x0, x1);   d0 = AE_SUB32S(x0, x1);
          s1 = AE_ADD32S(x2, x3);   d1 = AE_SUB32S(x2, x3);
          s2 = AE_ADD32S(x4, x5);   d2 = AE_SUB32S(x4, x5);
          s3 = AE_ADD32S(x6, x7);   d3 = AE_SUB32S(x6, x7);

          x0 = s0;        x4 = d0;
          x1 = s1;        x5 = d1;
          x2 = s2;        x6 = d2;
          x3 = s3;        x7 = d3;
        }

        AE_S32X2_XP(x0, p_y0, stride_ * sizeof(ae_int32x2));
        AE_S32X2_XP(x1, p_y0, stride_ * sizeof(ae_int32x2));
        AE_S32X2_XP(x2, p_y0, stride_ * sizeof(ae_int32x2));
        AE_S32X2_XP(x3, p_y0, stride_ * sizeof(ae_int32x2));
        AE_S32X2_XP(x4, p_y0, stride_ * sizeof(ae_int32x2));
        AE_S32X2_XP(x5, p_y0, stride_ * sizeof(ae_int32x2));
        AE_S32X2_XP(x6, p_y0, stride_ * sizeof(ae_int32x2));
        AE_S32X2_XP(x7, p_y0, (-7 * stride_ + 1)* sizeof(ae_int32x2));
      }
    }
#endif
    #undef stride
    #undef N
}
