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

#define FIRST_STAGE_SCALE 3
   
#define FFT_BUTTERFLY_R2(idx)              \
{                                          \
    vA1l = AE_L32X2_I(p_x0, 8);            \
    AE_L32X2_IP(vA0l,p_x0, 16);            \
    vB0l = AE_SEL32_HH(vA0l, vA1l);        \
    vB1l = AE_SEL32_LL(vA0l, vA1l);        \
    vA2s = AE_MOVINT16X4_FROMF32X2(vB0l);  \
    vA3s = AE_MOVINT16X4_FROMF32X2(vB1l);  \
    vA0s = AE_ADD16S(vA2s, vA3s);          \
    vA1s = AE_SUB16S(vA2s, vA3s);          \
    vA0l = AE_MOVINT32X2_FROMINT16X4(vA0s);\
    vA1l = AE_MOVINT32X2_FROMINT16X4(vA1s);\
    AE_S32_L_X(vA0l, p_y1, (idx));         \
    vA0l = AE_SEL32_LH(vA0l, vA0l);        \
    AE_S32_L_X(vA1l, p_y3, (idx));         \
    AE_S32_L_X(vA0l, p_y0, (idx));         \
    vA1l = AE_SEL32_LH(vA1l, vA1l);        \
    AE_S32_L_X(vA1l, p_y2, (idx));         \
}

#define FFT_BUTTERFLY_R4(idx)              \
{                                          \
    vA1l = AE_L32X2_I(p_x0, 8);            \
    AE_L32X2_IP(vA0l,p_x0, 16);            \
    vA2s = AE_MOVINT16X4_FROMF32X2(vA0l);  \
    vA3s = AE_MOVINT16X4_FROMF32X2(vA1l);  \
    vA0s = AE_ADD16S(vA2s, vA3s);          \
    vA1s = AE_SUB16S(vA2s, vA3s);          \
    vA2l=AE_MULC32X16_L(vA3l,vA1s);        \
    vB1s=AE_CVT16X4(vA2l,vA2l);            \
    vA2l=AE_MOVINT32X2_FROMINT16X4(vB1s);  \
    vA0l = AE_MOVINT32X2_FROMINT16X4(vA0s);\
    vA1l = AE_MOVINT32X2_FROMINT16X4(vA1s);\
    vA1l = AE_SEL32_HH(vA0l,vA1l);         \
    vA0l = AE_SEL32_LL(vA0l,vA2l);         \
    vA2s = AE_MOVINT16X4_FROMF32X2(vA1l);  \
    vA3s = AE_MOVINT16X4_FROMF32X2(vA0l);  \
    vA0s = AE_ADD16S(vA2s, vA3s);          \
    vA1s = AE_SUB16S(vA2s, vA3s);          \
    vA0l = AE_MOVINT32X2_FROMINT16X4(vA0s);\
    vA1l = AE_MOVINT32X2_FROMINT16X4(vA1s);\
    AE_S32_L_X(vA0l, p_y3, (idx));         \
    AE_S32_L_X(vA1l, p_y1, (idx));         \
    vA0l = AE_SEL32_LH(vA0l, vA0l);        \
    AE_S32_L_X(vA0l, p_y0, (idx));         \
    vA1l = AE_SEL32_LH(vA1l, vA1l);        \
    AE_S32_L_X(vA1l, p_y2, (idx));         \
}

static const int16_t ALIGN(8) fft_twd_r8[12] = {
    (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000,
    (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x5a82,(int16_t)0xA57E,
    (int16_t)0x0000,(int16_t)0x8000,(int16_t)0xa57e,(int16_t)0xA57E
};

/*
   scaled fft with reordering
   N=16 - size of FFT
   NOTE: y is input and output, x - temporary
*/
void fft16_16x16(int16_t *y, int16_t *x, const int16_t *ptwd)
{
    const int N = 16;
    int k;
    ae_int16x4 * px0;
    ae_int64   * px1;
    ae_int16x4 * restrict py0;
    ae_int16x4 * restrict p4_x0,
               * restrict p4_x1,
               * restrict p4_x2,
               * restrict p4_x3;
    const ae_int16x4 * restrict p_twd;

    ae_int32 * restrict p_y0;
    ae_int32 * restrict p_y1;
    ae_int32 * restrict p_y2;
    ae_int32 * restrict p_y3;
    ae_int32x2 * restrict p_x0;
    ae_int16x4  vA0s, vA1s, vA2s, vA3s, vB1s, vB2s, vB3s,
                vC1s, vC2s, vC3s, vT0s, vT1s, vT2s;
    ae_int32x2  vA0l, vA1l, vA2l, vA3l, vB0l, vB1l, vB2l, vB3l, vC0l, vC1l, vC2l, vC3l;
    ae_int16x4  vA0, vA1, vB0, vB1;
    ae_int64    vB_64;

    {
        px0  = (ae_int16x4 *)(x);
        px1  = (ae_int64   *)(x+2*N-4);
        py0  = (ae_int16x4 *)(y);

        /* reordering */
    	for (k=0; k<(N>>2); k++) 
        { 
            AE_L16X4_IP(vA0, py0, sizeof(ae_int16x4));
            AE_L16X4_IP(vA1, py0, sizeof(ae_int16x4));
            vB0 = AE_SEL16I    (vA0, vA1, 7);/* even numbers */
            vB1 = AE_SEL16_6420(vA0, vA1);
            vB_64 = AE_MOVINT64_FROMINT16X4(vB1);
            AE_S16X4_IP(vB0  , px0,       sizeof(ae_int16x4));
            AE_S64_IP  (vB_64, px1, -(int)sizeof(ae_int16x4));
        }
    }
    /* fft of half-size */
    __Pragma("no_reorder");
    p_twd = (const ae_int16x4 *)ptwd;

    { // first iteration for pre-scaling
      p4_x0 = (ae_int16x4 *)x;
      p4_x1 = p4_x0 + (N>>3);
      p4_x2 = p4_x1 + (N>>3);
      p4_x3 = p4_x2 + (N>>3);

      //-----------------------------------------------------------------
      // Set up offsets to access "N/4", "N/2", "3N/4" complex point or  
      // "N/2", "N", "3N/2" half word                                    
      //-----------------------------------------------------------------

      vA0s = AE_L16X4_I(p4_x0, 0);
      vA1s = AE_L16X4_I(p4_x1, 0);
      vA2s = AE_L16X4_I(p4_x2, 0);
      vA3s = AE_L16X4_I(p4_x3, 0);

      vA0s = AE_SRAI16(vA0s, FIRST_STAGE_SCALE);
      vA1s = AE_SRAI16(vA1s, FIRST_STAGE_SCALE);
      vA2s = AE_SRAI16(vA2s, FIRST_STAGE_SCALE);
      vA3s = AE_SRAI16(vA3s, FIRST_STAGE_SCALE);

      vB2s = AE_SUB16S(vA0s, vA2s);
      vB3s = AE_SUB16S(vA1s, vA3s);
      vA0s = AE_ADD16S(vA0s, vA2s);
      vA1s = AE_ADD16S(vA1s, vA3s);

      vT2s = AE_L16X4_I(p_twd, 16);
      vB1s = AE_SUB16S(vA0s, vA1s);
      vC3s = AE_ADD16S(vA0s, vA1s);
      vT1s = AE_L16X4_I(p_twd, 8);
      AE_L16X4_IP(vT0s, p_twd, 24);
      vC0l = AE_CVT32X2F16_32(vB1s);
      vC1l = AE_CVT32X2F16_10(vB1s);
      vC3s = AE_SRAI16(vC3s, 2);
      vC0l = AE_SRAI32(vC0l, 2);
      AE_S16X4_IP(vC3s, p4_x0, 8);
      vC0l = AE_MULFC32X16RAS_L(vC0l, vT0s);
      vC1l = AE_SRAI32(vC1l, 2);
      vA0s = AE_L16X4_I(p4_x0, 0);
      vB1l = AE_CVT32X2F16_32(vB3s);
      vC1l = AE_MULFC32X16RAS_H(vC1l, vT2s);
      vA2s = AE_L16X4_I(p4_x2, 8);
      vB0l = AE_CVT32X2F16_32(vB2s);
      vB1l = AE_SEL32_LH(vB1l, vB1l);
      vC1s = AE_ROUND16X4F32SASYM(vC0l, vC1l);
      vB3l = AE_CVT32X2F16_10(vB3s);
      vC0l = AE_SUBADD32S(vB0l, vB1l);
      vB3l = AE_SEL32_LH(vB3l, vB3l);
      vB2l = AE_CVT32X2F16_10(vB2s);
      vC1l = AE_ADDSUB32S(vB0l, vB1l);
      vA0s = AE_SRAI16(vA0s, FIRST_STAGE_SCALE);
      vC2l = AE_SUBADD32S(vB2l, vB3l);
      vA2s = AE_SRAI16(vA2s, FIRST_STAGE_SCALE);
      vC3l = AE_ADDSUB32S(vB2l, vB3l);
      vC0l = AE_SRAI32(vC0l, 2);
      vC0l = AE_MULFC32X16RAS_H(vC0l, vT1s);
      vC2l = AE_SRAI32(vC2l, 2);
      vA1s = AE_L16X4_I(p4_x1, 8);
      vC2l = AE_MULFC32X16RAS_L(vC2l, vT2s);
      vC1l = AE_SRAI32(vC1l, 2);
      vC1l = AE_MULFC32X16RAS_H(vC1l, vT0s);
      vC3l = AE_SRAI32(vC3l, 2);
      vC3l = AE_MULFC32X16RAS_L(vC3l, vT1s);
      vA1s = AE_SRAI16(vA1s, FIRST_STAGE_SCALE);
      vA3s = AE_L16X4_I(p4_x3, 8);
      vC2s = AE_ROUND16X4F32SASYM(vC0l, vC2l);
      vB2s = AE_SUB16S(vA0s, vA2s);
      vC3s = AE_ROUND16X4F32SASYM(vC1l, vC3l);
      AE_S16X4_IP(vC2s, p4_x3, 8);
      vA0s = AE_ADD16S(vA0s, vA2s);
      vA3s = AE_SRAI16(vA3s, FIRST_STAGE_SCALE);
      AE_S16X4_IP(vC1s, p4_x1, 8);
      AE_S16X4_IP(vC3s, p4_x2, 8);
      vB3s = AE_SUB16S(vA1s, vA3s);
      vA1s = AE_ADD16S(vA1s, vA3s);

      // last iteration
      vT2s = AE_L16X4_I(p_twd, 16);
      vB1s = AE_SUB16S(vA0s, vA1s);
      vC3s = AE_ADD16S(vA0s, vA1s);
      vT1s = AE_L16X4_I(p_twd, 8);
      AE_L16X4_IP(vT0s, p_twd, 24);
      vC0l = AE_CVT32X2F16_32(vB1s);
      vC1l = AE_CVT32X2F16_10(vB1s);
      vC3s = AE_SRAI16(vC3s, 2);
      vC0l = AE_SRAI32(vC0l, 2);
      AE_S16X4_IP(vC3s, p4_x0, 8);
      vC0l = AE_MULFC32X16RAS_L(vC0l, vT0s);
      vC1l = AE_SRAI32(vC1l, 2);
      vB1l = AE_CVT32X2F16_32(vB3s);
      vC1l = AE_MULFC32X16RAS_H(vC1l, vT2s);
      vB0l = AE_CVT32X2F16_32(vB2s);
      vB1l = AE_SEL32_LH(vB1l, vB1l);
      vC1s = AE_ROUND16X4F32SASYM(vC0l, vC1l);
      vB3l = AE_CVT32X2F16_10(vB3s);
      vC0l = AE_SUBADD32S(vB0l, vB1l);
      vB3l = AE_SEL32_LH(vB3l, vB3l);
      vB2l = AE_CVT32X2F16_10(vB2s);
      vC1l = AE_ADDSUB32S(vB0l, vB1l);
      vA0s = AE_SRAI16(vA0s, FIRST_STAGE_SCALE);
      vC2l = AE_SUBADD32S(vB2l, vB3l);
      vA2s = AE_SRAI16(vA2s, FIRST_STAGE_SCALE);
      vC3l = AE_ADDSUB32S(vB2l, vB3l);
      vC0l = AE_SRAI32(vC0l, 2);
      vC0l = AE_MULFC32X16RAS_H(vC0l, vT1s);
      vC2l = AE_SRAI32(vC2l, 2);
      vC2l = AE_MULFC32X16RAS_L(vC2l, vT2s);
      vC1l = AE_SRAI32(vC1l, 2);
      vC1l = AE_MULFC32X16RAS_H(vC1l, vT0s);
      vC3l = AE_SRAI32(vC3l, 2);
      vC3l = AE_MULFC32X16RAS_L(vC3l, vT1s);
      vA1s = AE_SRAI16(vA1s, FIRST_STAGE_SCALE);
      vC2s = AE_ROUND16X4F32SASYM(vC0l, vC2l);
      vC3s = AE_ROUND16X4F32SASYM(vC1l, vC3l);
      AE_S16X4_IP(vC2s, p4_x3, 8);
      vA3s = AE_SRAI16(vA3s, FIRST_STAGE_SCALE);
      AE_S16X4_IP(vC1s, p4_x1, 8);
      AE_S16X4_IP(vC3s, p4_x2, 8);
    }

    __Pragma("no_reorder");
    {
      p_y0=(ae_int32 *)(y);
      p_y1=(ae_int32 *)(p_y0 + (N >> 2));
      p_y2=(ae_int32 *)(p_y1 + (N >> 2));
      p_y3=(ae_int32 *)(p_y2 + (N >> 2));         
      p_x0=(ae_int32x2 *)(x);
      int32_t i0,i1,i2,i3;

      vA3l=AE_MOVDA32X2(0,-1);
      //--------------------------------------------------------------------------
      // last stage is RADIX4 !!!
      //--------------------------------------------------------------------------
      i0 = 0;
      i1 = 8;
      i2 = 4;
      i3 = 12;
      FFT_BUTTERFLY_R4(i0);
      FFT_BUTTERFLY_R4(i1);
      FFT_BUTTERFLY_R4(i2);
      FFT_BUTTERFLY_R4(i3);
    }
}




/*
   scaled fft with reordering
   N=32 - size of FFT
   NOTE: y is input and output, x - temporary
*/
void fft32_16x16(int16_t *y, int16_t *x, const int16_t *ptwd)
{
    const int N = 32;
    int k;
    ae_int16x4 * px0;
    ae_int64   * px1;
    ae_int16x4 * restrict py0;
    ae_int16x4 * restrict p4_x0,
               * restrict p4_x1,
               * restrict p4_x2,
               * restrict p4_x3;
    const ae_int16x4 * restrict p_twd;

    ae_int32 * restrict p_y0;
    ae_int32 * restrict p_y1;
    ae_int32 * restrict p_y2;
    ae_int32 * restrict p_y3;
    ae_int32x2 * restrict p_x0;
    ae_int16x4  vA0s, vA1s, vA2s, vA3s, vB1s, vB2s, vB3s,
                vC1s, vC2s, vC3s, vT0s, vT1s, vT2s;
    ae_int32x2  vA0l, vA1l, vB0l, vB1l, vB2l, vB3l, vC0l, vC1l, vC2l, vC3l;

    int stride;
    ae_int16x4  vA0, vA1, vB0, vB1;
    ae_int64    vB_64;

    {
        px0  = (ae_int16x4 *)(x);
        px1  = (ae_int64   *)(x+2*N-4);
        py0  = (ae_int16x4 *)(y);

        /* reordering */
        __Pragma("loop_count min=2");
    	for (k=0; k<(N>>2); k++) 
        { 
            AE_L16X4_IP(vA0, py0, sizeof(ae_int16x4));
            AE_L16X4_IP(vA1, py0, sizeof(ae_int16x4));
            vB0 = AE_SEL16I    (vA0, vA1, 7);/* even numbers */
            vB1 = AE_SEL16_6420(vA0, vA1);
            vB_64 = AE_MOVINT64_FROMINT16X4(vB1);
            AE_S16X4_IP(vB0  , px0,       sizeof(ae_int16x4));
            AE_S64_IP  (vB_64, px1, -(int)sizeof(ae_int16x4));
        }
    }
    /* fft of half-size */
    __Pragma("no_reorder");
    p_twd = (const ae_int16x4 *)ptwd;
    stride = N; // The stride is quartered with every iteration of the outer loop.

    { // first iteration for pre-scaling
      int lc = (N>>3)-1;

      p4_x0 = (ae_int16x4 *)x;
      p4_x1 = p4_x0 + (stride>>3);
      p4_x2 = p4_x1 + (stride>>3);
      p4_x3 = p4_x2 + (stride>>3);

      //-----------------------------------------------------------------
      // Set up offsets to access "N/4", "N/2", "3N/4" complex point or  
      // "N/2", "N", "3N/2" half word                                    
      //-----------------------------------------------------------------

      vA0s = AE_L16X4_I(p4_x0, 0);
      vA1s = AE_L16X4_I(p4_x1, 0);
      vA2s = AE_L16X4_I(p4_x2, 0);
      vA3s = AE_L16X4_I(p4_x3, 0);

      vA0s = AE_SRAI16(vA0s, FIRST_STAGE_SCALE);
      vA1s = AE_SRAI16(vA1s, FIRST_STAGE_SCALE);
      vA2s = AE_SRAI16(vA2s, FIRST_STAGE_SCALE);
      vA3s = AE_SRAI16(vA3s, FIRST_STAGE_SCALE);

      vB2s = AE_SUB16S(vA0s, vA2s);
      vB3s = AE_SUB16S(vA1s, vA3s);
      vA0s = AE_ADD16S(vA0s, vA2s);
      vA1s = AE_ADD16S(vA1s, vA3s);

      do
      {
        vT2s = AE_L16X4_I(p_twd, 16);
        vB1s = AE_SUB16S(vA0s, vA1s);
        vC3s = AE_ADD16S(vA0s, vA1s);
        vT1s = AE_L16X4_I(p_twd, 8);
        AE_L16X4_IP(vT0s, p_twd, 24);
        vC0l = AE_CVT32X2F16_32(vB1s);
        vC1l = AE_CVT32X2F16_10(vB1s);
        vC3s = AE_SRAI16(vC3s, 2);
        vC0l = AE_SRAI32(vC0l, 2);
        AE_S16X4_IP(vC3s, p4_x0, 8);
        vC0l = AE_MULFC32X16RAS_L(vC0l, vT0s);
        vC1l = AE_SRAI32(vC1l, 2);
        vA0s = AE_L16X4_I(p4_x0, 0);
        vB1l = AE_CVT32X2F16_32(vB3s);
        vC1l = AE_MULFC32X16RAS_H(vC1l, vT2s);
        vA2s = AE_L16X4_I(p4_x2, 8);
        vB0l = AE_CVT32X2F16_32(vB2s);
        vB1l = AE_SEL32_LH(vB1l, vB1l);
        vC1s = AE_ROUND16X4F32SASYM(vC0l, vC1l);
        vB3l = AE_CVT32X2F16_10(vB3s);
        vC0l = AE_SUBADD32S(vB0l, vB1l);
        vB3l = AE_SEL32_LH(vB3l, vB3l);
        vB2l = AE_CVT32X2F16_10(vB2s);
        vC1l = AE_ADDSUB32S(vB0l, vB1l);
        vA0s = AE_SRAI16(vA0s, FIRST_STAGE_SCALE);
        vC2l = AE_SUBADD32S(vB2l, vB3l);
        vA2s = AE_SRAI16(vA2s, FIRST_STAGE_SCALE);
        vC3l = AE_ADDSUB32S(vB2l, vB3l);
        vC0l = AE_SRAI32(vC0l, 2);
        vC0l = AE_MULFC32X16RAS_H(vC0l, vT1s);
        vC2l = AE_SRAI32(vC2l, 2);
        vA1s = AE_L16X4_I(p4_x1, 8);
        vC2l = AE_MULFC32X16RAS_L(vC2l, vT2s);
        vC1l = AE_SRAI32(vC1l, 2);
        vC1l = AE_MULFC32X16RAS_H(vC1l, vT0s);
        vC3l = AE_SRAI32(vC3l, 2);
        vC3l = AE_MULFC32X16RAS_L(vC3l, vT1s);
        vA1s = AE_SRAI16(vA1s, FIRST_STAGE_SCALE);
        vA3s = AE_L16X4_I(p4_x3, 8);
        vC2s = AE_ROUND16X4F32SASYM(vC0l, vC2l);
        vB2s = AE_SUB16S(vA0s, vA2s);
        vC3s = AE_ROUND16X4F32SASYM(vC1l, vC3l);
        AE_S16X4_IP(vC2s, p4_x3, 8);
        vA0s = AE_ADD16S(vA0s, vA2s);
        vA3s = AE_SRAI16(vA3s, FIRST_STAGE_SCALE);
        AE_S16X4_IP(vC1s, p4_x1, 8);
        AE_S16X4_IP(vC3s, p4_x2, 8);
        vB3s = AE_SUB16S(vA1s, vA3s);
        vA1s = AE_ADD16S(vA1s, vA3s);
      } while(--lc);
      // last iteration
      vT2s = AE_L16X4_I(p_twd, 16);
      vB1s = AE_SUB16S(vA0s, vA1s);
      vC3s = AE_ADD16S(vA0s, vA1s);
      vT1s = AE_L16X4_I(p_twd, 8);
      AE_L16X4_IP(vT0s, p_twd, 24);
      vC0l = AE_CVT32X2F16_32(vB1s);
      vC1l = AE_CVT32X2F16_10(vB1s);
      vC3s = AE_SRAI16(vC3s, 2);
      vC0l = AE_SRAI32(vC0l, 2);
      AE_S16X4_IP(vC3s, p4_x0, 8);
      vC0l = AE_MULFC32X16RAS_L(vC0l, vT0s);
      vC1l = AE_SRAI32(vC1l, 2);
      vB1l = AE_CVT32X2F16_32(vB3s);
      vC1l = AE_MULFC32X16RAS_H(vC1l, vT2s);
      vB0l = AE_CVT32X2F16_32(vB2s);
      vB1l = AE_SEL32_LH(vB1l, vB1l);
      vC1s = AE_ROUND16X4F32SASYM(vC0l, vC1l);
      vB3l = AE_CVT32X2F16_10(vB3s);
      vC0l = AE_SUBADD32S(vB0l, vB1l);
      vB3l = AE_SEL32_LH(vB3l, vB3l);
      vB2l = AE_CVT32X2F16_10(vB2s);
      vC1l = AE_ADDSUB32S(vB0l, vB1l);
      vA0s = AE_SRAI16(vA0s, FIRST_STAGE_SCALE);
      vC2l = AE_SUBADD32S(vB2l, vB3l);
      vA2s = AE_SRAI16(vA2s, FIRST_STAGE_SCALE);
      vC3l = AE_ADDSUB32S(vB2l, vB3l);
      vC0l = AE_SRAI32(vC0l, 2);
      vC0l = AE_MULFC32X16RAS_H(vC0l, vT1s);
      vC2l = AE_SRAI32(vC2l, 2);
      vC2l = AE_MULFC32X16RAS_L(vC2l, vT2s);
      vC1l = AE_SRAI32(vC1l, 2);
      vC1l = AE_MULFC32X16RAS_H(vC1l, vT0s);
      vC3l = AE_SRAI32(vC3l, 2);
      vC3l = AE_MULFC32X16RAS_L(vC3l, vT1s);
      vA1s = AE_SRAI16(vA1s, FIRST_STAGE_SCALE);
      vC2s = AE_ROUND16X4F32SASYM(vC0l, vC2l);
      vC3s = AE_ROUND16X4F32SASYM(vC1l, vC3l);
      AE_S16X4_IP(vC2s, p4_x3, 8);
      vA3s = AE_SRAI16(vA3s, FIRST_STAGE_SCALE);
      AE_S16X4_IP(vC1s, p4_x1, 8);
      AE_S16X4_IP(vC3s, p4_x2, 8);

      stride  >>= 2;
    }

    {
      int lc = (N>>3)-1;

      p4_x0 = (ae_int16x4 *)x;
      p4_x1 = p4_x0 + (stride>>3);
      p4_x2 = p4_x1 + (stride>>3);
      p4_x3 = p4_x2 + (stride>>3);

      //-----------------------------------------------------------------
      // Set up offsets to access "N/4", "N/2", "3N/4" complex point or  
      // "N/2", "N", "3N/2" half word                                    
      //-----------------------------------------------------------------
      vA0s = AE_L16X4_I(p4_x0, 0);
      vA1s = AE_L16X4_I(p4_x1, 0);
      vA2s = AE_L16X4_I(p4_x2, 0);
      vA3s = AE_L16X4_I(p4_x3, 0);

      vB2s = AE_SUB16S(vA0s, vA2s);
      vB3s = AE_SUB16S(vA1s, vA3s);
      vA0s = AE_ADD16S(vA0s, vA2s);
      vA1s = AE_ADD16S(vA1s, vA3s);

      vT0s = AE_L16X4_I(((const ae_int16x4*)fft_twd_r8), 0);
      vT1s = AE_L16X4_I(((const ae_int16x4*)fft_twd_r8), 8);
      vT2s = AE_L16X4_I(((const ae_int16x4*)fft_twd_r8), 16);

      do
      {
        vB1s = AE_SUB16S(vA0s, vA1s);
        vC3s = AE_ADD16S(vA0s, vA1s);
        vC0l = AE_CVT32X2F16_32(vB1s);
        vC1l = AE_CVT32X2F16_10(vB1s);
        vC3s = AE_SRAI16(vC3s, 1);
        vC0l = AE_SRAI32(vC0l, 1);
        AE_S16X4_IP(vC3s, p4_x0, 32);
        vC0l = AE_MULFC32X16RAS_L(vC0l, vT0s);
        vC1l = AE_SRAI32(vC1l, 1);
        vA0s = AE_L16X4_I(p4_x0, 0);
        vB1l = AE_CVT32X2F16_32(vB3s);
        vC1l = AE_MULFC32X16RAS_H(vC1l, vT2s);
        vA2s = AE_L16X4_I(p4_x2, 32);
        vB0l = AE_CVT32X2F16_32(vB2s);
        vB1l = AE_SEL32_LH(vB1l, vB1l);
        vC1s = AE_ROUND16X4F32SASYM(vC0l, vC1l);
        vB3l = AE_CVT32X2F16_10(vB3s);
        vC0l = AE_SUBADD32S(vB0l, vB1l);
        vB3l = AE_SEL32_LH(vB3l, vB3l);
        vB2l = AE_CVT32X2F16_10(vB2s);
        vC1l = AE_ADDSUB32S(vB0l, vB1l);
        vC2l = AE_SUBADD32S(vB2l, vB3l);
        vC3l = AE_ADDSUB32S(vB2l, vB3l);
        vC0l = AE_SRAI32(vC0l, 1);
        vC0l = AE_MULFC32X16RAS_H(vC0l, vT1s);
        vC2l = AE_SRAI32(vC2l, 1);
        vA1s = AE_L16X4_I(p4_x1, 32);
        vC2l = AE_MULFC32X16RAS_L(vC2l, vT2s);
        vC1l = AE_SRAI32(vC1l, 1);
        vC1l = AE_MULFC32X16RAS_H(vC1l, vT0s);
        vC3l = AE_SRAI32(vC3l, 1);
        vC3l = AE_MULFC32X16RAS_L(vC3l, vT1s);
        vA3s = AE_L16X4_I(p4_x3, 32);
        vC2s = AE_ROUND16X4F32SASYM(vC0l, vC2l);
        vB2s = AE_SUB16S(vA0s, vA2s);
        vC3s = AE_ROUND16X4F32SASYM(vC1l, vC3l);
        AE_S16X4_XP(vC2s, p4_x3, 32);
        vA0s = AE_ADD16S(vA0s, vA2s);
        AE_S16X4_IP(vC1s, p4_x1, 32);
        AE_S16X4_IP(vC3s, p4_x2, 32);
        vB3s = AE_SUB16S(vA1s, vA3s);
        vA1s = AE_ADD16S(vA1s, vA3s);
      } while(--lc);
      vB1s = AE_SUB16S(vA0s, vA1s);
      vC3s = AE_ADD16S(vA0s, vA1s);
      vC0l = AE_CVT32X2F16_32(vB1s);
      vC1l = AE_CVT32X2F16_10(vB1s);
      vC3s = AE_SRAI16(vC3s, 1);
      vC0l = AE_SRAI32(vC0l, 1);
      AE_S16X4_IP(vC3s, p4_x0, 32);
      vC0l = AE_MULFC32X16RAS_L(vC0l, vT0s);
      vC1l = AE_SRAI32(vC1l, 1);
      vB1l = AE_CVT32X2F16_32(vB3s);
      vC1l = AE_MULFC32X16RAS_H(vC1l, vT2s);
      vB0l = AE_CVT32X2F16_32(vB2s);
      vB1l = AE_SEL32_LH(vB1l, vB1l);
      vC1s = AE_ROUND16X4F32SASYM(vC0l, vC1l);
      vB3l = AE_CVT32X2F16_10(vB3s);
      vC0l = AE_SUBADD32S(vB0l, vB1l);
      vB3l = AE_SEL32_LH(vB3l, vB3l);
      vB2l = AE_CVT32X2F16_10(vB2s);
      vC1l = AE_ADDSUB32S(vB0l, vB1l);
      vC2l = AE_SUBADD32S(vB2l, vB3l);
      vC3l = AE_ADDSUB32S(vB2l, vB3l);
      vC0l = AE_SRAI32(vC0l, 1);
      vC0l = AE_MULFC32X16RAS_H(vC0l, vT1s);
      vC2l = AE_SRAI32(vC2l, 1);
      vC2l = AE_MULFC32X16RAS_L(vC2l, vT2s);
      vC1l = AE_SRAI32(vC1l, 1);
      vC1l = AE_MULFC32X16RAS_H(vC1l, vT0s);
      vC3l = AE_SRAI32(vC3l, 1);
      vC3l = AE_MULFC32X16RAS_L(vC3l, vT1s);
      vC2s = AE_ROUND16X4F32SASYM(vC0l, vC2l);
      vC3s = AE_ROUND16X4F32SASYM(vC1l, vC3l);
      AE_S16X4_XP(vC2s, p4_x3, 32);
      AE_S16X4_IP(vC1s, p4_x1, 32);
      AE_S16X4_IP(vC3s, p4_x2, 32);

      stride >>= 2;
    }

    {
      p_y0=(ae_int32 *)(y);
      p_y1=(ae_int32 *)(p_y0 + (N >> 2));
      p_y2=(ae_int32 *)(p_y1 + (N >> 2));
      p_y3=(ae_int32 *)(p_y2 + (N >> 2));         
      p_x0=(ae_int32x2 *)(x);
      int i;
      int32_t i0,i1,i2,i3,ai;
      i = NSA(N)+2;
      ai=((int32_t)0x1)<<i;
      i0=0;

      {
        //--------------------------------------------------------------------------
        // last stage is RADIX2 !!!
        //--------------------------------------------------------------------------
        for (i = 0; i < (N>>4); i++) 
        {
          i1 = AE_ADDBRBA32(i0, ai);
          i2 = AE_ADDBRBA32(i1, ai);
          i3 = AE_ADDBRBA32(i2, ai);
          FFT_BUTTERFLY_R2(i0);
          FFT_BUTTERFLY_R2(i1);
          FFT_BUTTERFLY_R2(i2);
          FFT_BUTTERFLY_R2(i3);
          i0 = AE_ADDBRBA32(i3, ai);
        }
      } 
    }
}
