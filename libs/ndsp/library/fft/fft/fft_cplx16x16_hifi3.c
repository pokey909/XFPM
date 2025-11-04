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
    C code optimized for HiFi3/HiFi3z
*/
/*===========================================================================
  Fast Fourier Transforms:
  fft_cplx             FFT on Complex Data
  fft_real             FFT on Real Data
  ifft_cplx            Inverse FFT on Complex Data
  ifft_real            Inverse FFT on Real Data
  fft_cplx<prec>_ie    FFT on Complex Data with Optimized Memory Usage
  fft_real<prec>_ie    FFT on Real Data with Optimized Memory Usage
  ifft_cplx<prec>_ie   Inverse FFT on Complex Data with Optimized Memory Usage
  ifft_real<prec>_ie   Inverse FFT on Real Data with Optimized Memory Usage
  dct/dct4             Discrete Cosine Transform, Type II/Type IV
  mdct                 Modified Discrete Cosine Transform
  imdct                Inverse Modified Discrete Cosine Transform
  dct2d                2-D Discrete Cosine Transform
  idct2d               Inverse 2-D Discrete Cosine Transform

  There are limited combinations of precision and scaling options available:
  ----------------+---------------------------------------------------------------
      FFT/IFFT    | Scaling options                        | Restrictions on the
                  |                                        | input dynamic range
  ----------------+---------------------------------------------------------------
  cplx24x24,      | 0 - no scaling                         | input signal < 2^23/(2*N),
                  |                                        | N-fft-size
  real24x24       | 1 - 24-bit scaling                     |        none
                  | 2 - 32-bit scaling on the first stage  |        none
                  | and 24-bit scaling later               |        none
                  | 3 - fixed scaling before each stage    |        none
------------------------------------------------------------------------------------
  cplx32x16       | 3 - fixed scaling before each stage    |        none
------------------------------------------------------------------------------------
  cplx16x16       | 2 - 16-bit dynamic scaling             |        none
                  | 3 - fixed scaling before each stage    |        none
  cplx32x32       | 2 - 32-bit dynamic scaling             |        none
                  | 3 - fixed scaling before each stage    |        none
  cplx32x32_ie    | 2 - 32-bit dynamic scaling             |        none
                  | 3 - fixed scaling before each stage    |        none
------------------------------------------------------------------------------------
  cplx16x16_ie    | 2 - 16-bit dynamic scaling             |        none
------------------------------------------------------------------------------------
  cplx32x16_ie    | 3 - fixed scaling before each stage    |        none
  cplx24x24_ie    | 3 - fixed scaling before each stage    |        none
  real32x16       | 3 - fixed scaling before each stage    |        none
------------------------------------------------------------------------------------
  real16x16       | 2 - 16-bit dynamic scaling             |        none
                  | 3 - fixed scaling before each stage    |        none
  real32x32       | 2 - 32-bit dynamic scaling             |        none
                  | 3 - fixed scaling before each stage    |        none
  real32x32_ie    | 2 - 32-bit dynamic scaling             |        none
                  | 3 - fixed scaling before each stage    |        none
------------------------------------------------------------------------------------
  real16x16_ie    | 2 - 16-bit dynamic scaling             |        none
------------------------------------------------------------------------------------
  real32x16_ie    | 3 - fixed scaling before each stage    |        none
  real24x24_ie    | 3 - fixed scaling before each stage    |        none
  real32x16_ie_24p| 3 - fixed scaling before each stage    |        none
  ----------------+---------------------------------------------------------------
  real24x24_ie_24p| 1 - 24-bit scaling                     |        none
  ----------------+---------------------------------------------------------------
  DCT:            |
  ----------------+---------------------------------------------------------------
  dct_24x24       | 3 - fixed scaling before each stage    |        none
  dct_32x16       | 3 - fixed scaling before each stage    |        none
  dct_32x32       | 3 - fixed scaling before each stage    |        none
  dct_16x16       | 3 - fixed scaling before each stage    |        none
  dct4_24x24      | 3 - fixed scaling before each stage    |        none
  dct4_32x16      | 3 - fixed scaling before each stage    |        none
  dct4_32x32      | 3 - fixed scaling before each stage    |        none
  mdct_24x24      | 3 - fixed scaling before each stage    |        none
  mdct_32x16      | 3 - fixed scaling before each stage    |        none
  mdct_32x32      | 3 - fixed scaling before each stage    |        none
  imdct_24x24     | 3 - fixed scaling before each stage    |        none
  imdct_32x16     | 3 - fixed scaling before each stage    |        none
  imdct_32x32     | 3 - fixed scaling before each stage    |        none
  ----------------+---------------------------------------------------------------
  dct2d           | 0 - no scaling                         |        none
  ----------------+---------------------------------------------------------------
===========================================================================*/
#include "NatureDSP_Signal_fft.h"
#include "common.h"
#include "fft_cplx_twiddles.h"

#define FIRST_STAGE_SCALE 3
   
#define DFT4XI2_HIFI3(x0, x1, x2, x3, shift)/* output x0, x3, x1, x2*/\
{\
    ae_int16x4 t0, t1, t2, t3;                                                               \
    xtbool4 mask = 0x5;                                                                      \
    x0 = AE_SRAA16RS(x0, shift);                                                             \
    x1 = AE_SRAA16RS(x1, shift);                                                             \
    x2 = AE_SRAA16RS(x2, shift);                                                             \
    x3 = AE_SRAA16RS(x3, shift);                                                             \
    t0 = AE_ADD16S(x0, x2); t2 = AE_SUB16S(x0, x2); /*AE_ADDANDSUBRNG16RAS_S1(x0, x2); */    \
    t1 = AE_ADD16S(x1, x3); t3 = AE_SUB16S(x1, x3); /*AE_ADDANDSUBRNG16RAS_S1(x1, x3); */    \
    x0 = t0; x1 = t1; x2 = t2; x3 = t3;                                                      \
     AE_MOVT16X4(x3, AE_NEG16S(x3), mask); /*  x3 = AE_CONJ16S(x3); */                        \
    x3 = AE_SHORTSWAP(x3);   x3 = AE_SEL16_5432(x3, x3);                                     \
                                                                                             \
    t0 = AE_ADD16S(x0, x1); t1 = AE_SUB16S(x0, x1); /*AE_ADDANDSUBRNG16RAS_S2(x0, x1); */    \
    t2 = AE_ADD16S(x2, x3); t3 = AE_SUB16S(x2, x3); /*AE_ADDANDSUBRNG16RAS_S2(x2, x3); */    \
    x0 = t0; x1 = t1; x2 = t2; x3 = t3;                                                      \
    t1 = x3;                                                                                 \
    t2 = x1;                                                                                 \
    t3 = x2;                                                                                 \
    x1 = t1;                                                                                 \
    x2 = t2;                                                                                 \
    x3 = t3;                                                                                 \
}


#if !(XCHAL_HAVE_HIFI3Z)

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

static ALIGN(16) const int16_t fft_twd_r8[12] = {
    (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000,
    (int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x5a82,(int16_t)0xA57E,
    (int16_t)0x0000,(int16_t)0x8000,(int16_t)0xa57e,(int16_t)0xA57E
};

static int fft_core_s3(int16_t *y, int16_t *x, const tFftDescr *pDescr)
{
  int N;
  ae_int16x4  * restrict p4_x0, * restrict p4_x1,
              * restrict p4_x2, * restrict p4_x3;
  const ae_int16x4 * restrict p_twd;
  const ae_int16 * restrict p_inc;

  ae_int32 * restrict p_y0;
  ae_int32 * restrict p_y1;
  ae_int32 * restrict p_y2;
  ae_int32 * restrict p_y3;
  ae_int32x2 * restrict p_x0;

  ae_int16x4  vS;
  ae_int16x4  vA0s, vA1s, vA2s, vA3s, vB1s, vB2s, vB3s,
              vC1s, vC2s, vC3s, vT0s, vT1s, vT2s;
  ae_int32x2  vA0l, vA1l, vA2l, vA3l, vB0l, vB1l, vB2l, vB3l, vC0l, vC1l, vC2l, vC3l;
  const cint32_ptr * restrict seq;
  int step_circ;
  int stride,scale=0;

  {
    N=pDescr->N;
    p_inc = (const ae_int16 *)(pDescr->inc);
    seq=(const cint32_ptr * )pDescr->twd;
    p_twd = (const ae_int16x4 *)*seq++;
  }
  stride     =   N; // The stride is quartered with every iteration of the outer loop.

  if (stride > 4)
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
    /* 17 cycles per pipeline stage in steady state with unroll=1,  DFT4XI2 */
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

    scale += 2;
    stride  >>= 2;
  }

  while (stride > 8) 
  {
    int lc = N>>4;
    int step;

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

    p_twd=(const ae_int16x4 *)(*seq++);
    WUR_AE_CBEGIN0( (uintptr_t)p_twd );
    WUR_AE_CEND0  ( (uintptr_t)(p_twd)+stride*3 );
    step_circ=24;
    /* 32 cycles per pipeline stage in steady state with unroll=1 */
    do
    {
      AE_L16_IP(vS, p_inc, 2);
      vS = AE_SLAI16S(vS, 1);
      step = AE_MOVAD16_0(vS)+4;
      vT2s = AE_L16X4_I(p_twd, 16);
      vB1s = AE_SUB16S(vA0s, vA1s);
      vC3s = AE_ADD16S(vA0s, vA1s);
      vT1s = AE_L16X4_I(p_twd, 8);
      AE_L16X4_XC(vT0s, p_twd, step_circ);
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
      vC2l = AE_SUBADD32S(vB2l, vB3l);
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
      vA3s = AE_L16X4_I(p4_x3, 8);
      vC2s = AE_ROUND16X4F32SASYM(vC0l, vC2l);
      vB2s = AE_SUB16S(vA0s, vA2s);
      vC3s = AE_ROUND16X4F32SASYM(vC1l, vC3l);
      AE_S16X4_IP(vC2s, p4_x3, 8);
      vA0s = AE_ADD16S(vA0s, vA2s);
      AE_S16X4_IP(vC1s, p4_x1, 8);
      AE_S16X4_IP(vC3s, p4_x2, 8);
      vB3s = AE_SUB16S(vA1s, vA3s);
      vA1s = AE_ADD16S(vA1s, vA3s);

      vT2s = AE_L16X4_I(p_twd, 16);
      vB1s = AE_SUB16S(vA0s, vA1s);
      vC3s = AE_ADD16S(vA0s, vA1s);
      vT1s = AE_L16X4_I(p_twd, 8);
      AE_L16X4_XC(vT0s, p_twd, step_circ);
      vC0l = AE_CVT32X2F16_32(vB1s);
      vC1l = AE_CVT32X2F16_10(vB1s);
      vC3s = AE_SRAI16(vC3s, 2);
      vC0l = AE_SRAI32(vC0l, 2);
      AE_S16X4_XP(vC3s, p4_x0, step);
      vC0l = AE_MULFC32X16RAS_L(vC0l, vT0s);
      vC1l = AE_SRAI32(vC1l, 2);
      vA0s = AE_L16X4_I(p4_x0, 0);
      vB1l = AE_CVT32X2F16_32(vB3s);
      vC1l = AE_MULFC32X16RAS_H(vC1l, vT2s);
      vA2s = AE_L16X4_X(p4_x2, step);
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
      vC0l = AE_SRAI32(vC0l, 2);
      vC0l = AE_MULFC32X16RAS_H(vC0l, vT1s);
      vC2l = AE_SRAI32(vC2l, 2);
      vA1s = AE_L16X4_X(p4_x1, step);
      vC2l = AE_MULFC32X16RAS_L(vC2l, vT2s);
      vC1l = AE_SRAI32(vC1l, 2);
      vC1l = AE_MULFC32X16RAS_H(vC1l, vT0s);
      vC3l = AE_SRAI32(vC3l, 2);
      vC3l = AE_MULFC32X16RAS_L(vC3l, vT1s);
      vA3s = AE_L16X4_X(p4_x3, step);
      vC2s = AE_ROUND16X4F32SASYM(vC0l, vC2l);
      vB2s = AE_SUB16S(vA0s, vA2s);
      vC3s = AE_ROUND16X4F32SASYM(vC1l, vC3l);
      AE_S16X4_XP(vC2s, p4_x3, step);
      vA0s = AE_ADD16S(vA0s, vA2s);
      AE_S16X4_XP(vC1s, p4_x1, step);
      AE_S16X4_XP(vC3s, p4_x2, step);
      vB3s = AE_SUB16S(vA1s, vA3s);
      vA1s = AE_ADD16S(vA1s, vA3s);
    } while(--lc);

    scale += 2;
    stride >>= 2;
  }

  if (stride==8) 
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
    /* 15 cycles per pipeline stage in steady state with unroll=1  */
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
    scale += 1;
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
        i = AE_NSAZ32_L(N)+2;
        ai=((int32_t)0x1)<<i;
        i0=0;
        if (stride==2) 
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
        else 
        {
            vA3l=AE_MOVDA32X2(0,0xFFFFFFFF);
            //--------------------------------------------------------------------------
            // last stage is RADIX4 !!!
            //--------------------------------------------------------------------------
           for (i = 0; i < (N>>4); i++) 
           {
                NASSERT_ALIGN8(p_y0);
                NASSERT_ALIGN8(p_y1);
                NASSERT_ALIGN8(p_y2);
                NASSERT_ALIGN8(p_y3);
                i1 = AE_ADDBRBA32(i0, ai);
                i2 = AE_ADDBRBA32(i1, ai);
                i3 = AE_ADDBRBA32(i2, ai);
                NASSERT_ALIGN8(p_x0);
                FFT_BUTTERFLY_R4(i0);
                FFT_BUTTERFLY_R4(i1);
                FFT_BUTTERFLY_R4(i2);
                FFT_BUTTERFLY_R4(i3);
                i0 = AE_ADDBRBA32(i3, ai);
            }
        }
    }

  return scale;
}

#endif //#if !(XCHAL_HAVE_HIFI3Z)




#if (XCHAL_HAVE_HIFI3Z)

#ifdef COMPILER_XTENSA
#define ATTRIBUTE_ALWAYS_INLINE __attribute__((always_inline))
#define ATTRIBUTE_NEVER_INLINE  __attribute__((noinline))
#else
#define ATTRIBUTE_ALWAYS_INLINE
#define ATTRIBUTE_NEVER_INLINE
#endif

#define FIRST_STAGE_SCALE 3

#define DFT4XI2(x0, x1, x2, x3)/* output x0, x3, x1, x2*/\
{\
    ae_int16x4 t1, t2, t3;\
    AE_ADDANDSUBRNG16RAS_S1(x0, x2);\
    AE_ADDANDSUBRNG16RAS_S1(x1, x3);\
\
    x3 = AE_MOVINT16X4_FROMF16X4(AE_MUL16JS(AE_MOVF16X4_FROMINT16X4(x3)));\
\
    AE_ADDANDSUBRNG16RAS_S2(x0, x1);\
    AE_ADDANDSUBRNG16RAS_S2(x2, x3);\
\
    t1 = x3; \
    t2 = x1; \
    t3 = x2; \
    x1 = t1; \
    x2 = t2;\
    x3 = t3;\
}


inline_ void MULCx2(ae_int16x4 *x, const ae_int16x4 *t)
{
    ae_f16x4 f;
    f = AE_MULFC16RAS(AE_MOVF16X4_FROMINT16X4(*x), AE_MOVF16X4_FROMINT16X4(*t));   
    *x = AE_MOVINT16X4_FROMF16X4(f);
}

#if 1
static void stage_first_DFT4xI2_inpl(int16_t *x, int16_t *twd, int N)
{
    ae_int16x4 * restrict px = (ae_int16x4 *)x;
    ae_int16x4 * restrict py = (ae_int16x4 *)x;
    ae_int32x2 * restrict ptw32x2 = (ae_int32x2 *)twd;
    ae_int32x2 r1, r2, r3;
    int stride = N / 4 * sizeof(int16_t) * 2;

    int i = N >> 3;
    ae_int16x4 x0, x1, x2, x3;
    ae_int16x4 w1, w2, w3;

    WUR_AE_SAR(((FIRST_STAGE_SCALE - 1) << 1) + 1);
    __Pragma("loop_count min=2");
    
    /* hifi3z: 7 cycles per pipeline stage in steady state with unroll=1*/
    for (i = 0; i < (N >> 3); i++)
    {

        x1 = AE_L16X4_X(px, stride);
        x3 = AE_L16X4_X(px, stride * 3);
        x2 = AE_L16X4_X(px, stride * 2);
        AE_L16X4_XP(x0, px, 8);
        /*zzzzzzzzzzzzzzzzzzzzzzzzzzzzzz HiFi3Z version zzzzzzzzzzzzzzzzzzzzzzzzzzz*/
        DFT4XI2(x0, x1, x2, x3);

        /* Use AE_L32X2_XP for reading twiddles ! */
        AE_L32X2_XP(r1, ptw32x2, 8);
        AE_L32X2_XP(r2, ptw32x2, 8);
        AE_L32X2_XP(r3, ptw32x2, 8);

        w1 = AE_MOVINT16X4_FROMF32X2(r1);
        w2 = AE_MOVINT16X4_FROMF32X2(r2);
        w3 = AE_MOVINT16X4_FROMF32X2(r3);

        MULCx2(&x1, &w1);
        MULCx2(&x2, &w2);
        MULCx2(&x3, &w3);


        AE_S16X4RNG_XP(x0, py, stride);
        AE_S16X4RNG_XP(x2, py, stride);
        AE_S16X4RNG_XP(x1, py, stride);
        AE_S16X4RNG_XP(x3, py, 8 - 3 * stride);
    }//for (i = 0; i < (N >> 3); i++)
} //stage_first_DFT4xI2_inpl
#else

static void stage_first_DFT4xI2_inpl(int16_t *x, int16_t *twd, int N)
{
    int stride = N / 4 * sizeof(int16_t)* 2;
    int i;
    ae_int16x4 * restrict px = (ae_int16x4 *)x;
    ae_int16x4 * restrict py = (ae_int16x4 *)x;
    ae_p16x2s * restrict  ptw16x2 = (ae_p16x2s *)twd;
    ae_int16x4 x0, x1, x2, x3;
    ae_int32x2  t10, t20, t30, t11, t21, t31; 

    __Pragma("loop_count min=2");
    for (i = 0; i < (N >> 3); i++)
    {
        x1 = AE_L16X4_X(px, stride);
        x3 = AE_L16X4_X(px, stride * 3);
        x2 = AE_L16X4_X(px, stride * 2);
        AE_L16X4_XP(x0, px, 8);

        DFT4XI2_HIFI3(x0, x1, x2, x3, 3); 
        /* ============================ Hifi3 version ====================================*/
        t10 =  AE_L16X2M_X(ptw16x2, 0);
        t20 =  AE_L16X2M_X(ptw16x2, 8);
        t30 =  AE_L16X2M_X(ptw16x2, 16);
        t11 =  AE_L16X2M_X(ptw16x2, 4);
        t21 =  AE_L16X2M_X(ptw16x2, 12);
        t31 =  AE_L16X2M_X(ptw16x2, 20);
        ptw16x2 += 6;

        t10 = AE_SLAI32(t10, 8);
        t20 = AE_SLAI32(t20, 8);
        t30 = AE_SLAI32(t30, 8);
        t11 = AE_SLAI32(t11, 8);
        t21 = AE_SLAI32(t21, 8);
        t31 = AE_SLAI32(t31, 8);

        t10 = AE_SEL32_LH(t10, t10);
        t20 = AE_SEL32_LH(t20, t20);
        t30 = AE_SEL32_LH(t30, t30);
        t11 = AE_SEL32_LH(t11, t11);
        t21 = AE_SEL32_LH(t21, t21);
        t31 = AE_SEL32_LH(t31, t31);

        t10 = AE_MULFC32X16RAS_H(t10,  x1); 
        t11 = AE_MULFC32X16RAS_L(t11,  x1);
        t20 = AE_MULFC32X16RAS_H(t20,  x2);
        t21 = AE_MULFC32X16RAS_L(t21,  x2);
        t30 = AE_MULFC32X16RAS_H(t30,  x3);
        t31 = AE_MULFC32X16RAS_L(t31,  x3);

        x1 = AE_ROUND16X4F32SASYM(t10, t11);
        x2 = AE_ROUND16X4F32SASYM(t20, t21);
        x3 = AE_ROUND16X4F32SASYM(t30, t31);

        AE_S16X4_XP(x0, py, stride);
        AE_S16X4_XP(x2, py, stride);
        AE_S16X4_XP(x1, py, stride);
        AE_S16X4_XP(x3, py, 8 - 3 * stride);
    }//for (i = 0; i < (N >> 3); i++)
} //stage_first_DFT4xI2_inpl
#endif



static void stage_last_DFT4xI2_inpl(int16_t *x, int16_t *y, int N)
{

    int i = N >> 3;
    const int shift = 2;
    ae_int16x4 * restrict px0 = (ae_int16x4 *)x;

    WUR_AE_SAR(((shift - 1) << 1) + 1);

    ae_int16x4 * restrict y10 = (ae_int16x4 *)(y);
    ae_int16x4 * restrict y11 = (ae_int16x4 *)(y + N );
    ae_int16x4 * restrict y12 = (ae_int16x4 *)(y + N / 2);
    ae_int16x4 * restrict y13 = (ae_int16x4 *)(y + 3 * N / 2);

    int ibitrev = 0;
    int ibitrev_off = 0x80000000 >> (29 - NSA(N));
    __Pragma("loop_count min=2");
    for (i=0;i<(N/4);i+=2)
    {
        ae_int16x4 e0, e1, e3, e4, e2, e5, h2, h3, h4, h5;
        
        /*Hifi3z: 13 cycles per pipeline stage in steady state with unroll=2 */
        e0 = px0[i];
        e1 = px0[i+1];
        e3 = px0[i+N/4];
        e4 = px0[i+1+N/4];

#if 1
        // Hifi3z zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
        AE_ADDANDSUBRNG16RAS_S1(e0, e1);
        e2 = e1; 

        AE_ADDANDSUBRNG16RAS_S1(e3, e4);
        e5 = e4; 

        h4 = AE_SEL16_7632(e0, e3);
        h5 = AE_SEL16_5410(e0, e3);
        AE_ADDANDSUBRNG16RAS_S2(h4, h5);

        h3 = AE_SEL16_7632(e1, e4);
        h2 = AE_SEL16_5410(e2, e5);
        h2 = AE_MUL16JS(h2);
        AE_ADDANDSUBRNG16RAS_S2(h3, h2);

        /* Store the 2 bfly outputs together  */
        AE_S16X4RNG_X(h4 , y10, ibitrev);
        AE_S16X4RNG_X(h5 , y11, ibitrev);
        AE_S16X4RNG_X(h2 , y12, ibitrev);
        AE_S16X4RNG_X(h3 , y13, ibitrev);
#else
        // ===================== Hifi3 ===========================
        ae_int16x4 tmp; 

        e0 = AE_SRAA16RS(e0, 2);
        e1 = AE_SRAA16RS(e1, 2);
        e3 = AE_SRAA16RS(e3, 2);
        e4 = AE_SRAA16RS(e4, 2);

        tmp = AE_ADD16S(e0, e1);            //AE_ADDANDSUBRNG16RAS_S1(e0, e1);
        e1 =  AE_SUB16S(e0, e1);
        e0 = tmp;

        e2 = e1;

        tmp = AE_ADD16S(e3, e4);            //AE_ADDANDSUBRNG16RAS_S1(e3, e4);
        e4  = AE_SUB16S(e3, e4);
        e3 = tmp;

        e5 = e4;

        h4 = AE_SEL16_7632(e0, e3);
        h5 = AE_SEL16_5410(e0, e3);

        tmp = AE_ADD16S(h4, h5);            //AE_ADDANDSUBRNG16RAS_S2(h4, h5);
        h5 =  AE_SUB16S(h4, h5);
        h4 = tmp;

        h3 = AE_SEL16_7632(e1, e4);
        h2 = AE_SEL16_5410(e2, e5);

        h2 = AE_SHORTSWAP(AE_CONJ16S(h2));  // h2 = AE_MUL16JS(h2);
        h2 = AE_SEL16_5432(h2, h2);

        tmp = AE_ADD16S(h3, h2);            //AE_ADDANDSUBRNG16RAS_S2(h3, h2);
        h2  = AE_SUB16S(h3, h2);
        h3 = tmp;

        /* Store the 2 bfly outputs together  */
        AE_S16X4_X(h4, y10, ibitrev);
        AE_S16X4_X(h5, y11, ibitrev);
        AE_S16X4_X(h2, y12, ibitrev);
        AE_S16X4_X(h3, y13, ibitrev);       
#endif

        ibitrev = AE_ADDBRBA32(ibitrev,ibitrev_off);
    }

} //stage_last_DFT4xI2_inpl




/*
    Inner stage with static scaling
*/
static void stage_inner_DFT4xI2_inpl_merged(int16_t *x, int16_t *y, const int16_t *twd, int N, int stride)
{
    const int shift = 2;
    ae_int16x4 * restrict px = (ae_int16x4 *)x;
    ae_int16x4 * restrict py = (ae_int16x4 *)y;
    ae_int32x2 * restrict ptw = (ae_int32x2 *)twd;

    uint32_t flag = 0;
    uint32_t flag_inc = 0x80000000 >> (NSA(stride) - NSA(N) - 1);

    int i;
    ae_int16x4 x0, x1, x2, x3;// , t;
    ae_int32x2 r1, r2, r3;
    ae_int16x4 w1, w2, w3;

    WUR_AE_CBEGIN0((unsigned)x);
    WUR_AE_CEND0((unsigned)(x + 2 * N - 4));
    WUR_AE_SAR(((shift - 1) << 1) + 1);

    __Pragma("loop_count min=2");
    for (i = 0; i < N / 8; i++)
    {
     
        /*Hifi3:   19 cycles per pipeline stage in steady state with unroll=2 */
        int tw_inc = 0;
        flag += flag_inc;
        XT_MOVEQZ(tw_inc, 24, flag);

        r2 = AE_L32X2_I(ptw, 8);
        r3 = AE_L32X2_I(ptw, 16);
        AE_L32X2_XP(r1, ptw, tw_inc);

        w1 = AE_MOVINT16X4_FROMF32X2(r1);
        w2 = AE_MOVINT16X4_FROMF32X2(r2);
        w3 = AE_MOVINT16X4_FROMF32X2(r3);

        AE_L16X4_XP(x0, px, stride);
        AE_L16X4_XP(x1, px, stride);
        AE_L16X4_XP(x2, px, stride);
        AE_L16X4_XC(x3, px, stride);

        DFT4XI2(x0, x1, x2, x3);

        MULCx2(&x1, &w1);
        MULCx2(&x2, &w2);
        MULCx2(&x3, &w3);

        AE_S16X4RNG_XP(x0, py, stride);
        AE_S16X4RNG_XP(x2, py, stride);
        AE_S16X4RNG_XP(x1, py, stride);
        AE_S16X4_XC(x3, py, stride);
    }
} //stage_inner_DFT4xI2_inpl_merged


#if 1
static void stage_inner_DFT4xI2_inpl_ref(int16_t *x, int16_t *y, const int16_t *twd, int N, int stride)
 {
     const int shift = 2;
     ae_int16x4 * restrict px = (ae_int16x4 *)x;
     ae_int16x4 * restrict py = (ae_int16x4 *)y;
     ae_int32x2 * restrict ptw = (ae_int32x2 *)twd;

     int M = 1 << (NSA(stride) - NSA(N)); // (N / stride);

     int i;
     int j;

     ae_int16x4 x0, x1, x2, x3;
     ae_int32x2 r1, r2, r3;
     ae_int16x4 w1, w2, w3;

     WUR_AE_SAR(((shift - 1) << 1) + 1);
     __Pragma("loop_count min=2"); 
    for (i = 0; i < stride /8; i++)
    {   /* ~ 20 cycles */
        py = (ae_int16x4 *)(8 * i + (uintptr_t)y);
        px = (ae_int16x4 *)(8 * i + (uintptr_t)x);

        AE_L32X2_XP(r1, ptw, 8);
        AE_L32X2_XP(r2, ptw, 8);
        AE_L32X2_XP(r3, ptw, 8);

        w1 = AE_MOVINT16X4_FROMF32X2(r1);
        w2 = AE_MOVINT16X4_FROMF32X2(r2);
        w3 = AE_MOVINT16X4_FROMF32X2(r3);

        /* Hifi3z: 7 cycles per pipeline stage in steady state with unroll=1 */
        __Pragma("loop_count min=2");
        for (j = 0; j < M ; j++)
        {

            AE_L16X4_XP(x0, px, stride);
            AE_L16X4_XP(x1, px, stride);
            AE_L16X4_XP(x2, px, stride);
            AE_L16X4_XP(x3, px, stride);
            /* zzzzzzzzzzzzzzzzzzzz Hifi3z zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz */
            DFT4XI2(x0, x1, x2, x3);

            MULCx2(&x1, &w1);
            MULCx2(&x2, &w2);
            MULCx2(&x3, &w3);

            AE_S16X4RNG_XP(x0, py, stride);
            AE_S16X4RNG_XP(x2, py, stride);
            AE_S16X4RNG_XP(x1, py, stride);
            AE_S16X4RNG_XP(x3, py, stride);
        }
    }
 } //stage_inner_DFT4xI2_inpl_ref
#else
static void stage_inner_DFT4xI2_inpl_ref(int16_t *x, int16_t *y, const int16_t *twd, int N, int stride)
{
    const int shift = 2;
    ae_int16x4 * restrict px = (ae_int16x4 *)x;
    ae_int16x4 * restrict py = (ae_int16x4 *)y;
    
    ae_p16x2s * restrict  ptw16x2 = (ae_p16x2s *)twd;

    int M = 1 << (NSA(stride) - NSA(N)); // (N / stride);

    int i;
    int j;

    ae_int16x4 x0, x1, x2, x3;
    ae_int32x2  t10, t20, t30, t11, t21, t31;
    ae_int32x2  x10, x20, x30, x11, x21, x31;

    WUR_AE_SAR(((shift - 1) << 1) + 1);
    __Pragma("loop_count min=2");
    for (i = 0; i < stride / 8; i++)
    {   /* ~ 20 cycles */
        py = (ae_int16x4 *)(8 * i + (uintptr_t)y);
        px = (ae_int16x4 *)(8 * i + (uintptr_t)x);

        t10 = AE_L16X2M_X(ptw16x2, 0);
        t20 = AE_L16X2M_X(ptw16x2, 8);
        t30 = AE_L16X2M_X(ptw16x2, 16);
        t11 = AE_L16X2M_X(ptw16x2, 4);
        t21 = AE_L16X2M_X(ptw16x2, 12);
        t31 = AE_L16X2M_X(ptw16x2, 20);
        ptw16x2 += 6;

        t10 = AE_SLAI32(t10, 8);
        t20 = AE_SLAI32(t20, 8);
        t30 = AE_SLAI32(t30, 8);
        t11 = AE_SLAI32(t11, 8);
        t21 = AE_SLAI32(t21, 8);
        t31 = AE_SLAI32(t31, 8);

        t10 = AE_SEL32_LH(t10, t10);
        t20 = AE_SEL32_LH(t20, t20);
        t30 = AE_SEL32_LH(t30, t30);
        t11 = AE_SEL32_LH(t11, t11);
        t21 = AE_SEL32_LH(t21, t21);
        t31 = AE_SEL32_LH(t31, t31);

        /* Hifi3z: 7 cycles per pipeline stage in steady state with unroll=1 */
        __Pragma("loop_count min=2");
        for (j = 0; j < M; j++)
        {

            AE_L16X4_XP(x0, px, stride);
            AE_L16X4_XP(x1, px, stride);
            AE_L16X4_XP(x2, px, stride);
            AE_L16X4_XP(x3, px, stride);

            /*======================= Hifi3 ==============================*/
            DFT4XI2_HIFI3(x0, x1, x2, x3, 2);

            x10 = AE_MULFC32X16RAS_H(t10, x1);
            x11 = AE_MULFC32X16RAS_L(t11, x1);
            x20 = AE_MULFC32X16RAS_H(t20, x2);
            x21 = AE_MULFC32X16RAS_L(t21, x2);
            x30 = AE_MULFC32X16RAS_H(t30, x3);
            x31 = AE_MULFC32X16RAS_L(t31, x3);

            x1 = AE_ROUND16X4F32SASYM(x10, x11);
            x2 = AE_ROUND16X4F32SASYM(x20, x21);
            x3 = AE_ROUND16X4F32SASYM(x30, x31);

            AE_S16X4_XP(x0, py, stride);
            AE_S16X4_XP(x2, py, stride);
            AE_S16X4_XP(x1, py, stride);
            AE_S16X4_XP(x3, py, stride);
        }
    }
} //stage_inner_DFT4xI2_inpl_ref

#endif
/*
    Inner stage with dynamic scaling
*/



#if 1
static void stage_last_DFT2xI2_inpl(int16_t *x, int16_t *y, int N)
{
    ae_int16x4 * restrict px = (ae_int16x4 *)x;
    int i;
    WUR_AE_SAR(1);

    ae_int16x4  * restrict y10,
        *restrict y11,
        *restrict y12,
        *restrict y13;
    ae_int16x4 h2, h3, h4, h5;
    ae_int16x4 e0, e1, e2, e3;
    int ibitrev_off = 0x80000000 >> (29 - NSA(N));
    int ibitrev = 0;

    y10 = (ae_int16x4 *)(y);
    y11 = (ae_int16x4 *)(y + N);
    y12 = (ae_int16x4 *)(y + N / 2);
    y13 = (ae_int16x4 *)(y + 3 * N / 2);

    //last stage with bit-reverse
    __Pragma("loop_count min=2");
    for (i = 0; i<(N / 4); i += 2)
    {
       
        /*Hifi3z:  6 cycles per pipeline stage in steady state with unroll=1 */
        e0 = px[i];
        e1 = px[i + 1];
        e2 = px[i + N / 4];
        e3 = px[i + 1 + N / 4];

        h5 = AE_SEL16_5410(e0, e2);
        h4 = AE_SEL16_7632(e0, e2);
        /*zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz HIFI3z zzzzzzzzzzzzzzzzzzzzzzzzzz */
        AE_ADDANDSUBRNG16RAS_S1(h4, h5);

        h3 = AE_SEL16_5410(e1, e3);
        h2 = AE_SEL16_7632(e1, e3);

        AE_ADDANDSUBRNG16RAS_S1(h2, h3);

        /* Store the 2 bfly outputs together */
        AE_S16X4RNG_X(h4, y10, ibitrev);
        AE_S16X4RNG_X(h5, y11, ibitrev);
        AE_S16X4RNG_X(h2, y12, ibitrev);
        AE_S16X4RNG_X(h3, y13, ibitrev);

        ibitrev = AE_ADDBRBA32(ibitrev, ibitrev_off);
    }
} //stage_last_DFT2xI2_inpl
#else
static void stage_last_DFT2xI2_inpl(int16_t *x, int16_t *y, int N)
{
    ae_int16x4 * restrict px = (ae_int16x4 *)x;
    int i;
    WUR_AE_SAR(1);

    ae_int16x4  * restrict y10,
        *restrict y11,
        *restrict y12,
        *restrict y13;
    ae_int16x4 h2, h3, h4, h5;
    ae_int16x4 e0, e1, e2, e3;
    int ibitrev_off = 0x80000000 >> (29 - NSA(N));
    int ibitrev = 0;

    y10 = (ae_int16x4 *)(y);
    y11 = (ae_int16x4 *)(y + N);
    y12 = (ae_int16x4 *)(y + N / 2);
    y13 = (ae_int16x4 *)(y + 3 * N / 2);

    //last stage with bit-reverse
    __Pragma("loop_count min=2");
    for (i = 0; i<(N / 4); i += 2)
    {
        ae_int16x4 tmp; 

        e0 = px[i];
        e1 = px[i + 1];
        e2 = px[i + N / 4];
        e3 = px[i + 1 + N / 4];

        e0 = AE_SRAA16RS(e0, 1);
        e1 = AE_SRAA16RS(e1, 1);
        e2 = AE_SRAA16RS(e2, 1);
        e3 = AE_SRAA16RS(e3, 1);

        h5 = AE_SEL16_5410(e0, e2);
        h4 = AE_SEL16_7632(e0, e2);
        /*========================= Hifi3 ====================================*/
        tmp = AE_ADD16S(h4, h5); //AE_ADDANDSUBRNG16RAS_S1(h4, h5);
        h5  = AE_SUB16S(h4, h5); 
        h4 = tmp; 

        h3 = AE_SEL16_5410(e1, e3);
        h2 = AE_SEL16_7632(e1, e3);

        tmp = AE_ADD16S(h2, h3);    //AE_ADDANDSUBRNG16RAS_S1(h2, h3);
        h3  = AE_SUB16S(h2, h3);
        h2 = tmp; 

        /* Store the 2 bfly outputs together */
        AE_S16X4_X(h4, y10, ibitrev);
        AE_S16X4_X(h5, y11, ibitrev);
        AE_S16X4_X(h2, y12, ibitrev);
        AE_S16X4_X(h3, y13, ibitrev);

        ibitrev = AE_ADDBRBA32(ibitrev, ibitrev_off);
    }
} //stage_last_DFT2xI2_inpl
#endif
#endif

#if (XCHAL_HAVE_HIFI3Z)
 static  int stage_second_DFT4xI2_inpl_ref_ds(int16_t *x, int16_t *y, const int16_t *twd, int N, int stride)
{
    int shift;
    ae_int16x4 * restrict px = (ae_int16x4 *)x;
    ae_int16x4 * restrict py = (ae_int16x4 *)y;
    ae_int32x2 * restrict ptw = (ae_int32x2 *)twd;
    int i;

    ae_int16x4 x0, x1, x2, x3;
    ae_int32x2 r1, r2, r3;
    ae_int16x4 w1, w2, w3;

    shift = AE_CALCRNG3();
    /* 28 cycles per pipeline stage in steady state with unroll=1 */
    __Pragma("loop_count min=1");
    for (i = 0; i < stride / 8; i++)
    {   /* ~ 20 cycles */
        //py = (ae_int16x4 *)(8 * i + (uintptr_t)y);
        //px = (ae_int16x4 *)(8 * i + (uintptr_t)x);

        AE_L32X2_XP(r1, ptw, 8);
        AE_L32X2_XP(r2, ptw, 8);
        AE_L32X2_XP(r3, ptw, 8);

        w1 = AE_MOVINT16X4_FROMF32X2(r1);
        w2 = AE_MOVINT16X4_FROMF32X2(r2);
        w3 = AE_MOVINT16X4_FROMF32X2(r3);
        
        
        /* Hifi3z: 7 cycles per pipeline stage in steady state with unroll=1 */
        //__Pragma("loop_count min=2");
        //for (j = 0; j < M; j++)
        {
            AE_L16X4_XP(x0, px, stride);
            AE_L16X4_XP(x1, px, stride);
            AE_L16X4_XP(x2, px, stride);
            AE_L16X4_XP(x3, px, stride);

            DFT4XI2(x0, x1, x2, x3);
            /* zzzzzzzzzzzzzzzzzzzzzzz Hifi3z zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz*/
            MULCx2(&x1, &w1);
            MULCx2(&x2, &w2);
            MULCx2(&x3, &w3);

            AE_S16X4RNG_XP(x0, py, stride);
            AE_S16X4RNG_XP(x2, py, stride);
            AE_S16X4RNG_XP(x1, py, stride);
            AE_S16X4RNG_XP(x3, py, stride);
        }
        {

            AE_L16X4_XP(x0, px, stride);
            AE_L16X4_XP(x1, px, stride);
            AE_L16X4_XP(x2, px, stride);
            AE_L16X4_XP(x3, px, stride);

            DFT4XI2(x0, x1, x2, x3);
            /* zzzzzzzzzzzzzzzzzzzzzzz Hifi3z zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz*/
            MULCx2(&x1, &w1);
            MULCx2(&x2, &w2);
            MULCx2(&x3, &w3);

            AE_S16X4RNG_XP(x0, py, stride);
            AE_S16X4RNG_XP(x2, py, stride);
            AE_S16X4RNG_XP(x1, py, stride);
            AE_S16X4RNG_XP(x3, py, stride);
        }
        {

            AE_L16X4_XP(x0, px, stride);
            AE_L16X4_XP(x1, px, stride);
            AE_L16X4_XP(x2, px, stride);
            AE_L16X4_XP(x3, px, stride);

            DFT4XI2(x0, x1, x2, x3);
            /* zzzzzzzzzzzzzzzzzzzzzzz Hifi3z zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz*/
            MULCx2(&x1, &w1);
            MULCx2(&x2, &w2);
            MULCx2(&x3, &w3);

            AE_S16X4RNG_XP(x0, py, stride);
            AE_S16X4RNG_XP(x2, py, stride);
            AE_S16X4RNG_XP(x1, py, stride);
            AE_S16X4RNG_XP(x3, py, stride);
        }
        {

            AE_L16X4_XP(x0, px, stride);
            AE_L16X4_XP(x1, px, stride);
            AE_L16X4_XP(x2, px, stride);
            AE_L16X4_XP(x3, px, sizeof(x0)-15 * stride);

            DFT4XI2(x0, x1, x2, x3);
            /* zzzzzzzzzzzzzzzzzzzzzzz Hifi3z zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz*/
            MULCx2(&x1, &w1);
            MULCx2(&x2, &w2);
            MULCx2(&x3, &w3);

            AE_S16X4RNG_XP(x0, py, stride);
            AE_S16X4RNG_XP(x2, py, stride);
            AE_S16X4RNG_XP(x1, py, stride);
            AE_S16X4RNG_XP(x3, py, sizeof(x0)-15 * stride);
        }
    }
    return shift;
} //stage_second_DFT4xI2_inpl_ref_ds
#endif

#if (XCHAL_HAVE_HIFI3Z)
static int stage_inner_DFT4xI2_inpl_ref_ds(int16_t *x, int16_t *y, const int16_t *twd, int N, int stride)
{
    int shift;
    ae_int16x4 * restrict px = (ae_int16x4 *)x;
    ae_int16x4 * restrict py = (ae_int16x4 *)y;
    ae_int32x2 * restrict ptw = (ae_int32x2 *)twd;

    int M = 1 << (NSA(stride) - NSA(N)); // (N / stride);
    int i;
    int j;
    ae_int16x4 x0, x1, x2, x3;
    ae_int32x2 r1, r2, r3;
    ae_int16x4 w1, w2, w3;

    shift = AE_CALCRNG3();
 

    __Pragma("loop_count min=2");
    for (i = 0; i < stride / 8; i++)
    {   /* ~ 20 cycles */
        py = (ae_int16x4 *)(8 * i + (uintptr_t)y);
        px = (ae_int16x4 *)(8 * i + (uintptr_t)x);

        AE_L32X2_XP(r1, ptw, 8);
        AE_L32X2_XP(r2, ptw, 8);
        AE_L32X2_XP(r3, ptw, 8);

        w1 = AE_MOVINT16X4_FROMF32X2(r1);
        w2 = AE_MOVINT16X4_FROMF32X2(r2);
        w3 = AE_MOVINT16X4_FROMF32X2(r3);

        ASSERT(M>=4); 
        /* Hifi3z: 7 cycles per pipeline stage in steady state with unroll=1 */
        __Pragma("loop_count min=2");
        for (j = 0; j < M; j++)
        {

            AE_L16X4_XP(x0, px, stride);
            AE_L16X4_XP(x1, px, stride);
            AE_L16X4_XP(x2, px, stride);
            AE_L16X4_XP(x3, px, stride);

            DFT4XI2(x0, x1, x2, x3);
            /* zzzzzzzzzzzzzzzzzzzzzzz Hifi3z zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz*/
            MULCx2(&x1, &w1);
            MULCx2(&x2, &w2);
            MULCx2(&x3, &w3);

            AE_S16X4RNG_XP(x0, py, stride);
            AE_S16X4RNG_XP(x2, py, stride);
            AE_S16X4RNG_XP(x1, py, stride);
            AE_S16X4RNG_XP(x3, py, stride);
        }
    }
    return shift;
} //stage_inner_DFT4xI2_inpl_ref_ds
#else  //#if (XCHAL_HAVE_HIFI3Z)
static int stage_inner_DFT4xI2_inpl_ref_ds(int16_t *x, int16_t *y, const int16_t *twd, int N, int stride)
{
    int shift;
    ae_int16x4 * restrict px0 = (ae_int16x4 *)x;
    ae_int16x4 * restrict py = (ae_int16x4 *)y;
    ae_int32x2 * restrict px1 = (ae_int32x2 *)x;

    ae_p16x2s * restrict  ptw16x2 = (ae_p16x2s *)twd;

    int M = 1 << (NSA(stride) - NSA(N)); // (N / stride);

    int i;
    int j;

    ae_int16x4 x0, x1, x2, x3;
    ae_int32x2  t10, t20, t30, t11, t21, t31;
    ae_int32x2  x10, x20, x30, x11, x21, x31;
    ae_int16x4 acc16 = AE_MOVINT16X4_FROMF32X2( AE_MOVI(0) );

    int sar = RUR_AE_SAR();
    int bexp = NSA((sar << (30 - 5)) | 1);

    if (bexp > 3) bexp = 3;


    shift = 3 - bexp;
    ptw16x2--; 

    __Pragma("loop_count min=2");
    for (i = 0; i < stride / 8; i++)
    {   /* ~ 20 cycles */
        py = (ae_int16x4 *)(8 * i + (uintptr_t)y);
        px0 = (ae_int16x4 *)(8 * i + (uintptr_t)x);
        px1 = (ae_int32x2 *)(8 * i + stride + (uintptr_t)x);

        AE_L16X2M_XU(t10, ptw16x2, 4);
        AE_L16X2M_XU(t20, ptw16x2, 4);
        AE_L16X2M_XU(t30, ptw16x2, 4);
        AE_L16X2M_XU(t11, ptw16x2, 4);
        AE_L16X2M_XU(t21, ptw16x2, 4);
        AE_L16X2M_XU(t31, ptw16x2, 4);
   

        t10 = AE_SLAI32(t10, 8);
        t20 = AE_SLAI32(t20, 8);
        t30 = AE_SLAI32(t30, 8);
        t11 = AE_SLAI32(t11, 8);
        t21 = AE_SLAI32(t21, 8);
        t31 = AE_SLAI32(t31, 8);



        /* 17 cycles per pipeline stage in steady state with unroll=1
           14 cycles lower bound required by resources*/
        __Pragma("loop_count min=2");
        for (j = 0; j < M; j++)
        {
            ae_int32x2 tmp1, tmp3; 
#if 1
            AE_L16X4_XP(x0, px0, 2 * stride);
            AE_L16X4_XP(x2, px0, 2 * stride);
            AE_L32X2_XP(tmp1, px1, 2 * stride);
            AE_L32X2_XP(tmp3, px1, 2 * stride);

            x1 = AE_MOVINT16X4_FROMINT32X2(tmp1); 
            x3 = AE_MOVINT16X4_FROMINT32X2(tmp3);
#else
            AE_L16X4_XP(x1, px, stride);
            AE_L16X4_XP(x3, px, stride);
            AE_L16X4_XP(x0, px, stride);
            AE_L16X4_XP(x2, px, stride);

#endif
            /*======================= Hifi3 ==============================*/
           // DFT4XI2_HIFI3(x0, x1, x2, x3, shift);
            {
                ae_int16x4 t0, t1, t2, t3;                                                               
                xtbool4 mask = 0x5;                                                                      
                x0 = AE_SRAA16RS(x0, shift);                                                             
                x1 = AE_SRAA16RS(x1, shift);                                                             
                x2 = AE_SRAA16RS(x2, shift);                                                             
                x3 = AE_SRAA16RS(x3, shift);                                                             
                t0 = AE_ADD16S(x0, x2); t2 = AE_SUB16S(x0, x2); /*AE_ADDANDSUBRNG16RAS_S1(x0, x2); */    
                t1 = AE_ADD16S(x1, x3); t3 = AE_SUB16S(x1, x3); /*AE_ADDANDSUBRNG16RAS_S1(x1, x3); */    

                x0 = t0; x1 = t1; x2 = t2; x3 = t3;   
#if 1                                                                   
                x1 = AE_SHORTSWAP(x1);
                x1 = AE_SEL16_5432(x1, x1);


                AE_MOVF16X4(x3, AE_NEG16S(x3), mask); /*  x3 = AE_CONJ16S(x3); */                        
#else
                AE_MOVT16X4(x3, AE_NEG16S(x3), mask); /*  x3 = AE_CONJ16S(x3); */                        
                x3 = AE_SEL16_5432(x3, x3);
                x3 = AE_SHORTSWAP(x3);

#endif
               
                t0 = AE_ADD16S(x0, x1); t1 = AE_SUB16S(x0, x1); /*AE_ADDANDSUBRNG16RAS_S2(x0, x1); */    
                t2 = AE_ADD16S(x2, x3); t3 = AE_SUB16S(x2, x3); /*AE_ADDANDSUBRNG16RAS_S2(x2, x3); */    
                x0 = t0; x1 = t1; x2 = t2; x3 = t3;                                                      
                t1 = x3;                                                                                 
                t2 = x1;                                                                                 
                t3 = x2;                                                                                 
                x1 = t1;                                                                                 
                x2 = t2;                                                                                 
                x3 = t3;                                                                                 
            }

            x10 = AE_MULFC32X16RAS_H(t10, x1);
            x11 = AE_MULFC32X16RAS_L(t11, x1);
            x20 = AE_MULFC32X16RAS_H(t20, x2);
            x21 = AE_MULFC32X16RAS_L(t21, x2);
            x30 = AE_MULFC32X16RAS_H(t30, x3);
            x31 = AE_MULFC32X16RAS_L(t31, x3);

            x1 = AE_ROUND16X4F32SASYM(x10, x11);
            x2 = AE_ROUND16X4F32SASYM(x20, x21);
            x3 = AE_ROUND16X4F32SASYM(x30, x31);

            AE_S16X4_XP(x0, py, stride);
            AE_S16X4_XP(x2, py, stride);
            AE_S16X4_XP(x1, py, stride);
            AE_S16X4_XP(x3, py, stride);

            x0 = AE_ABS16S(x0);
            x2 = AE_ABS16S(x2);
            x1 = AE_ABS16S(x1);
            x3 = AE_ABS16S(x3);

            acc16 = AE_OR16(acc16, x0);
            acc16 = AE_OR16(acc16, x2);
            acc16 = AE_OR16(acc16, x1);
            acc16 = AE_OR16(acc16, x3);

        }
    }

    ae_int16x4 tmp0 = AE_SEL16_5432(acc16, acc16);  //AE_INTSWAP
    acc16 = AE_OR16(acc16, tmp0);
    acc16 = AE_OR16(acc16, AE_SHORTSWAP(acc16));
    int  tmpi = AE_MOVAD16_0(acc16) >> (15 - 3);
    WUR_AE_SAR(tmpi << 3);

    return shift;
} //stage_inner_DFT4xI2_inpl_ref_ds
#endif //#if (XCHAL_HAVE_HIFI3Z) ... #else ...

/*
fist stage with  dynamic scaling
*/
#if (XCHAL_HAVE_HIFI3Z)
static int stage_first_DFT4xI2_inpl_ds(int16_t *x, int16_t *twd, int N, int bexp)
{
    ae_int16x4 * restrict px = (ae_int16x4 *)x;
    ae_int16x4 * restrict py = (ae_int16x4 *)x;
    ae_int32x2 * restrict ptw32x2 = (ae_int32x2 *)twd;
    ae_int32x2 r1, r2, r3;
    int stride = N / 4 * sizeof(int16_t)* 2;

    int i = N >> 3;
    ae_int16x4 x0, x1, x2, x3;
    ae_int16x4 w1, w2, w3;

    WUR_AE_SAR(((FIRST_STAGE_SCALE - 1) << 1) + 1);
    __Pragma("loop_count min=2");

    /* 11 cycles per pipeline stage in steady state with unroll=1*/
    for (i = 0; i < (N >> 3); i++)
    {
        x1 = AE_L16X4_X(px, stride);
        x3 = AE_L16X4_X(px, stride * 3);
        x2 = AE_L16X4_X(px, stride * 2);
        AE_L16X4_XP(x0, px, 8);

        // Normalize input data
        x0 = AE_SLAA16S(x0, bexp);
        x1 = AE_SLAA16S(x1, bexp);
        x2 = AE_SLAA16S(x2, bexp);
        x3 = AE_SLAA16S(x3, bexp);
        /*zzzzzzzzzzzzzzzzzzzz Hifi3z zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz*/
        DFT4XI2(x0, x1, x2, x3);

        /* Use AE_L32X2_XP for reading twiddles ! */
        AE_L32X2_XP(r1, ptw32x2, 8);
        AE_L32X2_XP(r2, ptw32x2, 8);
        AE_L32X2_XP(r3, ptw32x2, 8);

        w1 = AE_MOVINT16X4_FROMF32X2(r1);
        w2 = AE_MOVINT16X4_FROMF32X2(r2);
        w3 = AE_MOVINT16X4_FROMF32X2(r3);

        MULCx2(&x1, &w1);
        MULCx2(&x2, &w2);
        MULCx2(&x3, &w3);

        AE_S16X4RNG_XP(x0, py, stride);
        AE_S16X4RNG_XP(x2, py, stride);
        AE_S16X4RNG_XP(x1, py, stride);
        AE_S16X4RNG_XP(x3, py, 8 - 3 * stride);
    }//for (i = 0; i < (N >> 3); i++)

    return 3 - bexp;
} //stage_first_DFT4xI2_inpl_ds
#else //#if (XCHAL_HAVE_HIFI3Z)


static int stage_first_DFT4xI2_inpl_ds(int16_t *x, int16_t *twd, int N, int bexp)
{
    int stride = N / 4 * sizeof(int16_t)* 2;
    int i;

    ae_int16x4 * restrict px = (ae_int16x4 *)x;
    ae_int16x4 * restrict py = (ae_int16x4 *)x;
    ae_p16x2s * restrict  ptw16x2 = (ae_p16x2s *)twd;
    ae_int16x4 x0, x1, x2, x3;
    ae_int16x4 tmp0;
    ae_int16x4 acc16 = AE_MOVINT16X4_FROMF32X2(AE_MOVI(0)); 
    ae_int32x2  t10, t20, t30, t11, t21, t31;

    WUR_AE_SAR(0);

    __Pragma("loop_count min=2");
    /* 20 cycles per pipeline stage in steady state with unroll=1*/
       
       
    for (i = 0; i < (N >> 3); i++)
    {

        AE_L16X4_XP(x0, px, stride);
        AE_L16X4_XP(x1, px, stride);
        AE_L16X4_XP(x2, px, stride);
        AE_L16X4_XP(x3, px, 8-3*stride);

        DFT4XI2_HIFI3(x0, x1, x2, x3, 3-bexp);

        /* ============================ Hifi3 version ====================================*/
        t10 = AE_L16X2M_X(ptw16x2, 0);
        t20 = AE_L16X2M_X(ptw16x2, 4);
        t30 = AE_L16X2M_X(ptw16x2, 8);
        t11 = AE_L16X2M_X(ptw16x2, 12);
        t21 = AE_L16X2M_X(ptw16x2, 16);
        t31 = AE_L16X2M_X(ptw16x2, 20);
        ptw16x2 += 6;

        t10 = AE_SLAI32(t10, 8);
        t20 = AE_SLAI32(t20, 8);
        t30 = AE_SLAI32(t30, 8);
        t11 = AE_SLAI32(t11, 8);
        t21 = AE_SLAI32(t21, 8);
        t31 = AE_SLAI32(t31, 8);

        t10 = AE_MULFC32X16RAS_H(t10, x1);
        t11 = AE_MULFC32X16RAS_L(t11, x1);
        t20 = AE_MULFC32X16RAS_H(t20, x2);
        t21 = AE_MULFC32X16RAS_L(t21, x2);
        t30 = AE_MULFC32X16RAS_H(t30, x3);
        t31 = AE_MULFC32X16RAS_L(t31, x3);

        x1 = AE_ROUND16X4F32SASYM(t10, t11);
        x2 = AE_ROUND16X4F32SASYM(t20, t21);
        x3 = AE_ROUND16X4F32SASYM(t30, t31);

        AE_S16X4_XP(x0, py, stride);
        AE_S16X4_XP(x2, py, stride);
        AE_S16X4_XP(x1, py, stride);
        AE_S16X4_XP(x3, py, 8 - 3 * stride);
        /*
        tmp0  = AE_OR16( AE_ABS16S(x0), AE_ABS16S(x2) );
        tmp1  = AE_OR16( AE_ABS16S(x1), AE_ABS16S(x3) );
        tmp0  = AE_OR16( tmp0, tmp1);
        acc16 = AE_OR16(tmp0, acc16);
        */
        x0 = AE_ABS16S(x0);
        x1 = AE_ABS16S(x1);
        x2 = AE_ABS16S(x2);
        x3 = AE_ABS16S(x3);

        acc16 = AE_OR16(acc16, x0);
        acc16 = AE_OR16(acc16, x1);
        acc16 = AE_OR16(acc16, x2);
        acc16 = AE_OR16(acc16, x3);

    }//for (i = 0; i < (N >> 3); i++)

    tmp0 = AE_SEL16_5432(acc16, acc16);  //AE_INTSWAP
    acc16 = AE_OR16(acc16, tmp0);
    acc16 = AE_OR16(acc16, AE_SHORTSWAP(acc16));
    int tmpi = AE_MOVAD16_0(acc16) >> (15 - 3);


    WUR_AE_SAR(tmpi << 3);

    return 3 - bexp;
} //stage_first_DFT4xI2_inpl_ds
#endif //#if (XCHAL_HAVE_HIFI3Z)

#if (XCHAL_HAVE_HIFI3Z)
static int stage_last_DFT4xI2_inpl_ds(int16_t *x, int16_t *y, int N)
{

    int i = N >> 3;
    int shift = 0;
    ae_int16x4 * restrict px0 = (ae_int16x4 *)x;

    ae_int16x4 * restrict y10 = (ae_int16x4 *)(y);
    ae_int16x4 * restrict y11 = (ae_int16x4 *)(y + N);
    ae_int16x4 * restrict y12 = (ae_int16x4 *)(y + N / 2);
    ae_int16x4 * restrict y13 = (ae_int16x4 *)(y + 3 * N / 2);

    int ibitrev = 0;
    int ibitrev_off = 0x80000000 >> (29 - NSA(N));


    shift = AE_CALCRNG3();

    __Pragma("loop_count min=2");
    for (i = 0; i<(N / 4); i += 2)
    {
        ae_int16x4 e0, e1, e3, e4, e2, e5, h2, h3, h4, h5;
        
        /*Hifi3z: 13 cycles per pipeline stage in steady state with unroll=2 */
        e0 = px0[i];
        e1 = px0[i + 1];
        e3 = px0[i + N / 4];
        e4 = px0[i + 1 + N / 4];
        /* zzzzzzzzzzzzzzzzzzzzzzzzzzzzzz Hifi3z zzzzzzzzzzzzzzzzzzzzzzzzzz*/
        AE_ADDANDSUBRNG16RAS_S1(e0, e1);
        e2 = e1;

        AE_ADDANDSUBRNG16RAS_S1(e3, e4);
        e5 = e4;

        h4 = AE_SEL16_7632(e0, e3);
        h5 = AE_SEL16_5410(e0, e3);
        AE_ADDANDSUBRNG16RAS_S2(h4, h5);

        h3 = AE_SEL16_7632(e1, e4);
        h2 = AE_SEL16_5410(e2, e5);
        h2 = AE_MUL16JS(h2);
        AE_ADDANDSUBRNG16RAS_S2(h3, h2);

        /* Store the 2 bfly outputs together  */
        AE_S16X4RNG_X(h4, y10, ibitrev);
        AE_S16X4RNG_X(h5, y11, ibitrev);
        AE_S16X4RNG_X(h2, y12, ibitrev);
        AE_S16X4RNG_X(h3, y13, ibitrev);

        ibitrev = AE_ADDBRBA32(ibitrev, ibitrev_off);
    }
    return shift;
} //stage_last_DFT4xI2_inpl_ds
#else //#if (XCHAL_HAVE_HIFI3Z)
static int stage_last_DFT4xI2_inpl_ds(int16_t *x, int16_t *y, int N)
{

    int i = N >> 3;
    int shift = 0;
    ae_int16x4 * restrict px0 = (ae_int16x4 *)x;
    int sar = RUR_AE_SAR();
    int bexp = NSA((sar << (30 - 5)) | 1);

    if (bexp > 3) bexp = 3;
    shift = 3 - bexp;


    ae_int16x4 * restrict y10 = (ae_int16x4 *)(y);
    ae_int16x4 * restrict y11 = (ae_int16x4 *)(y + N);
    ae_int16x4 * restrict y12 = (ae_int16x4 *)(y + N / 2);
    ae_int16x4 * restrict y13 = (ae_int16x4 *)(y + 3 * N / 2);

    int ibitrev = 0;
    int ibitrev_off = 0x80000000 >> (29 - NSA(N));

    __Pragma("loop_count min=2");
    /* 13 cycles per pipeline stage in steady state with unroll=1
        10 cycles lower bound required by resources */
    for (i = 0; i<(N / 4); i += 2)
    {
        xtbool4 mask = 0x5;
        ae_int16x4 e0, e1, e3, e4, e2, e5, h2, h3, h4, h5;
        e0 = px0[i];
        e1 = px0[i + 1];
        e3 = px0[i + N / 4];
        e4 = px0[i + 1 + N / 4];
        // ===================== Hifi3 ===========================
        ae_int16x4 tmp;

        e0 = AE_SRAA16RS(e0, shift);
        e1 = AE_SRAA16RS(e1, shift);
        e3 = AE_SRAA16RS(e3, shift);
        e4 = AE_SRAA16RS(e4, shift);

        tmp = AE_ADD16S(e0, e1);            //AE_ADDANDSUBRNG16RAS_S1(e0, e1);
        e1 = AE_SUB16S(e0, e1);
        e0 = tmp;

        e2 = e1;

        tmp = AE_ADD16S(e3, e4);            //AE_ADDANDSUBRNG16RAS_S1(e3, e4);
        e4 = AE_SUB16S(e3, e4);
        e3 = tmp;

        e5 = e4;

        h4 = AE_SEL16_7632(e0, e3);
        h5 = AE_SEL16_5410(e0, e3);

        tmp = AE_ADD16S(h4, h5);            //AE_ADDANDSUBRNG16RAS_S2(h4, h5);
        h5 = AE_SUB16S(h4, h5);
        h4 = tmp;

        h3 = AE_SEL16_7632(e1, e4);
        h2 = AE_SEL16_5410(e2, e5);

        /* h2 = AE_MUL16JS(h2); */
        AE_MOVT16X4(h2, AE_NEG16S(h2), mask);  // conj(h2)
        h2 = AE_SHORTSWAP(h2);   
        h2 = AE_SEL16_5432(h2, h2);

        tmp = AE_ADD16S(h3, h2);            //AE_ADDANDSUBRNG16RAS_S2(h3, h2);
        h2 = AE_SUB16S(h3, h2);
        h3 = tmp;

        /* Store the 2 bfly outputs together  */
        AE_S16X4_X(h4, y10, ibitrev);
        AE_S16X4_X(h5, y11, ibitrev);
        AE_S16X4_X(h2, y12, ibitrev);
        AE_S16X4_X(h3, y13, ibitrev);

        ibitrev = AE_ADDBRBA32(ibitrev, ibitrev_off);
    }
    return shift;
} //stage_last_DFT4xI2_inpl_ds 

#endif //#if (XCHAL_HAVE_HIFI3Z) ... #else...


#if (XCHAL_HAVE_HIFI3Z)
static int stage_last_DFT2xI2_inpl_ds(int16_t *x, int16_t *y, int N)
{

    int sar = RUR_AE_SAR();
    ae_int16x4 * restrict px = (ae_int16x4 *)x;
    int i;

    sar = 1 & (sar >> 5);


    WUR_AE_SAR(sar);

    ae_int16x4  * restrict y10,
        *restrict y11,
        *restrict y12,
        *restrict y13;
    ae_int16x4 h2, h3, h4, h5;
    ae_int16x4 e0, e1, e2, e3;
    int ibitrev_off = 0x80000000 >> (29 - NSA(N));
    int ibitrev = 0;

    y10 = (ae_int16x4 *)(y);
    y11 = (ae_int16x4 *)(y + N);
    y12 = (ae_int16x4 *)(y + N / 2);
    y13 = (ae_int16x4 *)(y + 3 * N / 2);

    //last stage with bit-reverse
    __Pragma("loop_count min=2");
    for (i = 0; i<(N / 4); i += 2)
    {
        /*Hifi3z:  6 cycles per pipeline stage in steady state with unroll=1 */
        e0 = px[i];
        e1 = px[i + 1];
        e2 = px[i + N / 4];
        e3 = px[i + 1 + N / 4];

        h5 = AE_SEL16_5410(e0, e2);
        h4 = AE_SEL16_7632(e0, e2);
        /* zzzzzzzzzzzzzzzzzzzzzzzzzzzzzz Hifi3z zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz */
        AE_ADDANDSUBRNG16RAS_S1(h4, h5);

        h3 = AE_SEL16_5410(e1, e3);
        h2 = AE_SEL16_7632(e1, e3);

        AE_ADDANDSUBRNG16RAS_S1(h2, h3);

        /* Store the 2 bfly outputs together */
        AE_S16X4RNG_X(h4, y10, ibitrev);
        AE_S16X4RNG_X(h5, y11, ibitrev);
        AE_S16X4RNG_X(h2, y12, ibitrev);
        AE_S16X4RNG_X(h3, y13, ibitrev);

        ibitrev = AE_ADDBRBA32(ibitrev, ibitrev_off);
    }
    return sar;
} //stage_last_DFT2xI2_inpl_ds
#else //#if (XCHAL_HAVE_HIFI3Z)
static int stage_last_DFT2xI2_inpl_ds(int16_t *x, int16_t *y, int N)
{

    int sar = RUR_AE_SAR();
    ae_int16x4 * restrict px = (ae_int16x4 *)x;
    int i;

    sar = 1 & (sar >> 5);


    WUR_AE_SAR(sar);

    ae_int16x4  * restrict y10,
        *restrict y11,
        *restrict y12,
        *restrict y13;
    ae_int16x4 h2, h3, h4, h5;
    ae_int16x4 e0, e1, e2, e3;
    int ibitrev_off = 0x80000000 >> (29 - NSA(N));
    int ibitrev = 0;

    y10 = (ae_int16x4 *)(y);
    y11 = (ae_int16x4 *)(y + N);
    y12 = (ae_int16x4 *)(y + N / 2);
    y13 = (ae_int16x4 *)(y + 3 * N / 2);

    /*  13 cycles per pipeline stage in steady state with unroll=1
        9 cycles lower bound required by resources  */
    __Pragma("loop_count min=2");
    for (i = 0; i<(N / 4); i += 2)
    {
        ae_int16x4 tmp;
        e0 = px[i];
        e1 = px[i + 1];
        e2 = px[i + N / 4];
        e3 = px[i + 1 + N / 4];

        e0 = AE_SRAA16RS(e0, sar);
        e1 = AE_SRAA16RS(e1, sar);
        e2 = AE_SRAA16RS(e2, sar);
        e3 = AE_SRAA16RS(e3, sar);

        h5 = AE_SEL16_5410(e0, e2);
        h4 = AE_SEL16_7632(e0, e2);
        /*========================= Hifi3 ====================================*/
        tmp = AE_ADD16S(h4, h5); //AE_ADDANDSUBRNG16RAS_S1(h4, h5);
        h5 = AE_SUB16S(h4, h5);
        h4 = tmp;

        h3 = AE_SEL16_5410(e1, e3);
        h2 = AE_SEL16_7632(e1, e3);

        tmp = AE_ADD16S(h2, h3);    //AE_ADDANDSUBRNG16RAS_S1(h2, h3);
        h3 = AE_SUB16S(h2, h3);
        h2 = tmp;

        /* Store the 2 bfly outputs together */
        AE_S16X4_X(h4, y10, ibitrev);
        AE_S16X4_X(h5, y11, ibitrev);
        AE_S16X4_X(h2, y12, ibitrev);
        AE_S16X4_X(h3, y13, ibitrev);

        ibitrev = AE_ADDBRBA32(ibitrev, ibitrev_off);
    }
    return sar;
} //stage_last_DFT2xI2_inpl_ds

#endif //#if (XCHAL_HAVE_HIFI3Z)

/*-------------------------------------------------------------------------
  FFT on Complex Data
  These functions make FFT on complex data.
    Scaling: 
      +-------------------+----------------------------------------+
      |      Function     |           Scaling options              |
      +-------------------+----------------------------------------+
      |  fft_cplx16x16    |  2 - 16-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  fft_cplx32x32    |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  fft_cplx32x16    |  3 - fixed scaling before each stage   | 
      |  fft_cplx24x24    |  0 - no scaling                        | 
      |                   |  1 - 24-bit scaling                    |
      |                   |  2 - 32-bit scaling on the first stage |
      |                   |  and 24-bit scaling later              |
      |                   |  3 - fixed scaling before each stage   |
      +-------------------+----------------------------------------+
  NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED after 
     the call
  3. 32x32 FFTs support mixed radix transforms 
  4. N - FFT size

  Precision: 
  24x24  24-bit input/outputs, 24-bit twiddles
  32x16  32-bit input/outputs, 16-bit twiddles
  32x32  32-bit input/outputs, 32-bit twiddles
  16x16  16-bit input/outputs, 16-bit twiddles
 
  Input:
  x[2*N]     complex input signal. Real and imaginary data are interleaved 
             and real data goes first
  scalingOpt scaling option (see table above)
  Output:
  y[2*N]     output spectrum. Real and imaginary data are interleaved and 
             real data goes first

  Returned value: total number of right shifts occurred during scaling 
                  procedure

  Restrictions:
  x,y        should not overlap
  x,y        aligned on a 8-bytes boundary

-------------------------------------------------------------------------*/
int fft_cplx16x16( int16_t* y,int16_t* x,fft_handle_t h,int scalingOption)
#if !(XCHAL_HAVE_HIFI3Z)
{
    int scale;
    ae_int16x4 * restrict px;

    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT(scalingOption == 3 || scalingOption == 2);

    if (scalingOption == 3)
    {
        scale=FIRST_STAGE_SCALE;
        scale+=fft_core_s3 (y, x, (const tFftDescr*)h);
        return scale;
    }
    else
    {
        int N;
        const tFftDescr *pDescr = (const tFftDescr*)h;
        int stride;
        cint32_ptr *p = (cint32_ptr *)pDescr->twd;
        int16_t *p_twd = (int16_t*)*p++;
        
        N = pDescr->N;
        stride = N; // The stride is quartered with every iteration of the outer loop.

        int shift;
        int bexp;        
        {
            int i;
            ae_int16x4 acc = AE_MOVINT16X4_FROMINT32X2(AE_MOVI(0)), tmp;
            px = (ae_int16x4*)x;

            /* 1 cycles per pipeline stage in steady state with unroll=1*/
            __Pragma("loop_count min=4 factor=4");
            for (i = 0; i < (N >> 1); i++)   
            {
                AE_L16X4_IP(tmp, px, sizeof(*px));
                tmp = AE_ABS16S(tmp); 
                acc = AE_OR16(acc, tmp);
            }
            acc = AE_OR16(acc, AE_SEL16_5432(acc, acc));
            acc = AE_OR16(acc, AE_SHORTSWAP(acc));

            i = AE_MOVAD16_0(acc);
            bexp = NSA(i) - 16;
            XT_MOVEQZ(bexp, 0, i);
        }

        shift = stage_first_DFT4xI2_inpl_ds(x, (int16_t*)p_twd, N, bexp);
        stride >>= 2;

        while (stride > 4)
        {
            p_twd = (int16_t*)*p++;
            shift += stage_inner_DFT4xI2_inpl_ref_ds(x, x, p_twd, N, stride);
            stride >>= 2;
        }
        if (stride == 4)
        {
            shift += stage_last_DFT4xI2_inpl_ds(x, y, N);
        }
        else
        {
            shift += stage_last_DFT2xI2_inpl_ds(x, y, N);
        }

        return  shift;
      }
    
}
#else
{
   
    int N;
    const tFftDescr *pDescr = (const tFftDescr*)h;
    int stride;
    cint32_ptr *p = (cint32_ptr *)pDescr->twd;
    ae_int16x4 * restrict px; 
    int16_t *p_twd = (int16_t*)*p++;

    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT(scalingOption == 3 || scalingOption == 2);
    N = pDescr->N;
    stride = N; // The stride is quartered with every iteration of the outer loop.
   // AE_CALCRNG3()
    if (scalingOption == 3)
    {

        stage_first_DFT4xI2_inpl(x, (int16_t*)p_twd, N);
        stride >>= 2;
        if (stride > 4)
        {
            p_twd = (int16_t*)*p++;


            stage_inner_DFT4xI2_inpl_merged(x, x, p_twd, N, stride);
            stride >>= 2;
            while (stride > 4)
            {
                p_twd = (int16_t*)*p++;
                stage_inner_DFT4xI2_inpl_ref(x, x, p_twd, N, stride);
                stride >>= 2;
            }
        }
        if (stride == 4)
        {
            stage_last_DFT4xI2_inpl(x, y, N);
        }
        else
        {
            stage_last_DFT2xI2_inpl(x, y, N);
        }
        return  31 - NSA(N);
    }
    else // if (scalingOption == 3)
    {
        int shift; 
        int bexp; 
        
        {
            int i; 
            ae_int16x4 acc = AE_MOVINT16X4_FROMINT32X2( AE_MOVI(0) ), tmp; 

            __Pragma("loop_count min=4 factor=4"); 
            px = (ae_int16x4*)x;
            for (i = 0; i < (N >> 1); i++)
            {
                AE_L16X4_IP(tmp, px, sizeof(*px) ); 
                acc = AE_MAXABS16S(acc, tmp); 
            }
            acc = AE_MAX16(acc, AE_SEL16_5432(acc, acc));
            acc = AE_MAX16(acc, AE_SHORTSWAP(acc) );
           
            i = AE_MOVAD16_0(acc);
            bexp = NSA(i)-16;
            XT_MOVEQZ(bexp, 0, i); 
        }

        ASSERT(bexp >= 0); 
       
        shift = stage_first_DFT4xI2_inpl_ds(x, (int16_t*)p_twd, N, bexp);
        stride >>= 2;
        //todo stage_second need here
        if (stride > 4)
        {
            p_twd = (int16_t*)*p++;
            shift += stage_second_DFT4xI2_inpl_ref_ds(x, x, p_twd, N, stride);
            stride >>= 2;
            while (stride > 4)
            {
                p_twd = (int16_t*)*p++;
                shift += stage_inner_DFT4xI2_inpl_ref_ds(x, x, p_twd, N, stride);
                stride >>= 2;
            }
        }
        if (stride == 4)
        {
            shift += stage_last_DFT4xI2_inpl_ds(x, y, N);
        }
        else
        {
           shift += stage_last_DFT2xI2_inpl_ds(x, y, N);
        }

       
        return  shift;
    } //if (scalingOption == 3) else
}
#endif


