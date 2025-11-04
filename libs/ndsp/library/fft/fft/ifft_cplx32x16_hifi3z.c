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
  C code optimized for HiFi3
	Integrit, 2006-2017
*/
#include "NatureDSP_Signal_fft.h"
#include "fft_cplx_twiddles.h"
#include "common.h"

#if ( XCHAL_HAVE_HIFI3Z)
#define FIRST_STAGE_SCALE 3

#define IFFT_BUTTERFLY_S3(_T,_step,_X,twd_step)\
      {                                         \
      vT1 = AE_L16X4_I(p_twd, 8);               \
      vB0 = AE_ADD32S(vA0, vA2);                \
      vB3 = AE_SUB32S(vA1, vA3);                \
      vT2 = AE_L16X4_I(p_twd, 16);              \
      vB2 = AE_SUB32S(vA0, vA2);                \
      vB1 = AE_ADD32S(vA1, vA3);                \
      AE_L16X4_##_X(vTL, p_twd, twd_step); \
      vT0 = (vTL);                              \
      vB3 = AE_SEL32_LH(vB3, vB3);              \
      vC1 = AE_SUB32S(vB0, vB1);                \
      vA0 = AE_L32X2_I(p_x0, 8);                \
      vC0 = AE_ADD32S(vB0, vB1);                \
      vC2 = AE_ADDSUB32S(vB2, vB3);             \
      vA1 = AE_L32X2_I(p_x1, 8);                \
      vC1 = AE_MULFC32X16RAS_L(vC1, vT0);       \
      vC3 = AE_SUBADD32S(vB2, vB3);             \
      vA2 = AE_L32X2_I(p_x2, 8);                \
      vC2 = AE_MULFC32X16RAS_H(vC2, vT1);       \
      vC0 = AE_SRAS32(vC0);                     \
      vA3 = AE_L32X2_I(p_x3, 8);                \
      vC3 = AE_MULFC32X16RAS_H(vC3, vT0);       \
      vC1 = AE_SRAS32(vC1);                     \
      AE_S32X2_IP(vC0, p_x0, 8);                \
      vB0 = AE_ADD32S(vA0, vA2);                \
      vC2 = AE_SRAS32(vC2);                     \
      AE_S32X2_IP(vC1, p_x1, 8);                \
      vB2 = AE_SUB32S(vA0, vA2);                \
      vC3 = AE_SRAS32(vC3);                     \
      AE_S32X2_IP(vC2, p_x3, 8);                \
      vB1 = AE_ADD32S(vA1, vA3);                \
      vB3 = AE_SUB32S(vA1, vA3);                \
      AE_S32X2_IP(vC3, p_x2, 8);                \
      vC1 = AE_SUB32S(vB0, vB1);                \
      vB3 = AE_SEL32_LH(vB3, vB3);              \
      vA0 = AE_L32X2_##_T(p_x0, _step);         \
      vC0 = AE_ADD32S(vB0, vB1);                \
      vC2 = AE_ADDSUB32S(vB2, vB3);             \
      vA1 = AE_L32X2_##_T(p_x1, _step);         \
      vC1 = AE_MULFC32X16RAS_H(vC1, vT2);       \
      vC3 = AE_SUBADD32S(vB2, vB3);             \
      vA2 = AE_L32X2_##_T(p_x2, _step);         \
      vC2 = AE_MULFC32X16RAS_L(vC2, vT2);       \
      vC0 = AE_SRAS32(vC0);                     \
      vA3 = AE_L32X2_##_T(p_x3, _step);         \
      vC3 = AE_MULFC32X16RAS_L(vC3, vT1);       \
      vC1 = AE_SRAS32(vC1);                     \
      AE_S32X2_##_T##P(vC0, p_x0, _step);       \
      vC2 = AE_SRAS32(vC2);                     \
      AE_S32X2_##_T##P(vC1, p_x1, _step);       \
      vC3 = AE_SRAS32(vC3);                     \
      AE_S32X2_##_T##P(vC2, p_x3, _step);       \
      AE_S32X2_##_T##P(vC3, p_x2, _step);       \
      }

#define IFFT_BUTTERFLY_S3_LAST(_T,_step,_X,twd_step)\
      {                                         \
      vT1 = AE_L16X4_I(p_twd, 8);               \
      vB0 = AE_ADD32S(vA0, vA2);                \
      vB3 = AE_SUB32S(vA1, vA3);                \
      vT2 = AE_L16X4_I(p_twd, 16);              \
      vB2 = AE_SUB32S(vA0, vA2);                \
      vB1 = AE_ADD32S(vA1, vA3);                \
      AE_L16X4_##_X(vTL, p_twd, twd_step);      \
      vT0 = (vTL);                              \
      vB3 = AE_SEL32_LH(vB3, vB3);              \
      vC1 = AE_SUB32S(vB0, vB1);                \
      vA0 = AE_L32X2_I(p_x0, 8);                \
      vC0 = AE_ADD32S(vB0, vB1);                \
      vC2 = AE_ADDSUB32S(vB2, vB3);             \
      vA1 = AE_L32X2_I(p_x1, 8);                \
      vC1 = AE_MULFC32X16RAS_L(vC1, vT0);       \
      vC3 = AE_SUBADD32S(vB2, vB3);             \
      vA2 = AE_L32X2_I(p_x2, 8);                \
      vC2 = AE_MULFC32X16RAS_H(vC2, vT1);       \
      vC0 = AE_SRAS32(vC0);                     \
      vA3 = AE_L32X2_I(p_x3, 8);                \
      vC3 = AE_MULFC32X16RAS_H(vC3, vT0);       \
      vC1 = AE_SRAS32(vC1);                     \
      AE_S32X2_IP(vC0, p_x0, 8);                \
      vB0 = AE_ADD32S(vA0, vA2);                \
      vC2 = AE_SRAS32(vC2);                     \
      AE_S32X2_IP(vC1, p_x1, 8);                \
      vB2 = AE_SUB32S(vA0, vA2);                \
      vC3 = AE_SRAS32(vC3);                     \
      AE_S32X2_IP(vC2, p_x3, 8);                \
      vB1 = AE_ADD32S(vA1, vA3);                \
      vB3 = AE_SUB32S(vA1, vA3);                \
      AE_S32X2_IP(vC3, p_x2, 8);                \
      vC1 = AE_SUB32S(vB0, vB1);                \
      vB3 = AE_SEL32_LH(vB3, vB3);              \
      vC0 = AE_ADD32S(vB0, vB1);                \
      vC2 = AE_ADDSUB32S(vB2, vB3);             \
      vC1 = AE_MULFC32X16RAS_H(vC1, vT2);       \
      vC3 = AE_SUBADD32S(vB2, vB3);             \
      vC2 = AE_MULFC32X16RAS_L(vC2, vT2);       \
      vC0 = AE_SRAS32(vC0);                     \
      vC3 = AE_MULFC32X16RAS_L(vC3, vT1);       \
      vC1 = AE_SRAS32(vC1);                     \
      AE_S32X2_##_T##P(vC0, p_x0, _step);       \
      vC2 = AE_SRAS32(vC2);                     \
      AE_S32X2_##_T##P(vC1, p_x1, _step);       \
      vC3 = AE_SRAS32(vC3);                     \
      AE_S32X2_##_T##P(vC2, p_x3, _step);       \
      AE_S32X2_##_T##P(vC3, p_x2, _step);       \
      }

#define IFFT_BUTTERFLY_S3_T3F(_T, _step)        \
      {                                         \
      vA0 = AE_L32X2_I(p_x0, 8);                \
      vC0 = AE_ADD32S(vB0, vB1);                \
      vC2 = AE_ADDSUB32S(vB2, vB3);             \
      vA1 = AE_L32X2_I(p_x1, 8);                \
      vC1 = AE_MULFC32X16RAS_L(vC1, vT0);       \
      vC3 = AE_SUBADD32S(vB2, vB3);             \
      vA2 = AE_L32X2_I(p_x2, 8);                \
      vC2 = AE_MULFC32X16RAS_H(vC2, vT1);       \
      vC0 = AE_SRAS32(vC0);                     \
      vA3 = AE_L32X2_I(p_x3, 8);                \
      vC3 = AE_MULFC32X16RAS_H(vC3, vT0);       \
      vC1 = AE_SRAS32(vC1);                     \
      AE_S32X2_IP(vC0, p_x0, 8);                \
      vB0 = AE_ADD32S(vA0, vA2);                \
      vC2 = AE_SRAS32(vC2);                     \
      AE_S32X2_IP(vC1, p_x1, 8);                \
      vB2 = AE_SUB32S(vA0, vA2);                \
      vC3 = AE_SRAS32(vC3);                     \
      AE_S32X2_IP(vC2, p_x3, 8);                \
      vB1 = AE_ADD32S(vA1, vA3);                \
      vB3 = AE_SUB32S(vA1, vA3);                \
      AE_S32X2_IP(vC3, p_x2, 8);                \
      vC1 = AE_SUB32S(vB0, vB1);                \
      vB3 = AE_SEL32_LH(vB3, vB3);              \
      vA0 = AE_L32X2_##_T(p_x0, _step);         \
      vC0 = AE_ADD32S(vB0, vB1);                \
      vC2 = AE_ADDSUB32S(vB2, vB3);             \
      vA1 = AE_L32X2_##_T(p_x1, _step);         \
      vC1 = AE_MULFC32X16RAS_H(vC1, vT2);       \
      vC3 = AE_SUBADD32S(vB2, vB3);             \
      vA2 = AE_L32X2_##_T(p_x2, _step);         \
      vC2 = AE_MULFC32X16RAS_L(vC2, vT2);       \
      vC0 = AE_SRAS32(vC0);                     \
      vA3 = AE_L32X2_##_T(p_x3, _step);         \
      vC3 = AE_MULFC32X16RAS_L(vC3, vT1);       \
      vC1 = AE_SRAS32(vC1);                     \
      AE_S32X2_##_T##P(vC0, p_x0, _step);       \
      vB0 = AE_ADD32S(vA0, vA2);                \
      vC2 = AE_SRAS32(vC2);                     \
      AE_S32X2_##_T##P(vC1, p_x1, _step);       \
      vB3 = AE_SUB32S(vA1, vA3);                \
      vC3 = AE_SRAS32(vC3);                     \
      AE_S32X2_##_T##P(vC2, p_x3, _step);       \
      vB2 = AE_SUB32S(vA0, vA2);                \
      vB1 = AE_ADD32S(vA1, vA3);                \
      AE_S32X2_##_T##P(vC3, p_x2, _step);       \
      vB3 = AE_SEL32_LH(vB3, vB3);              \
      vC1 = AE_SUB32S(vB0, vB1);                \
      }

#define IFFT_BUTTERFLY_S3_T3F_LAST(_T, _step)   \
      {                                         \
      vA0 = AE_L32X2_I(p_x0, 8);                \
      vC0 = AE_ADD32S(vB0, vB1);                \
      vC2 = AE_ADDSUB32S(vB2, vB3);             \
      vA1 = AE_L32X2_I(p_x1, 8);                \
      vC1 = AE_MULFC32X16RAS_L(vC1, vT0);       \
      vC3 = AE_SUBADD32S(vB2, vB3);             \
      vA2 = AE_L32X2_I(p_x2, 8);                \
      vC2 = AE_MULFC32X16RAS_H(vC2, vT1);       \
      vC0 = AE_SRAS32(vC0);                     \
      vA3 = AE_L32X2_I(p_x3, 8);                \
      vC3 = AE_MULFC32X16RAS_H(vC3, vT0);       \
      vC1 = AE_SRAS32(vC1);                     \
      AE_S32X2_IP(vC0, p_x0, 8);                \
      vB0 = AE_ADD32S(vA0, vA2);                \
      vC2 = AE_SRAS32(vC2);                     \
      AE_S32X2_IP(vC1, p_x1, 8);                \
      vB2 = AE_SUB32S(vA0, vA2);                \
      vC3 = AE_SRAS32(vC3);                     \
      AE_S32X2_IP(vC2, p_x3, 8);                \
      vB1 = AE_ADD32S(vA1, vA3);                \
      vB3 = AE_SUB32S(vA1, vA3);                \
      AE_S32X2_IP(vC3, p_x2, 8);                \
      vC1 = AE_SUB32S(vB0, vB1);                \
      vB3 = AE_SEL32_LH(vB3, vB3);              \
      vC0 = AE_ADD32S(vB0, vB1);                \
      vC2 = AE_ADDSUB32S(vB2, vB3);             \
      vC1 = AE_MULFC32X16RAS_H(vC1, vT2);       \
      vC3 = AE_SUBADD32S(vB2, vB3);             \
      vC2 = AE_MULFC32X16RAS_L(vC2, vT2);       \
      vC0 = AE_SRAS32(vC0);                     \
      vC3 = AE_MULFC32X16RAS_L(vC3, vT1);       \
      vC1 = AE_SRAS32(vC1);                     \
      AE_S32X2_##_T##P(vC0, p_x0, _step);       \
      vC2 = AE_SRAS32(vC2);                     \
      AE_S32X2_##_T##P(vC1, p_x1, _step);       \
      vC3 = AE_SRAS32(vC3);                     \
      AE_S32X2_##_T##P(vC2, p_x3, _step);       \
      AE_S32X2_##_T##P(vC3, p_x2, _step);       \
      }

#define IFFT_BUTTERFLY_R2(idx)      \
    {                               \
      vA1 = AE_L32X2_I(p_x0, 8);    \
      vA2 = AE_L32X2_I(p_x0, 16);   \
      vA3 = AE_L32X2_I(p_x0, 24);   \
      AE_L32X2_XP(vA0, p_x0, 32);   \
      vB0 = AE_ADD32S(vA0, vA1);    \
      vB2 = AE_SUB32S(vA0, vA1);    \
      vB1 = AE_ADD32S(vA2, vA3);    \
      vB3 = AE_SUB32S(vA2, vA3);    \
      AE_S32X2_X(vB0, p_y0, idx);   \
      AE_S32X2_X(vB1, p_y1, idx);   \
      AE_S32X2_X(vB2, p_y2, idx);   \
      AE_S32X2_X(vB3, p_y3, idx);   \
    }

#define IFFT_BUTTERFLY_R4(idx)      \
    {                               \
      vA1 = AE_L32X2_I(p_x0, 8);    \
      vA2 = AE_L32X2_I(p_x0, 16);   \
      vA3 = AE_L32X2_I(p_x0, 24);   \
      AE_L32X2_XP(vA0, p_x0, 32);   \
      vB0 = AE_ADD32S(vA0, vA2);    \
      vB2 = AE_SUB32S(vA0, vA2);    \
      vB1 = AE_ADD32S(vA1, vA3);    \
      vB3 = AE_SUB32S(vA1, vA3);    \
      vB3 = AE_SEL32_LH(vB3, vB3);  \
      vA0 = AE_ADD32S(vB0, vB1);    \
      vA2 = AE_SUB32S(vB0, vB1);    \
      vA1 = AE_SUBADD32S(vB2, vB3); \
      vA3 = AE_ADDSUB32S(vB2, vB3); \
      AE_S32X2_X(vA0, p_y0, idx);   \
      AE_S32X2_X(vA1, p_y1, idx);   \
      AE_S32X2_X(vA2, p_y2, idx);   \
      AE_S32X2_X(vA3, p_y3, idx);   \
    }

#if !(XCHAL_HAVE_HIFI3Z)
static ALIGN(16) const int16_t  ifft_twd_r8[12]  = {
    (int16_t)0x7fff, (int16_t)0x0000, (int16_t)0x7fff, (int16_t)0x0000, 
    (int16_t)0x7fff, (int16_t)0x0000, (int16_t)0x5a82, (int16_t)0x5a82, 
    (int16_t)0x0000, (int16_t)0x7fff, (int16_t)0xa57e, (int16_t)0x5a82,
};

static int ifft_core_s3(int32_t *y, int32_t *x, const tFftDescr* pDescr)
{
  ae_int32x2 * restrict p_x0, * restrict p_x1,
             * restrict p_x2, * restrict p_x3;
  const ae_int16x4 * restrict p_twd;
  int step_circ, lc;

  ae_int32x2  vA0, vA1, vA2, vA3, vB0, vB1, vB2, vB3, vC0, vC1, vC2, vC3;
  ae_int16x4  vTL;
  ae_f16x4    vT0, vT1, vT2;
  const ae_int16 * restrict p_inc;
  int stride, scale = 0;
  int rsa;
  int N;
  const cint32_ptr * restrict seq;

  {
  	N=pDescr->N;
  	p_inc = (const ae_int16 *)(pDescr->inc);
    seq=(const cint32_ptr * )pDescr->twd;
    p_twd = (const ae_int16x4 *)*seq++;
  }

  stride     =   N; // The stride is quartered with every iteration of the outer loop.
  rsa = 2;

  if (stride > 4)
  {
    //-----------------------------------------------------------------
    // Set up pointers to access "N/4", "N/2", "3N/4" complex point or  
    // "N/2", "N", "3N/2" half word                                    
    //-----------------------------------------------------------------
    p_x0 = (ae_int32x2 *)x;
    p_x1 = p_x0 + (stride>>2);
    p_x2 = p_x1 + (stride>>2);
    p_x3 = p_x2 + (stride>>2);

    vA0 = AE_L32X2_I(p_x0, 0);
    vA1 = AE_L32X2_I(p_x1, 0);
    vA2 = AE_L32X2_I(p_x2, 0);
    vA3 = AE_L32X2_I(p_x3, 0);

    vA0 = AE_SRAI32(vA0, FIRST_STAGE_SCALE);
    vA1 = AE_SRAI32(vA1, FIRST_STAGE_SCALE);
    vA2 = AE_SRAI32(vA2, FIRST_STAGE_SCALE);
    vA3 = AE_SRAI32(vA3, FIRST_STAGE_SCALE);

    if (stride==8) rsa >>= 1;

    WUR_AE_SAR(rsa);

    vB0 = AE_ADD32S(vA0, vA2);
    vB2 = AE_SUB32S(vA0, vA2);
    vB1 = AE_ADD32S(vA1, vA3);
    vB3 = AE_SUB32S(vA1, vA3);

    vT1 = AE_L16X4_I(p_twd, 8);
    vB3 = AE_SEL32_LH(vB3, vB3);
    vC0 = AE_ADD32S(vB0, vB1);

    for(lc=0; lc<((N>>3)-1); lc++)
    {
      AE_L16X4_IP(vTL, p_twd, 24);
      vT0 = (vTL);
      vC1 = AE_SUB32S(vB0, vB1);
      vC2 = AE_ADDSUB32S(vB2, vB3);
      vA0 = AE_L32X2_I(p_x0, 8);
      vC3 = AE_SUBADD32S(vB2, vB3);
      vC0 = AE_SRAS32(vC0);
      vA2 = AE_L32X2_I(p_x2, 8);
      vC1 = AE_MULFC32X16RAS_L(vC1, vT0);
      vA0 = AE_SRAI32(vA0, FIRST_STAGE_SCALE);
      vA1 = AE_L32X2_I(p_x1, 8);
      vC2 = AE_MULFC32X16RAS_H(vC2, vT1);
      vA2 = AE_SRAI32(vA2, FIRST_STAGE_SCALE);
      vA3 = AE_L32X2_I(p_x3, 8);
      vC3 = AE_MULFC32X16RAS_H(vC3, vT0);
      vA1 = AE_SRAI32(vA1, FIRST_STAGE_SCALE);
      vT2 = AE_L16X4_I(p_twd, -8);
      vB0 = AE_ADD32S(vA0, vA2);
      vA3 = AE_SRAI32(vA3, FIRST_STAGE_SCALE);
      AE_S32X2_IP(vC0, p_x0, 8);
      vB2 = AE_SUB32S(vA0, vA2);
      vC1 = AE_SRAS32(vC1);
      AE_S32X2_IP(vC1, p_x1, 8);
      vB1 = AE_ADD32S(vA1, vA3);
      vC2 = AE_SRAS32(vC2);
      AE_S32X2_IP(vC2, p_x3, 8);
      vB3 = AE_SUB32S(vA1, vA3);
      vC3 = AE_SRAS32(vC3);
      AE_S32X2_IP(vC3, p_x2, 8);
      vB3 = AE_SEL32_LH(vB3, vB3);
      vC0 = AE_ADD32S(vB0, vB1);
      vA0 = AE_L32X2_I(p_x0, 8);
      vC1 = AE_SUB32S(vB0, vB1);
      vC0 = AE_SRAS32(vC0);
      vA1 = AE_L32X2_I(p_x1, 8);
      vC1 = AE_MULFC32X16RAS_H(vC1, vT2);
      vC2 = AE_ADDSUB32S(vB2, vB3);
      vA0 = AE_SRAI32(vA0, FIRST_STAGE_SCALE);
      vA2 = AE_L32X2_I(p_x2, 8);
      vC2 = AE_MULFC32X16RAS_L(vC2, vT2);
      vC3 = AE_SUBADD32S(vB2, vB3);
      vA1 = AE_SRAI32(vA1, FIRST_STAGE_SCALE);
      vA3 = AE_L32X2_I(p_x3, 8);
      vC3 = AE_MULFC32X16RAS_L(vC3, vT1);
      vA2 = AE_SRAI32(vA2, FIRST_STAGE_SCALE);
      AE_S32X2_IP(vC0, p_x0, 8);
      vB0 = AE_ADD32S(vA0, vA2);
      vC1 = AE_SRAS32(vC1);
      AE_S32X2_IP(vC1, p_x1, 8);
      vB2 = AE_SUB32S(vA0, vA2);
      vA3 = AE_SRAI32(vA3, FIRST_STAGE_SCALE);
      vT1 = AE_L16X4_I(p_twd, 8);
      vB3 = AE_SUB32S(vA1, vA3);
      vC2 = AE_SRAS32(vC2);
      AE_S32X2_IP(vC2, p_x3, 8);
      vB1 = AE_ADD32S(vA1, vA3);
      vC3 = AE_SRAS32(vC3);
      AE_S32X2_IP(vC3, p_x2, 8);
      vB3 = AE_SEL32_LH(vB3, vB3);
      vC0 = AE_ADD32S(vB0, vB1);
    };
    {
      AE_L16X4_IP(vTL, p_twd, 24);
      vT0 = (vTL);
      vC1 = AE_SUB32S(vB0, vB1);
      vC2 = AE_ADDSUB32S(vB2, vB3);
      vA0 = AE_L32X2_I(p_x0, 8);
      vC3 = AE_SUBADD32S(vB2, vB3);
      vC0 = AE_SRAS32(vC0);
      vA2 = AE_L32X2_I(p_x2, 8);
      vC1 = AE_MULFC32X16RAS_L(vC1, vT0);
      vA0 = AE_SRAI32(vA0, FIRST_STAGE_SCALE);
      vA1 = AE_L32X2_I(p_x1, 8);
      vC2 = AE_MULFC32X16RAS_H(vC2, vT1);
      vA2 = AE_SRAI32(vA2, FIRST_STAGE_SCALE);
      vA3 = AE_L32X2_I(p_x3, 8);
      vC3 = AE_MULFC32X16RAS_H(vC3, vT0);
      vA1 = AE_SRAI32(vA1, FIRST_STAGE_SCALE);
      vT2 = AE_L16X4_I(p_twd, -8);
      vB0 = AE_ADD32S(vA0, vA2);
      vA3 = AE_SRAI32(vA3, FIRST_STAGE_SCALE);
      AE_S32X2_IP(vC0, p_x0, 8);
      vB2 = AE_SUB32S(vA0, vA2);
      vC1 = AE_SRAS32(vC1);
      AE_S32X2_IP(vC1, p_x1, 8);
      vB1 = AE_ADD32S(vA1, vA3);
      vC2 = AE_SRAS32(vC2);
      AE_S32X2_IP(vC2, p_x3, 8);
      vB3 = AE_SUB32S(vA1, vA3);
      vC3 = AE_SRAS32(vC3);
      AE_S32X2_IP(vC3, p_x2, 8);
      vB3 = AE_SEL32_LH(vB3, vB3);
      vC0 = AE_ADD32S(vB0, vB1);
      vC1 = AE_SUB32S(vB0, vB1);
      vC0 = AE_SRAS32(vC0);
      vC1 = AE_MULFC32X16RAS_H(vC1, vT2);
      vC2 = AE_ADDSUB32S(vB2, vB3);
      vC2 = AE_MULFC32X16RAS_L(vC2, vT2);
      vC3 = AE_SUBADD32S(vB2, vB3);
      vC3 = AE_MULFC32X16RAS_L(vC3, vT1);
      AE_S32X2_IP(vC0, p_x0, 8);
      vC1 = AE_SRAS32(vC1);
      AE_S32X2_IP(vC1, p_x1, 8);
      vC2 = AE_SRAS32(vC2);
      AE_S32X2_IP(vC2, p_x3, 8);
      vC3 = AE_SRAS32(vC3);
      AE_S32X2_IP(vC3, p_x2, 8);
    }
    scale+= (stride!=8) ? 2:1;
    stride  >>= 2;
  }

  while (stride > 16)
  {
    int step;

    //-----------------------------------------------------------------
    // Set up pointers to access "N/4", "N/2", "3N/4" complex point or  
    // "N/2", "N", "3N/2" half word                                    
    //-----------------------------------------------------------------
    p_x0 = (ae_int32x2 *)x;
    p_x1 = p_x0 + (stride>>2);
    p_x2 = p_x1 + (stride>>2);
    p_x3 = p_x2 + (stride>>2);

    vA0 = AE_L32X2_I(p_x0, 0);
    vA1 = AE_L32X2_I(p_x1, 0);
    vA2 = AE_L32X2_I(p_x2, 0);
    vA3 = AE_L32X2_I(p_x3, 0);

    p_twd=(const ae_int16x4 *)(*seq++);
    WUR_AE_CBEGIN0( (uintptr_t)p_twd );
    WUR_AE_CEND0  ( (uintptr_t)(p_twd)+stride*3 );
    step_circ=24;

    __Pragma("loop_count min=1");
    for(lc=0; lc<(N>>4); lc++)
    {
      AE_L16_IP(vTL, p_inc, 2);
      vTL = AE_SLAI16S(vTL, 2);
      step = AE_MOVAD16_0(vTL);

      IFFT_BUTTERFLY_S3(I, 8,XC,step_circ);
      IFFT_BUTTERFLY_S3(X, step,XC,step_circ)
    };

    scale+= (stride!=8) ? 2:1;
    stride  >>= 2;
  }

  if (stride==16)
  {
    const int step_fwd = 26*4;
    int step_circ;

    //-----------------------------------------------------------------
    // Set up pointers to access "N/4", "N/2", "3N/4" complex point or  
    // "N/2", "N", "3N/2" half word                                    
    //-----------------------------------------------------------------
    p_x0 = (ae_int32x2 *)x;
    p_x1 = p_x0 + (stride>>2);
    p_x2 = p_x1 + (stride>>2);
    p_x3 = p_x2 + (stride>>2);

    vA0 = AE_L32X2_I(p_x0, 0);
    vA1 = AE_L32X2_I(p_x1, 0);
    vA2 = AE_L32X2_I(p_x2, 0);
    vA3 = AE_L32X2_I(p_x3, 0);
    p_twd=(const ae_int16x4 *)(*seq++);
    WUR_AE_CBEGIN0( (uintptr_t)p_twd );
    WUR_AE_CEND0  ( (uintptr_t)(p_twd)+4*12 );
    step_circ=24;

    for(lc=0; lc<((N>>4)-1); lc++)
    {
      IFFT_BUTTERFLY_S3(X, 8,XC,step_circ)
      IFFT_BUTTERFLY_S3(X, step_fwd,XC,step_circ)
    };
    {
      IFFT_BUTTERFLY_S3(X, 8,XC,step_circ)
      IFFT_BUTTERFLY_S3_LAST(X, step_fwd,XC,step_circ)
    }

    scale+= 2;
    stride  >>= 2;
  }


  if (stride==8) 
  {
    //-----------------------------------------------------------------
    // Set up pointers to access "N/4", "N/2", "3N/4" complex point or  
    // "N/2", "N", "3N/2" half word                                    
    //-----------------------------------------------------------------
    p_x0 = (ae_int32x2 *)x;
    p_x1 = p_x0 + (stride>>2);
    p_x2 = p_x1 + (stride>>2);
    p_x3 = p_x2 + (stride>>2);

    vA0 = AE_L32X2_I(p_x0, 0);
    vA1 = AE_L32X2_I(p_x1, 0);
    vA2 = AE_L32X2_I(p_x2, 0);
    vA3 = AE_L32X2_I(p_x3, 0);

    rsa >>= 1;
    WUR_AE_SAR(rsa);
    vT0 = AE_L16X4_I(((const ae_int16x4*)ifft_twd_r8), 0);
    vT1 = AE_L16X4_I(((const ae_int16x4*)ifft_twd_r8), 8);
    vT2 = AE_L16X4_I(((const ae_int16x4*)ifft_twd_r8), 16);

    vB0 = AE_ADD32S(vA0, vA2);
    vB3 = AE_SUB32S(vA1, vA3);
    vB2 = AE_SUB32S(vA0, vA2);
    vB1 = AE_ADD32S(vA1, vA3);
    vB3 = AE_SEL32_LH(vB3, vB3);
    vC1 = AE_SUB32S(vB0, vB1);

    for(lc=0; lc<((N>>3)-1); lc++)
    {
      IFFT_BUTTERFLY_S3_T3F(I, 56)
    };
    IFFT_BUTTERFLY_S3_T3F_LAST(I, 56)

    scale += 1;
    stride >>= 2;
  }

  {
    int32_t i,i0,i1,i2,i3,ai;
    ae_int32x2 *p_y0 = (ae_int32x2 *)(y);
    ae_int32x2 *p_y1 = (p_y0 + (N >> 2));
    ae_int32x2 *p_y2 = (p_y1 + (N >> 2));
    ae_int32x2 *p_y3 = (p_y2 + (N >> 2));
    p_x0 = (ae_int32x2 *)(x);

    i = AE_NSAZ32_L(N)+1;
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
          IFFT_BUTTERFLY_R2(i0);
          IFFT_BUTTERFLY_R2(i1);
          IFFT_BUTTERFLY_R2(i2);
          IFFT_BUTTERFLY_R2(i3);
          i0 = AE_ADDBRBA32(i3, ai);
        }
    } 
    else 
    {
        //--------------------------------------------------------------------------
        // last stage is RADIX4 !!!
        //--------------------------------------------------------------------------
        for (i = 0; i < (N>>4); i++) 
        {
          i1 = AE_ADDBRBA32(i0, ai);
          i2 = AE_ADDBRBA32(i1, ai);
          i3 = AE_ADDBRBA32(i2, ai);
          //----------------------------------------------------------------
          // Read eight inputs, and perform radix4 decomposition            
          //----------------------------------------------------------------
          IFFT_BUTTERFLY_R4(i0);
          IFFT_BUTTERFLY_R4(i1);
          IFFT_BUTTERFLY_R4(i2);
          IFFT_BUTTERFLY_R4(i3);
          i0 = AE_ADDBRBA32(i3, ai);
        }
    }
  }


  return scale;
}
#endif //#if XCHAL_HAVE_HIFI3Z
/*-------------------------------------------------------------------------
  Inverse FFT on Complex Data
  These functions make inverse FFT on complex data.
      Scaling: 
      +-------------------+----------------------------------------+
      |      Function     |           Scaling options              |
      +-------------------+----------------------------------------+
      |  ifft_cplx16x16   |  2 - 16-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  ifft_cplx32x32   |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  ifft_cplx32x16   |  3 - fixed scaling before each stage   | 
      |  ift_cplx24x24    |  0 - no scaling                        | 
      |                   |  1 - 24-bit scaling                    |
      |                   |  2 - 32-bit scaling on the first stage |
      |                   |  and 24-bit scaling later              |
      |                   |  3 - fixed scaling before each stage   |
      +-------------------+----------------------------------------+
  NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED after 
     the call.
  3. 32x32 FFTs support mixed radix transforms.
  4. N - FFT size.

  Precision: 
  32x32  32-bit input/outputs, 32-bit twiddles
  24x24  24-bit input/outputs, 24-bit twiddles
  32x16  32-bit input/outputs, 16-bit twiddles
  16x16  16-bit input/outputs, 16-bit twiddles

  Input:
  x[2*N]     complex input spectrum. Real and imaginary data are interleaved 
             and real data goes first
  scalingOpt scaling option (see table above)
  Output:
  y[2*N]     complex output signal. Real and imaginary data are interleaved 
             and real data goes first

  Returned value: total number of right shifts occurred during scaling 
                  procedure

  Restrictions:
  x,y        should not overlap
  x,y        aligned on a 8-bytes boundary
-------------------------------------------------------------------------*/

extern void fft_revorder_ie(int32_t *x, int N);
extern int fft_stage_last_ie(int32_t *x,
    int32_t *y,
    int N);


#define SWAP_PINT32(_x, _y) \
{ int32_t *_tmp = _x; _x = _y; _y=_tmp; }

int ifft_cplx32x16( int32_t* y,int32_t* x,fft_handle_t h,int scalingOption)
{

    int N;
    const cint32_ptr * restrict seq;
    int shift = 0;

    //int32_t *pdst = y; 
    const ae_int32x2 * restrict px0;
    const ae_int32x2 * restrict px1;
          ae_int32x2 * restrict py0;
          ae_int32x2 * restrict py1;
    const ae_int32x2 * restrict p16tw1;

    {
        N = ((const tFftDescr*)h)->N;
        seq = (const cint32_ptr *)((const tFftDescr*)h)->twd;
    }

    NASSERT(scalingOption==3);
    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);

    if(N==16)
    {
        int i;
        px0 = (ae_int32x2*)x; 
        py0 = (ae_int32x2*)y; 
 
        // memcpy(y, x, sizeof(*x)*2*N); 
        for(i = 0; i < 16; i++)
        {
            ae_int32x2 tmp; 
            AE_L32X2_IP(tmp, px0, 8); 
            AE_S32X2_IP(tmp, py0, 8); 
        }        
        SWAP_PINT32(x, y); 
    }

    {
        NASSERT_ALIGN8(x);
        NASSERT_ALIGN8(y);
        int s;
        int stride = N/4;     


        ae_int32x2  vA0, vA1, vA2, vA3,
            vB0, vB1, vB2, vB3,
            vC0, vC1, vC2, vC3;
        int i;
        int log2n = 0;

        WUR_AE_CBEGIN0((unsigned)x);
        WUR_AE_CEND0((unsigned)(uintptr_t)(((complex_fract32*)x) + (N - 2)));

        s = 3;
        WUR_AE_SAR(s);
        { // Start of the first stage

            p16tw1 = (ae_int32x2*)*seq++;
            /*
            unsigned acc = 0;
            const unsigned acc_inc = (log2n == 0) ? 0 : (0x80000000 >> (log2n - 1));

            tw_step = 1;
            */
            //!!!!
            log2n += 2;
            shift += s;
            
            px0 = (ae_int32x2*)x; // First stage is not inplace !!!
            px1 = px0 + 1;
            py0 = (ae_int32x2*)y;
            py1 = py0 + 1;
#if 1
            {     // first butterfly r4
                ae_int32x2 t1, t2, t3;
                ae_f16x4 t1_f16x4;
                ae_f16x4 t2_f16x4;
                ae_f16x4 t3_f16x4;
                /*
                int offset_inc = 0;
                acc += acc_inc;
                XT_MOVEQZ(offset_inc, 3 * 8, acc);
                */
                //                t3 = AE_L32X2_I(p16tw1, 8 * 2);
                //              t2 = AE_L32X2_I(p16tw1, 8 * 1);
                AE_L32X2_XP(t1, p16tw1, 8);
                AE_L32X2_XP(t2, p16tw1, 8);
                AE_L32X2_XP(t3, p16tw1, 8);

                t1_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
                t2_f16x4 = AE_MOVF16X4_FROMINT32X2(t2);
                t3_f16x4 = AE_MOVF16X4_FROMINT32X2(t3);
                //////////////

                vA3 = AE_L32X2_X(px0, 8 * (N - (0 * 2 + 3 * stride)));
                vA2 = AE_L32X2_X(px0, 8 * (N - (0 * 2 + 2 * stride)));
                vA1 = AE_L32X2_X(px0, 8 * (N - (0 * 2 + 1 * stride)));
                vA0 = AE_L32X2_X(px0, 0);

                vA3 = AE_SRAS32(vA3);
                vA2 = AE_SRAS32(vA2);
                vA1 = AE_SRAS32(vA1);
                vA0 = AE_SRAS32(vA0);

                vB0 = AE_ADD32S(vA0, vA2);
                vB2 = AE_SUB32S(vA0, vA2);
                vB1 = AE_ADD32S(vA1, vA3);
                vB3 = AE_SUB32S(vA1, vA3);

                vB3 = AE_SEL32_LH(vB3, vB3);

                vC0 = AE_ADD32S(vB0, vB1);
                vC2 = AE_SUB32S(vB0, vB1);
                vC3 = AE_SUBADD32S(vB2, vB3);
                vC1 = AE_ADDSUB32S(vB2, vB3);

                vC1 = AE_MULFC32X16RAS_H(vC1, t1_f16x4);
                vC2 = AE_MULFC32X16RAS_H(vC2, t2_f16x4);
                vC3 = AE_MULFC32X16RAS_H(vC3, t3_f16x4);

                AE_S32X2_X(vC3, py0, 3 * 8 * stride);
                AE_S32X2_X(vC2, py0, 1 * 8 * stride);
                AE_S32X2_X(vC1, py0, 2 * 8 * stride);
                AE_S32X2_XP(vC0, py0, 2 * 8);

                vA3 = AE_L32X2_X(px0, 8 * (N - (0 * 2 + 1 + 3 * stride)));
                vA2 = AE_L32X2_X(px0, 8 * (N - (0 * 2 + 1 + 2 * stride)));
                vA1 = AE_L32X2_X(px0, 8 * (N - (0 * 2 + 1 + 1 * stride)));
                vA0 = AE_L32X2_X(px0, 8 * (N - (0 * 2 + 1 + 0)));

                vA3 = AE_SRAS32(vA3);
                vA2 = AE_SRAS32(vA2);
                vA1 = AE_SRAS32(vA1);
                vA0 = AE_SRAS32(vA0);

                vB0 = AE_ADD32S(vA0, vA2);
                vB2 = AE_SUB32S(vA0, vA2);
                vB1 = AE_ADD32S(vA1, vA3);
                vB3 = AE_SUB32S(vA1, vA3);

                vB3 = AE_SEL32_LH(vB3, vB3);

                vC0 = AE_ADD32S(vB0, vB1);
                vC2 = AE_SUB32S(vB0, vB1);
                vC3 = AE_SUBADD32S(vB2, vB3);
                vC1 = AE_ADDSUB32S(vB2, vB3);

                vC1 = AE_MULFC32X16RAS_L(vC1, t1_f16x4);
                vC2 = AE_MULFC32X16RAS_L(vC2, t2_f16x4);
                vC3 = AE_MULFC32X16RAS_L(vC3, t3_f16x4);

                AE_S32X2_X(vC3, py1, 3 * 8 * stride);
                AE_S32X2_X(vC2, py1, 1 * 8 * stride);
                AE_S32X2_X(vC1, py1, 2 * 8 * stride);
                AE_S32X2_XP(vC0, py1, 2 * 8);
            } //  end first butterfly r4

            /*
            17 cycles per pipeline stage in steady state with unroll=1
            8.5 cycles per DFT4
            */
            px0 = (ae_int32x2*)(8 * (N - (1 * 2) - 3*stride) + (uintptr_t)x); 
            px1 = (ae_int32x2*)(8 * (N - (1 * 2 + 1) - 3 * stride) + (uintptr_t)x);
#endif
            __Pragma("loop_count min=1");
            __Pragma("concurrent")
            for (i = 0; i < (N>>3)-1; i++)
            {
                ae_int32x2 t1, t2, t3;
                ae_f16x4 t1_f16x4;
                ae_f16x4 t2_f16x4;
                ae_f16x4 t3_f16x4;

                AE_L32X2_IP(t1, p16tw1, 8);
                AE_L32X2_IP(t2, p16tw1, 8);
                AE_L32X2_IP(t3, p16tw1, 8);

                t1_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
                t2_f16x4 = AE_MOVF16X4_FROMINT32X2(t2);
                t3_f16x4 = AE_MOVF16X4_FROMINT32X2(t3);

                AE_L32X2_XP(vA3, px0, 8 * stride);
                AE_L32X2_XP(vA2, px0, 8 * stride);
                AE_L32X2_XP(vA1, px0, 8 * stride );
                AE_L32X2_XP(vA0, px0,  (-3*stride-2)*8 );
//                vA3=AE_L32X2_X(px0, 1*8 * stride);
//                vA2=AE_L32X2_X(px0, 2*8 * stride);
//                vA1=AE_L32X2_X(px0, 3*8 * stride );
//                AE_L32X2_XP(vA0, px0,  -2*8 );

                vA3 = AE_SRAS32(vA3);
                vA2 = AE_SRAS32(vA2);
                vA1 = AE_SRAS32(vA1);
                vA0 = AE_SRAS32(vA0);

                vB0 = AE_ADD32S(vA0, vA2);
                vB2 = AE_SUB32S(vA0, vA2);
                vB1 = AE_ADD32S(vA1, vA3);
                vB3 = AE_SUB32S(vA1, vA3);

                vB3 = AE_SEL32_LH(vB3, vB3);

                vC0 = AE_ADD32S(vB0, vB1);
                vC2 = AE_SUB32S(vB0, vB1);
                vC3 = AE_SUBADD32S(vB2, vB3);
                vC1 = AE_ADDSUB32S(vB2, vB3);

                vC1 = AE_MULFC32X16RAS_H(vC1, t1_f16x4);
                vC2 = AE_MULFC32X16RAS_H(vC2, t2_f16x4);
                vC3 = AE_MULFC32X16RAS_H(vC3, t3_f16x4);

                AE_S32X2_X(vC3, py0, 3 * 8 * stride);
                AE_S32X2_X(vC2, py0, 1 * 8 * stride);
                AE_S32X2_X(vC1, py0, 2 * 8 * stride);
                AE_S32X2_IP(vC0, py0, 2 * 8 );

                AE_L32X2_XP(vA3, px1, 8 * stride);
                AE_L32X2_XP(vA2, px1, 8 * stride);
                AE_L32X2_XP(vA1, px1, 8 * stride);
                AE_L32X2_XP(vA0, px1, (-3 * stride - 2) * 8);

                vA3 = AE_SRAS32(vA3);
                vA2 = AE_SRAS32(vA2);
                vA1 = AE_SRAS32(vA1);
                vA0 = AE_SRAS32(vA0);

                vB0 = AE_ADD32S(vA0, vA2);
                vB2 = AE_SUB32S(vA0, vA2);
                vB1 = AE_ADD32S(vA1, vA3);
                vB3 = AE_SUB32S(vA1, vA3);

                vB3 = AE_SEL32_LH(vB3, vB3);

                vC0 = AE_ADD32S(vB0, vB1);
                vC2 = AE_SUB32S(vB0, vB1);
                vC3 = AE_SUBADD32S(vB2, vB3);
                vC1 = AE_ADDSUB32S(vB2, vB3);

                vC1 = AE_MULFC32X16RAS_L(vC1, t1_f16x4);
                vC2 = AE_MULFC32X16RAS_L(vC2, t2_f16x4);
                vC3 = AE_MULFC32X16RAS_L(vC3, t3_f16x4);

                AE_S32X2_X(vC3, py1, 3 * 8 * stride);
                AE_S32X2_X(vC2, py1, 1 * 8 * stride);
                AE_S32X2_X(vC1, py1, 2 * 8 * stride);
                AE_S32X2_IP(vC0, py1, 2 * 8);
            } // for(i=1; i<N/8; i++);
            stride >>= 2;
        } // End of the first stage.

        s = 2;
        WUR_AE_SAR(s);

        WUR_AE_CBEGIN0((unsigned)y);
        WUR_AE_CEND0((unsigned)(uintptr_t)(((complex_fract32*)y) + (N - 2)));

        while (stride > 4)
        {

            
            p16tw1 = (const ae_int32x2*)*seq++;
            //!!!! new
            unsigned acc = 0;
            const unsigned acc_inc = (log2n == 0) ? 0 : (0x80000000 >> (log2n - 1));

            //!!!!
            log2n += 2;
            shift += s;
        
            py0 = (ae_int32x2*)y;
            py1 = py0 + 1;
            __Pragma("loop_count min=1");
            for (i = 0; i < (N>>3); i ++)
            {
                ae_int32x2 t1, t2, t3;
                ae_f16x4 t1_f16x4;
                ae_f16x4 t2_f16x4;
                ae_f16x4 t3_f16x4;

                int offset_inc = 0;
                acc += acc_inc;
                XT_MOVEQZ(offset_inc, 3 * 8, acc);

                t3 = AE_L32X2_I(p16tw1, 8 * 2);
                t2 = AE_L32X2_I(p16tw1, 8 * 1);
                AE_L32X2_XP(t1, p16tw1, offset_inc);

                t1_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
                t2_f16x4 = AE_MOVF16X4_FROMINT32X2(t2);
                t3_f16x4 = AE_MOVF16X4_FROMINT32X2(t3);

                px0=(const ae_int32x2*)py0;
                px1=(const ae_int32x2*)py1;
                vA3 = AE_L32X2_X(px0, 3 * 8 * stride);
                vA2 = AE_L32X2_X(px0, 2 * 8 * stride);
                vA1 = AE_L32X2_X(px0, 1 * 8 * stride);
                vA0=AE_L32X2_I     ( px0,0);

                vA3 = AE_SRAS32(vA3);
                vA2 = AE_SRAS32(vA2);
                vA1 = AE_SRAS32(vA1);
                vA0 = AE_SRAS32(vA0);

                vB0 = AE_ADD32S(vA0, vA2);
                vB2 = AE_SUB32S(vA0, vA2);
                vB1 = AE_ADD32S(vA1, vA3);
                vB3 = AE_SUB32S(vA1, vA3);

                vB3 = AE_SEL32_LH(vB3, vB3);

                vC0 = AE_ADD32S(vB0, vB1);
                vC2 = AE_SUB32S(vB0, vB1);
                vC3 = AE_SUBADD32S(vB2, vB3);
                vC1 = AE_ADDSUB32S(vB2, vB3);

                vC1 = AE_MULFC32X16RAS_H(vC1, t1_f16x4);
                vC2 = AE_MULFC32X16RAS_H(vC2, t2_f16x4);
                vC3 = AE_MULFC32X16RAS_H(vC3, t3_f16x4);

                AE_S32X2_X(vC3,  py0, 3 * 8 * stride);
                AE_S32X2_X(vC2,  py0, 1 * 8 * stride);
                AE_S32X2_X(vC1,  py0, 2 * 8 * stride);
                AE_S32X2_XC(vC0, py0, 4 * 8 * stride);
                vA3 = AE_L32X2_X(px1, 3 * 8 * stride);
                vA2 = AE_L32X2_X(px1, 2 * 8 * stride);
                vA1 = AE_L32X2_X(px1, 1 * 8 * stride);
                vA0=AE_L32X2_I(px1, 0);

                vA3 = AE_SRAS32(vA3);
                vA2 = AE_SRAS32(vA2);
                vA1 = AE_SRAS32(vA1);
                vA0 = AE_SRAS32(vA0);

                vB0 = AE_ADD32S(vA0, vA2);
                vB2 = AE_SUB32S(vA0, vA2);
                vB1 = AE_ADD32S(vA1, vA3);
                vB3 = AE_SUB32S(vA1, vA3);

                vB3 = AE_SEL32_LH(vB3, vB3);

                vC0 = AE_ADD32S(vB0, vB1);
                vC2 = AE_SUB32S(vB0, vB1);
                vC3 = AE_SUBADD32S(vB2, vB3);
                vC1 = AE_ADDSUB32S(vB2, vB3);

                vC1 = AE_MULFC32X16RAS_L(vC1, t1_f16x4);
                vC2 = AE_MULFC32X16RAS_L(vC2, t2_f16x4);
                vC3 = AE_MULFC32X16RAS_L(vC3, t3_f16x4);
                AE_S32X2_X( vC3,  py1, 3 * 8 * stride);
                AE_S32X2_X (vC2,  py1, 1 * 8 * stride);
                AE_S32X2_X (vC1,  py1, 2 * 8 * stride);
                AE_S32X2_XC(vC0,  py1, 4 * 8 * stride);
            } // for(i=0; i<N/8; i++);
            stride >>= 2;
        } // while (stride > 4)

        if (stride > 1)
        {
            int j;
            p16tw1 = (ae_int32x2*)*seq++;           
            shift += s;

            px0 = (ae_int32x2*)y;
            px1 = px0 + 1;
            py0 = (ae_int32x2*)x;
            py1 = py0 + 1;
            __Pragma("loop_count min=1");
            for (j = 0; j < stride ; j+=2)
            {
                ae_int32x2 t1, t2, t3;
                ae_f16x4 t1_f16x4;
                ae_f16x4 t2_f16x4;
                ae_f16x4 t3_f16x4;

                px0 = j + (ae_int32x2*)y;
                px1 = px0 + 1;
                py0 = j + (ae_int32x2*)x;
                py1 = py0 + 1;

                AE_L32X2_IP(t1, p16tw1, 8);
                AE_L32X2_IP(t2, p16tw1, 8);
                AE_L32X2_IP(t3, p16tw1, 8);   

                t1_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
                t2_f16x4 = AE_MOVF16X4_FROMINT32X2(t2);
                t3_f16x4 = AE_MOVF16X4_FROMINT32X2(t3);
                /*
                17 cycles per pipeline stage in steady state with unroll=1
                8.5 cycles per DFT4
                */
                __Pragma("loop_count min=1");
                for (i = 0; i < N / (4*stride); i++)
                {
                //    px = (ae_int32x2*)((i + j * 4 * stride) * sizeof(*px) + (uintptr_t)x);

                    vA3 = AE_L32X2_X(px0, 3 * 8 * stride);
                    vA2 = AE_L32X2_X(px0, 2 * 8 * stride);
                    vA1 = AE_L32X2_X(px0, 1 * 8 * stride);
                    AE_L32X2_XP(vA0, px0, 4 * 8 * stride);

                    vA3 = AE_SRAS32(vA3);
                    vA2 = AE_SRAS32(vA2);
                    vA1 = AE_SRAS32(vA1);
                    vA0 = AE_SRAS32(vA0);

                    vB0 = AE_ADD32S(vA0, vA2);
                    vB2 = AE_SUB32S(vA0, vA2);
                    vB1 = AE_ADD32S(vA1, vA3);
                    vB3 = AE_SUB32S(vA1, vA3);

                    vB3 = AE_SEL32_LH(vB3, vB3);

                    vC0 = AE_ADD32S(vB0, vB1);
                    vC2 = AE_SUB32S(vB0, vB1);
                    vC3 = AE_SUBADD32S(vB2, vB3);
                    vC1 = AE_ADDSUB32S(vB2, vB3);

                    vC1 = AE_MULFC32X16RAS_H(vC1, t1_f16x4);
                    vC2 = AE_MULFC32X16RAS_H(vC2, t2_f16x4);
                    vC3 = AE_MULFC32X16RAS_H(vC3, t3_f16x4);

                    AE_S32X2_X(vC3, py0, 3 * 8 * stride);
                    AE_S32X2_X(vC2, py0, 1 * 8 * stride);
                    AE_S32X2_X(vC1, py0, 2 * 8 * stride);
                    AE_S32X2_XP(vC0, py0, 4 * 8 * stride);

                    vA3 = AE_L32X2_X(px1, 3 * 8 * stride);
                    vA2 = AE_L32X2_X(px1, 2 * 8 * stride);
                    vA1 = AE_L32X2_X(px1, 1 * 8 * stride);
                    AE_L32X2_XP(vA0, px1, 4 * 8 * stride);

                    vA3 = AE_SRAS32(vA3);
                    vA2 = AE_SRAS32(vA2);
                    vA1 = AE_SRAS32(vA1);
                    vA0 = AE_SRAS32(vA0);

                    vB0 = AE_ADD32S(vA0, vA2);
                    vB2 = AE_SUB32S(vA0, vA2);
                    vB1 = AE_ADD32S(vA1, vA3);
                    vB3 = AE_SUB32S(vA1, vA3);

                    vB3 = AE_SEL32_LH(vB3, vB3);

                    vC0 = AE_ADD32S(vB0, vB1);
                    vC2 = AE_SUB32S(vB0, vB1);
                    vC3 = AE_SUBADD32S(vB2, vB3);
                    vC1 = AE_ADDSUB32S(vB2, vB3);

                    vC1 = AE_MULFC32X16RAS_L(vC1, t1_f16x4);
                    vC2 = AE_MULFC32X16RAS_L(vC2, t2_f16x4);
                    vC3 = AE_MULFC32X16RAS_L(vC3, t3_f16x4);

                    AE_S32X2_X(vC3, py1, 3 * 8 * stride);
                    AE_S32X2_X(vC2, py1, 1 * 8 * stride);
                    AE_S32X2_X(vC1, py1, 2 * 8 * stride);
                    AE_S32X2_XP(vC0, py1, 4 * 8 * stride);
                } //  (i = 0; i < N / (4*stride); i++)
            }
            stride >>= 2;
        }       
        if (N == 16)
        {
            shift += fft_stage_last_ie((int32_t*)y, (int32_t*)x, N);
        }
        else
            shift += fft_stage_last_ie((int32_t*)x, (int32_t*)y, N);

    }
   
    return shift;
}
#endif // #if XCHAL_HAVE_HIFI3Z
