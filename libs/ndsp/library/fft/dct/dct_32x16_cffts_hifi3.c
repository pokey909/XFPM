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

#define FFT_BUTTERFLY_S3_T3F(_T, _step)    \
      {                                    \
      vA0 = AE_L32X2_I(p_x0, 8);           \
      vC0 = AE_ADD32S(vB0, vB1);           \
      vC2 = AE_SUBADD32S(vB2, vB3);        \
      vA1 = AE_L32X2_I(p_x1, 8);           \
      vC1 = AE_MULFC32X16RAS_L(vC1, vT0);  \
      vC3 = AE_ADDSUB32S(vB2, vB3);        \
      vA2 = AE_L32X2_I(p_x2, 8);           \
      vC2 = AE_MULFC32X16RAS_H(vC2, vT1);  \
      vC0 = AE_SRAS32(vC0);                \
      vA3 = AE_L32X2_I(p_x3, 8);           \
      vC3 = AE_MULFC32X16RAS_H(vC3, vT0);  \
      vC1 = AE_SRAS32(vC1);                \
      AE_S32X2_IP(vC0, p_x0, 8);           \
      vB0 = AE_ADD32S(vA0, vA2);           \
      vC2 = AE_SRAS32(vC2);                \
      AE_S32X2_IP(vC1, p_x1, 8);           \
      vB2 = AE_SUB32S(vA0, vA2);           \
      vC3 = AE_SRAS32(vC3);                \
      AE_S32X2_IP(vC2, p_x3, 8);           \
      vB1 = AE_ADD32S(vA1, vA3);           \
      vB3 = AE_SUB32S(vA1, vA3);           \
      AE_S32X2_IP(vC3, p_x2, 8);           \
      vC1 = AE_SUB32S(vB0, vB1);           \
      vB3 = AE_SEL32_LH(vB3, vB3);         \
      vA0 = AE_L32X2_##_T(p_x0, _step);    \
      vC0 = AE_ADD32S(vB0, vB1);           \
      vC2 = AE_SUBADD32S(vB2, vB3);        \
      vA1 = AE_L32X2_##_T(p_x1, _step);    \
      vC1 = AE_MULFC32X16RAS_H(vC1, vT2);  \
      vC3 = AE_ADDSUB32S(vB2, vB3);        \
      vA2 = AE_L32X2_##_T(p_x2, _step);    \
      vC2 = AE_MULFC32X16RAS_L(vC2, vT2);  \
      vC0 = AE_SRAS32(vC0);                \
      vA3 = AE_L32X2_##_T(p_x3, _step);    \
      vC3 = AE_MULFC32X16RAS_L(vC3, vT1);  \
      vC1 = AE_SRAS32(vC1);                \
      AE_S32X2_##_T##P(vC0, p_x0, _step);  \
      vB0 = AE_ADD32S(vA0, vA2);           \
      vC2 = AE_SRAS32(vC2);                \
      AE_S32X2_##_T##P(vC1, p_x1, _step);  \
      vB3 = AE_SUB32S(vA1, vA3);           \
      vC3 = AE_SRAS32(vC3);                \
      AE_S32X2_##_T##P(vC2, p_x3, _step);  \
      vB2 = AE_SUB32S(vA0, vA2);           \
      vB1 = AE_ADD32S(vA1, vA3);           \
      AE_S32X2_##_T##P(vC3, p_x2, _step);  \
      vB3 = AE_SEL32_LH(vB3, vB3);         \
      vC1 = AE_SUB32S(vB0, vB1);           \
      }

#define FFT_BUTTERFLY_S3_T3F_LAST(_T, _step)    \
      {                                    \
      vA0 = AE_L32X2_I(p_x0, 8);           \
      vC0 = AE_ADD32S(vB0, vB1);           \
      vC2 = AE_SUBADD32S(vB2, vB3);        \
      vA1 = AE_L32X2_I(p_x1, 8);           \
      vC1 = AE_MULFC32X16RAS_L(vC1, vT0);  \
      vC3 = AE_ADDSUB32S(vB2, vB3);        \
      vA2 = AE_L32X2_I(p_x2, 8);           \
      vC2 = AE_MULFC32X16RAS_H(vC2, vT1);  \
      vC0 = AE_SRAS32(vC0);                \
      vA3 = AE_L32X2_I(p_x3, 8);           \
      vC3 = AE_MULFC32X16RAS_H(vC3, vT0);  \
      vC1 = AE_SRAS32(vC1);                \
      AE_S32X2_IP(vC0, p_x0, 8);           \
      vB0 = AE_ADD32S(vA0, vA2);           \
      vC2 = AE_SRAS32(vC2);                \
      AE_S32X2_IP(vC1, p_x1, 8);           \
      vB2 = AE_SUB32S(vA0, vA2);           \
      vC3 = AE_SRAS32(vC3);                \
      AE_S32X2_IP(vC2, p_x3, 8);           \
      vB1 = AE_ADD32S(vA1, vA3);           \
      vB3 = AE_SUB32S(vA1, vA3);           \
      AE_S32X2_IP(vC3, p_x2, 8);           \
      vC1 = AE_SUB32S(vB0, vB1);           \
      vB3 = AE_SEL32_LH(vB3, vB3);         \
      vC0 = AE_ADD32S(vB0, vB1);           \
      vC2 = AE_SUBADD32S(vB2, vB3);        \
      vC1 = AE_MULFC32X16RAS_H(vC1, vT2);  \
      vC3 = AE_ADDSUB32S(vB2, vB3);        \
      vC2 = AE_MULFC32X16RAS_L(vC2, vT2);  \
      vC0 = AE_SRAS32(vC0);                \
      vC3 = AE_MULFC32X16RAS_L(vC3, vT1);  \
      vC1 = AE_SRAS32(vC1);                \
      AE_S32X2_##_T##P(vC0, p_x0, _step);  \
      vC2 = AE_SRAS32(vC2);                \
      AE_S32X2_##_T##P(vC1, p_x1, _step);  \
      vC3 = AE_SRAS32(vC3);                \
      AE_S32X2_##_T##P(vC2, p_x3, _step);  \
      AE_S32X2_##_T##P(vC3, p_x2, _step);  \
      }

#define FFT_BUTTERFLY_R2(idx)       \
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
void fft16_32x16(int32_t *y, int32_t *x, const int16_t *ptwd)
{
  ae_int32x2  * p_x0, * p_x1, * restrict p_y0;
  const ae_int16x4  * restrict p_twd;
  ae_int32x2  vA0, vA1, vA2, vA3, vB0, vB1, vB2, vB3;
  ae_f16x4    vT1, vT2;
  ae_int16x4  vTT;

  const int step_right   = (int)sizeof(ae_int32x2);
  const int step_left    = -step_right;
  const int step_down    = 4*step_right;
  const int step_nextc   = step_right - 3*step_down;
  const int step_hupleft = 7*step_left;
  const int step_upleft  = 15*step_left;

  p_x0  = (ae_int32x2 *)(x);
  p_x1  = (ae_int32x2 *)(x+2*15);
  p_y0  = (ae_int32x2 *)(y);
  p_twd = (const ae_int16x4 *)(ptwd)+1;

  /*** reordering ***/
  AE_L32X2_XP(vA0, p_y0, step_right);
  AE_L32X2_XP(vA1, p_y0, step_right);
  AE_L32X2_XP(vA2, p_y0, step_right);
  AE_L32X2_XP(vA3, p_y0, step_right);

  vB0 = AE_SEL32_HH(vA0, vA1);
  vB3 = AE_SEL32_LL(vA1, vA0);
  AE_L32X2_XP(vA0, p_y0, step_right);
  AE_L32X2_XP(vA1, p_y0, step_right);
  vB1 = AE_SEL32_HH(vA2, vA3);
  vB2 = AE_SEL32_LL(vA3, vA2);
  AE_L32X2_XP(vA2, p_y0, step_right);
  AE_L32X2_XP(vA3, p_y0, step_right);
  AE_S32X2_XP(vB0, p_x0, step_right);
  AE_S32X2_XP(vB3, p_x1, step_left);
  AE_S32X2_XP(vB1, p_x0, step_right);
  AE_S32X2_XP(vB2, p_x1, step_left);

  vB0 = AE_SEL32_HH(vA0, vA1);
  vB3 = AE_SEL32_LL(vA1, vA0);
  AE_L32X2_XP(vA0, p_y0, step_right);
  AE_L32X2_XP(vA1, p_y0, step_right);
  vB1 = AE_SEL32_HH(vA2, vA3);
  vB2 = AE_SEL32_LL(vA3, vA2);
  AE_L32X2_XP(vA2, p_y0, step_right);
  AE_L32X2_XP(vA3, p_y0, step_right);
  AE_S32X2_XP(vB0, p_x0, step_right);
  AE_S32X2_XP(vB3, p_x1, step_left);
  AE_S32X2_XP(vB1, p_x0, step_right);
  AE_S32X2_XP(vB2, p_x1, step_left);

  vB0 = AE_SEL32_HH(vA0, vA1);
  vB3 = AE_SEL32_LL(vA1, vA0);
  AE_L32X2_XP(vA0, p_y0, step_right);
  AE_L32X2_XP(vA1, p_y0, step_right);
  vB1 = AE_SEL32_HH(vA2, vA3);
  vB2 = AE_SEL32_LL(vA3, vA2);
  AE_L32X2_XP(vA2, p_y0, step_right);
  AE_L32X2_XP(vA3, p_y0, step_upleft);
  AE_S32X2_XP(vB0, p_x0, step_right);
  AE_S32X2_XP(vB3, p_x1, step_left);
  AE_S32X2_XP(vB1, p_x0, step_right);
  AE_S32X2_XP(vB2, p_x1, step_left);

  vB0 = AE_SEL32_HH(vA0, vA1);
  vB3 = AE_SEL32_LL(vA1, vA0);
  vB1 = AE_SEL32_HH(vA2, vA3);
  vB2 = AE_SEL32_LL(vA3, vA2);
  AE_S32X2_XP(vB0, p_x0, step_right);
  AE_S32X2_XP(vB3, p_x1, step_left);
  AE_S32X2_XP(vB1, p_x0, step_hupleft);
  AE_S32X2_XP(vB2, p_x1, 0);

  /*** first radix4 stage ***/
  /* Elements 0, 4, 8, 12 */
  vA0 = AE_L32X2_I(p_x1, (0-8)*(int)sizeof(ae_int32x2));
  vA1 = AE_L32X2_I(p_x1, (4-8)*(int)sizeof(ae_int32x2));
  vA2 = AE_L32X2_I(p_x1, (8-8)*(int)sizeof(ae_int32x2));
  vA3 = AE_L32X2_I(p_x1, (12-8)*(int)sizeof(ae_int32x2));

  vA0 = AE_SRAI32(vA0, 3);
  vA1 = AE_SRAI32(vA1, 3);
  vA2 = AE_SRAI32(vA2, 3);
  vA3 = AE_SRAI32(vA3, 3);

  vB0 = AE_ADD32S(vA0, vA2);
  vB2 = AE_SUB32S(vA0, vA2);
  vB1 = AE_ADD32S(vA1, vA3);
  vA2 = AE_SUB32S(vA3, vA1);
  vB3 = AE_SUB32S(vA1, vA3);
  vB3 = AE_SEL32_LH(vA2, vB3);

  vA0 = AE_ADD32S(vB0, vB1);
  vA1 = AE_SUB32S(vB0, vB1);
  vA2 = AE_ADD32S(vB2, vB3);
  vA3 = AE_SUB32S(vB2, vB3);

  vB0 = AE_SRAI32(vA0, 2);
  vB1 = AE_SRAI32(vA1, 2);
  vB2 = AE_SRAI32(vA2, 2);
  vB3 = AE_SRAI32(vA3, 2);

  AE_S32X2_I(vB0, p_x1, (0-8)*(int)sizeof(ae_int32x2));
  AE_S32X2_I(vB3, p_x1, (4-8)*(int)sizeof(ae_int32x2));
  AE_S32X2_I(vB1, p_x1, (8-8)*(int)sizeof(ae_int32x2));
  AE_S32X2_I(vB2, p_x1, (12-8)*(int)sizeof(ae_int32x2));

  /* Elements 1, 5, 9, 13 */
  vA0 = AE_L32X2_I(p_x1, (1-8)*(int)sizeof(ae_int32x2));
  vA1 = AE_L32X2_I(p_x1, (5-8)*(int)sizeof(ae_int32x2));
  vA2 = AE_L32X2_I(p_x1, (9-8)*(int)sizeof(ae_int32x2));
  vA3 = AE_L32X2_I(p_x1, (13-8)*(int)sizeof(ae_int32x2));

  vA0 = AE_SRAI32(vA0, 3);
  vA1 = AE_SRAI32(vA1, 3);
  vA2 = AE_SRAI32(vA2, 3);
  vA3 = AE_SRAI32(vA3, 3);

  vB0 = AE_ADD32S(vA0, vA2);
  vB2 = AE_SUB32S(vA0, vA2);
  vB1 = AE_ADD32S(vA1, vA3);
  vA2 = AE_SUB32S(vA3, vA1);
  vB3 = AE_SUB32S(vA1, vA3);
  vB3 = AE_SEL32_LH(vA2, vB3);

  vA0 = AE_ADD32S(vB0, vB1);
  vA1 = AE_SUB32S(vB0, vB1);
  vA2 = AE_ADD32S(vB2, vB3);
  vA3 = AE_SUB32S(vB2, vB3);

  AE_L16X4_IP(vTT, p_twd, 8);
  vT1 = (vTT);
  AE_L16X4_IP(vTT, p_twd, 8);
  vT2 = (vTT);

  vA3 = AE_MULFC32X16RAS_L(vA3, vT1);
  vA1 = AE_MULFC32X16RAS_H(vA1, vT2);
  vA2 = AE_MULFC32X16RAS_L(vA2, vT2);

  vB0 = AE_SRAI32(vA0, 2);
  vB1 = AE_SRAI32(vA1, 2);
  vB2 = AE_SRAI32(vA2, 2);
  vB3 = AE_SRAI32(vA3, 2);

  AE_S32X2_I(vB0, p_x1, (1-8)*(int)sizeof(ae_int32x2));
  AE_S32X2_I(vB3, p_x1, (5-8)*(int)sizeof(ae_int32x2));
  AE_S32X2_I(vB1, p_x1, (9-8)*(int)sizeof(ae_int32x2));
  AE_S32X2_I(vB2, p_x1, (13-8)*(int)sizeof(ae_int32x2));

  /* Elements 2, 6, 10, 14 */
  vA0 = AE_L32X2_I(p_x1, (2-8)*(int)sizeof(ae_int32x2));
  vA1 = AE_L32X2_I(p_x1, (6-8)*(int)sizeof(ae_int32x2));
  vA2 = AE_L32X2_I(p_x1, (10-8)*(int)sizeof(ae_int32x2));
  vA3 = AE_L32X2_I(p_x1, (14-8)*(int)sizeof(ae_int32x2));

  vA0 = AE_SRAI32(vA0, 3);
  vA1 = AE_SRAI32(vA1, 3);
  vA2 = AE_SRAI32(vA2, 3);
  vA3 = AE_SRAI32(vA3, 3);

  vB0 = AE_ADD32S(vA0, vA2);
  vB2 = AE_SUB32S(vA0, vA2);
  vB1 = AE_ADD32S(vA1, vA3);
  vA2 = AE_SUB32S(vA3, vA1);
  vB3 = AE_SUB32S(vA1, vA3);
  vB3 = AE_SEL32_LH(vA2, vB3);

  vA0 = AE_ADD32S(vB0, vB1);
  vA1 = AE_SUB32S(vB0, vB1);
  vA2 = AE_ADD32S(vB2, vB3);
  vA3 = AE_SUB32S(vB2, vB3);

  AE_L16X4_IP(vTT, p_twd, 8);
  vT1 = (vTT);
  AE_L16X4_IP(vTT, p_twd, 8);
  vT2 = (vTT);

  vA3 = AE_MULFC32X16RAS_H(vA3, vT1);
  vA1 = AE_MULFC32X16RAS_L(vA1, vT1);
  vA2 = AE_MULFC32X16RAS_H(vA2, vT2);

  vB0 = AE_SRAI32(vA0, 2);
  vB1 = AE_SRAI32(vA1, 2);
  vB2 = AE_SRAI32(vA2, 2);
  vB3 = AE_SRAI32(vA3, 2);

  AE_S32X2_I(vB0, p_x1, (2-8)*(int)sizeof(ae_int32x2));
  AE_S32X2_I(vB3, p_x1, (6-8)*(int)sizeof(ae_int32x2));
  AE_S32X2_I(vB1, p_x1, (10-8)*(int)sizeof(ae_int32x2));
  AE_S32X2_I(vB2, p_x1, (14-8)*(int)sizeof(ae_int32x2));

  /* Elements 3, 7, 11, 15 */
  vA0 = AE_L32X2_I(p_x1, (3-8)*(int)sizeof(ae_int32x2));
  vA1 = AE_L32X2_I(p_x1, (7-8)*(int)sizeof(ae_int32x2));
  vA2 = AE_L32X2_I(p_x1, (11-8)*(int)sizeof(ae_int32x2));
  vA3 = AE_L32X2_I(p_x1, (15-8)*(int)sizeof(ae_int32x2));

  vA0 = AE_SRAI32(vA0, 3);
  vA1 = AE_SRAI32(vA1, 3);
  vA2 = AE_SRAI32(vA2, 3);
  vA3 = AE_SRAI32(vA3, 3);

  vB0 = AE_ADD32S(vA0, vA2);
  vB2 = AE_SUB32S(vA0, vA2);
  vB1 = AE_ADD32S(vA1, vA3);
  vA2 = AE_SUB32S(vA3, vA1);
  vB3 = AE_SUB32S(vA1, vA3);
  vB3 = AE_SEL32_LH(vA2, vB3);

  vA0 = AE_ADD32S(vB0, vB1);
  vA1 = AE_SUB32S(vB0, vB1);
  vA2 = AE_ADD32S(vB2, vB3);
  vA3 = AE_SUB32S(vB2, vB3);

  AE_L16X4_IP(vTT, p_twd, 8);
  vT1 = (vTT);

  vA3 = AE_MULFC32X16RAS_L(vA3, vT2);
  vA1 = AE_MULFC32X16RAS_H(vA1, vT1);
  vA2 = AE_MULFC32X16RAS_L(vA2, vT1);

  vB0 = AE_SRAI32(vA0, 2);
  vB1 = AE_SRAI32(vA1, 2);
  vB2 = AE_SRAI32(vA2, 2);
  vB3 = AE_SRAI32(vA3, 2);

  AE_S32X2_I(vB0, p_x1, (3-8)*(int)sizeof(ae_int32x2));
  AE_S32X2_I(vB3, p_x1, (7-8)*(int)sizeof(ae_int32x2));
  AE_S32X2_I(vB1, p_x1, (11-8)*(int)sizeof(ae_int32x2));
  AE_S32X2_I(vB2, p_x1, (15-8)*(int)sizeof(ae_int32x2));

  /*** last radix4 stage ***/
  AE_L32X2_XP(vA0, p_x0, step_right);
  AE_L32X2_XP(vA1, p_x0, step_right);
  AE_L32X2_XP(vA2, p_x0, step_right);
  AE_L32X2_XP(vA3, p_x0, step_right);

  vB0 = AE_ADD32S(vA0, vA2);
  vB2 = AE_SUB32S(vA0, vA2);
  vB1 = AE_ADD32S(vA1, vA3);
  vA2 = AE_SUB32S(vA3, vA1);
  vB3 = AE_SUB32S(vA1, vA3);
  vB3 = AE_SEL32_LH(vA2, vB3);

  vA0 = AE_ADD32S(vB0, vB1);
  vA2 = AE_SUB32S(vB0, vB1);
  vA1 = AE_SUB32S(vB2, vB3);
  vA3 = AE_ADD32S(vB2, vB3);

  AE_S32X2_XP(vA0, p_y0, step_down);
  AE_S32X2_XP(vA1, p_y0, step_down);
  AE_S32X2_XP(vA2, p_y0, step_down);
  AE_S32X2_XP(vA3, p_y0, step_nextc);

  AE_L32X2_XP(vA0, p_x0, step_right);
  AE_L32X2_XP(vA1, p_x0, step_right);
  AE_L32X2_XP(vA2, p_x0, step_right);
  AE_L32X2_XP(vA3, p_x0, step_right);

  vB0 = AE_ADD32S(vA0, vA2);
  vB2 = AE_SUB32S(vA0, vA2);
  vB1 = AE_ADD32S(vA1, vA3);
  vA2 = AE_SUB32S(vA3, vA1);
  vB3 = AE_SUB32S(vA1, vA3);
  vB3 = AE_SEL32_LH(vA2, vB3);

  vA0 = AE_ADD32S(vB0, vB1);
  vA2 = AE_SUB32S(vB0, vB1);
  vA1 = AE_SUB32S(vB2, vB3);
  vA3 = AE_ADD32S(vB2, vB3);

  AE_S32X2_XP(vA0, p_y0, step_down);
  AE_S32X2_XP(vA1, p_y0, step_down);
  AE_S32X2_XP(vA2, p_y0, step_down);
  AE_S32X2_XP(vA3, p_y0, step_nextc);

  AE_L32X2_XP(vA0, p_x0, step_right);
  AE_L32X2_XP(vA1, p_x0, step_right);
  AE_L32X2_XP(vA2, p_x0, step_right);
  AE_L32X2_XP(vA3, p_x0, step_right);

  vB0 = AE_ADD32S(vA0, vA2);
  vB2 = AE_SUB32S(vA0, vA2);
  vB1 = AE_ADD32S(vA1, vA3);
  vA2 = AE_SUB32S(vA3, vA1);
  vB3 = AE_SUB32S(vA1, vA3);
  vB3 = AE_SEL32_LH(vA2, vB3);

  vA0 = AE_ADD32S(vB0, vB1);
  vA2 = AE_SUB32S(vB0, vB1);
  vA1 = AE_SUB32S(vB2, vB3);
  vA3 = AE_ADD32S(vB2, vB3);

  AE_S32X2_XP(vA0, p_y0, step_down);
  AE_S32X2_XP(vA1, p_y0, step_down);
  AE_S32X2_XP(vA2, p_y0, step_down);
  AE_S32X2_XP(vA3, p_y0, step_nextc);

  AE_L32X2_XP(vA0, p_x0, step_right);
  AE_L32X2_XP(vA1, p_x0, step_right);
  AE_L32X2_XP(vA2, p_x0, step_right);
  AE_L32X2_XP(vA3, p_x0, step_right);

  vB0 = AE_ADD32S(vA0, vA2);
  vB2 = AE_SUB32S(vA0, vA2);
  vB1 = AE_ADD32S(vA1, vA3);
  vA2 = AE_SUB32S(vA3, vA1);
  vB3 = AE_SUB32S(vA1, vA3);
  vB3 = AE_SEL32_LH(vA2, vB3);

  vA0 = AE_ADD32S(vB0, vB1);
  vA2 = AE_SUB32S(vB0, vB1);
  vA1 = AE_SUB32S(vB2, vB3);
  vA3 = AE_ADD32S(vB2, vB3);

  AE_S32X2_XP(vA0, p_y0, step_down);
  AE_S32X2_XP(vA1, p_y0, step_down);
  AE_S32X2_XP(vA2, p_y0, step_down);
  AE_S32X2_I(vA3, p_y0, 0);
}

/*
   scaled fft with reordering
   N=32 - size of FFT
   NOTE: y is input and output, x - temporary
*/
void fft32_32x16(int32_t *y, int32_t *x, const int16_t *ptwd)
{
    int k;
    const int N = 32;
    ae_int32x2  * restrict p_x0, * restrict p_x1,
                * restrict p_x2, * restrict p_x3,
                * restrict p_y0;
    const ae_int16x4 * restrict p_twd;
    ae_int32x2  vA0, vA1, vA2, vA3, vB0, vB1, vB2, vB3;

    p_x0 = (ae_int32x2 *)(x);
    p_x1 = (ae_int32x2 *)(x+2*N-2);
    p_y0 = (ae_int32x2 *)(y);

    // reordering
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
#if 0
    fft_cplx32x16(y, x, cfft16_32, 3);
#else
    {
      int lc;

      ae_int32x2  vC0, vC1, vC2, vC3;
      ae_int16x4  vTL;
      ae_f16x4    vT0, vT1, vT2;
      int stride;
      int rsa;

      p_twd = (const ae_int16x4 *)ptwd;

      stride = N; // The stride is quartered with every iteration of the outer loop.
      rsa = 2;

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
          vC2 = AE_SUBADD32S(vB2, vB3);
          vA0 = AE_L32X2_I(p_x0, 8);
          vC3 = AE_ADDSUB32S(vB2, vB3);
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
          vC2 = AE_SUBADD32S(vB2, vB3);
          vA0 = AE_SRAI32(vA0, FIRST_STAGE_SCALE);
          vA2 = AE_L32X2_I(p_x2, 8);
          vC2 = AE_MULFC32X16RAS_L(vC2, vT2);
          vC3 = AE_ADDSUB32S(vB2, vB3);
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
          vC2 = AE_SUBADD32S(vB2, vB3);
          vA0 = AE_L32X2_I(p_x0, 8);
          vC3 = AE_ADDSUB32S(vB2, vB3);
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
          vC2 = AE_SUBADD32S(vB2, vB3);
          vC2 = AE_MULFC32X16RAS_L(vC2, vT2);
          vC3 = AE_ADDSUB32S(vB2, vB3);
          vC3 = AE_MULFC32X16RAS_L(vC3, vT1);
          AE_S32X2_IP(vC0, p_x0, 8);
          vC1 = AE_SRAS32(vC1);
          AE_S32X2_IP(vC1, p_x1, 8);
          vT1 = AE_L16X4_I(p_twd, 8);
          vC2 = AE_SRAS32(vC2);
          AE_S32X2_IP(vC2, p_x3, 8);
          vC3 = AE_SRAS32(vC3);
          AE_S32X2_IP(vC3, p_x2, 8);
        }
        stride  >>= 2;
      }

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
        vT0 = AE_L16X4_I(((const ae_int16x4*)fft_twd_r8), 0);
        vT1 = AE_L16X4_I(((const ae_int16x4*)fft_twd_r8), 8);
        vT2 = AE_L16X4_I(((const ae_int16x4*)fft_twd_r8), 16);

        vB0 = AE_ADD32S(vA0, vA2);
        vB3 = AE_SUB32S(vA1, vA3);
        vB2 = AE_SUB32S(vA0, vA2);
        vB1 = AE_ADD32S(vA1, vA3);
        vB3 = AE_SEL32_LH(vB3, vB3);
        vC1 = AE_SUB32S(vB0, vB1);

        for(lc=0; lc<((N>>3)-1); lc++)
        {
          FFT_BUTTERFLY_S3_T3F(I, 56)
        };
        FFT_BUTTERFLY_S3_T3F_LAST(I, 56)

        stride >>= 2;
      }

      {
        int32_t i,i0,i1,i2,i3,ai;
        ae_int32x2 *p_y0 = (ae_int32x2 *)(y);
        ae_int32x2 *p_y1 = (p_y0 + (N >> 2));
        ae_int32x2 *p_y2 = (p_y1 + (N >> 2));
        ae_int32x2 *p_y3 = (p_y2 + (N >> 2));
        p_x0 = (ae_int32x2 *)(x);

        i = NSA(N)+1;
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
#endif
}
