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

#define FFT_BUTTERFLY_S3_0(_T,_step,_X,twd_step)\
      {                                         \
        vT3 = AE_L32X2F24_I(p_twd, 16);         \
        vC0 = AE_ADD32S(vB0, vB1);      vC0 = AE_SRAI32R(vC0, 3);         \
        vB3 = AE_SEL32_LH(vB3, vB3);            \
        AE_L32X2F24_##_X(vT1, p_twd, twd_step); \
        vA2 = AE_SUBADD32S(vB2, vB3);   vA2 = AE_SRAI32R(vA2, 3);          \
        vF1 = AE_L32X2F24_##_T(p_x0, _step);    \
        vA0 = (vF1);                            \
        vA3 = AE_ADDSUB32S(vB2, vB3);   vA3 = AE_SRAI32R(vA3, 3);          \
                    \
        vF1 = AE_MOVF24X2_FROMINT32X2(vA1);     \
        vF0 = AE_MULFC24RA(vF1, vT2);           \
        vC1 = (vF0);                            \
        vF2 = AE_L32X2F24_##_T(p_x1, _step);    \
        vA1 = (vF2);                            \
        /*vA0 = AE_SRAI32(vA0, 3); */               \
        vF1 = AE_MOVF24X2_FROMINT32X2(vA2);     \
        vF0 = AE_MULFC24RA(vF1, vT3);           \
        vC2 = (vF0);                            \
        vF1 = AE_L32X2F24_##_T(p_x2, _step);    \
        vA2 = (vF1);                            \
        /*vA1 = AE_SRAI32(vA1, 3);     */           \
        vF1 = AE_MOVF24X2_FROMINT32X2(vA3);     \
        vF0 = AE_MULFC24RA(vF1, vT1);           \
        vC3 = (vF0);                            \
        vF2 = AE_L32X2F24_##_T(p_x3, _step);    \
        vA3 = (vF2);                            \
        /*vA2 = AE_SRAI32(vA2, 3);  */              \
        vT2 = AE_L32X2F24_I(p_twd, 8);          \
        vB0 = AE_ADD32S(vA0, vA2);              \
        /*vA3 = AE_SRAI32(vA3, 3); */               \
        vF1 = AE_MOVF24X2_FROMINT32X2(vC0);     \
        AE_S32X2F24_##_T##P(vF1, p_x0, _step);  \
        vB2 = AE_SUB32S(vA0, vA2);              \
                \
        vF2 = AE_MOVF24X2_FROMINT32X2(vC2);     \
        AE_S32X2F24_##_T##P(vF2, p_x3, _step);  \
        vB1 = AE_ADD32S(vA1, vA3);              \
                    \
        vF1 = AE_MOVF24X2_FROMINT32X2(vC1);     \
        AE_S32X2F24_##_T##P(vF1, p_x1, _step);  \
        vB3 = AE_SUB32S(vA1, vA3);              \
                 \
        vF2 = AE_MOVF24X2_FROMINT32X2(vC3);     \
        AE_S32X2F24_##_T##P(vF2, p_x2, _step);  \
        vA1 = AE_SUB32S(vB0, vB1);  vA1 = AE_SRAI32R(vA1, 3);             \
      }

#define FFT_BUTTERFLY_S3_0_LAST(_T,_step,_X,twd_step)\
      {                                         \
        vT3 = AE_L32X2F24_I(p_twd, 16);         \
        vC0 = AE_ADD32S(vB0, vB1);    vC0 = AE_SRAI32R(vC0, 3);            \
        vB3 = AE_SEL32_LH(vB3, vB3);            \
        AE_L32X2F24_##_X(vT1, p_twd, twd_step); \
        vA2 = AE_SUBADD32S(vB2, vB3);   vA2 = AE_SRAI32R(vA2, 3);          \
        vA3 = AE_ADDSUB32S(vB2, vB3);   vA3 = AE_SRAI32R(vA3, 3);        \
                   \
        vF1 = AE_MOVF24X2_FROMINT32X2(vA1);     \
        vF0 = AE_MULFC24RA(vF1, vT2);           \
        vC1 = (vF0);                            \
        vF1 = AE_MOVF24X2_FROMINT32X2(vA2);     \
        vF0 = AE_MULFC24RA(vF1, vT3);           \
        vC2 = (vF0);                            \
        vF1 = AE_MOVF24X2_FROMINT32X2(vA3);     \
        vF0 = AE_MULFC24RA(vF1, vT1);           \
        vC3 = (vF0);                            \
        vF1 = AE_MOVF24X2_FROMINT32X2(vC0);     \
        AE_S32X2F24_##_T##P(vF1, p_x0, _step);  \
                 \
        vF2 = AE_MOVF24X2_FROMINT32X2(vC2);     \
        AE_S32X2F24_##_T##P(vF2, p_x3, _step);  \
                \
        vF1 = AE_MOVF24X2_FROMINT32X2(vC1);     \
        AE_S32X2F24_##_T##P(vF1, p_x1, _step);  \
                    \
        vF2 = AE_MOVF24X2_FROMINT32X2(vC3);     \
        AE_S32X2F24_##_T##P(vF2, p_x2, _step);  \
      }

#define FFT_BUTTERFLY_S3_T3F(_T,_step) \
      {                                        \
        vB3 = AE_SEL32_LH(vB3, vB3);           \
        vA1 = AE_SUB32S(vB0, vB1);     vA1 = AE_SRAI32R(vA1, 2);           \
        vC0 = AE_ADD32S(vB0, vB1);             \
        vA2 = AE_SUBADD32S(vB2, vB3);  vA2 = AE_SRAI32R(vA2, 2);          \
        vF1 = AE_L32X2F24_##_T(p_x0, _step);   \
        vA0 = (vF1);                           \
        vA3 = AE_ADDSUB32S(vB2, vB3);  vA3 = AE_SRAI32R(vA3, 2);        \
        vF1 = AE_MOVF24X2_FROMINT32X2(vA1);    \
        vF0 = AE_MULFC24RA(vF1, vT2);          \
        vC1 = (vF0);                           \
        vF2 = AE_L32X2F24_##_T(p_x1, _step);   \
        vA1 = (vF2);                           \
        vF1 = AE_MOVF24X2_FROMINT32X2(vA2);    \
        vF0 = AE_MULFC24RA(vF1, vT3);          \
        vC2 = (vF0);                           \
        vF1 = AE_L32X2F24_##_T(p_x2, _step);   \
        vA2 = (vF1);                           \
        vF1 = AE_MOVF24X2_FROMINT32X2(vA3);    \
        vF0 = AE_MULFC24RA(vF1, vT1);          \
        vC3 = (vF0);                           \
        vF2 = AE_L32X2F24_##_T(p_x3, _step);   \
        vA3 = (vF2);                           \
        /*vC0 = AE_SRAI32(vC0, 2);*/       vC0 = AE_SRAI32R(vC0, 2);        \
        vF1 = AE_MOVF24X2_FROMINT32X2(vC0);    \
        AE_S32X2F24_##_T##P(vF1, p_x0, _step); \
        vB0 = AE_ADD32S(vA0, vA2);             \
        /*vC2 = AE_SRAI32(vC2, 2);  */             \
        vF2 = AE_MOVF24X2_FROMINT32X2(vC2);    \
        AE_S32X2F24_##_T##P(vF2, p_x3, _step); \
        vB2 = AE_SUB32S(vA0, vA2);             \
       /* vC1 = AE_SRAI32(vC1, 2);    */           \
        vF1 = AE_MOVF24X2_FROMINT32X2(vC1);    \
        AE_S32X2F24_##_T##P(vF1, p_x1, _step); \
        vB1 = AE_ADD32S(vA1, vA3);             \
       /* vC3 = AE_SRAI32(vC3, 2);      */         \
        vF2 = AE_MOVF24X2_FROMINT32X2(vC3);    \
        AE_S32X2F24_##_T##P(vF2, p_x2, _step); \
        vB3 = AE_SUB32S(vA1, vA3);             \
      }

#define FFT_BUTTERFLY_S3_T3F_LAST(_T,_step) \
      {                                        \
        vB3 = AE_SEL32_LH(vB3, vB3);           \
        vA1 = AE_SUB32S(vB0, vB1);      vA1 = AE_SRAI32R(vA1, 2);          \
        vC0 = AE_ADD32S(vB0, vB1);        vC0 = AE_SRAI32R(vC0, 2);         \
        vA2 = AE_SUBADD32S(vB2, vB3);    vA2 = AE_SRAI32R(vA2, 2);        \
        vA3 = AE_ADDSUB32S(vB2, vB3);    vA3 = AE_SRAI32R(vA3, 2);      \
        vF1 = AE_MOVF24X2_FROMINT32X2(vA1);    \
        vF0 = AE_MULFC24RA(vF1, vT2);          \
        vC1 = (vF0);                           \
        vF1 = AE_MOVF24X2_FROMINT32X2(vA2);  \
        vF0 = AE_MULFC24RA(vF1, vT3);          \
        vC2 = (vF0);                           \
        vF1 = AE_MOVF24X2_FROMINT32X2(vA3);    \
        vF0 = AE_MULFC24RA(vF1, vT1);          \
        vC3 = (vF0);                           \
        /*vC0 = AE_SRAI32(vC0, 2);*/             \
        vF1 = AE_MOVF24X2_FROMINT32X2(vC0);    \
        AE_S32X2F24_##_T##P(vF1, p_x0, _step); \
       /* vC2 = AE_SRAI32(vC2, 2);*/               \
        vF2 = AE_MOVF24X2_FROMINT32X2(vC2);    \
        AE_S32X2F24_##_T##P(vF2, p_x3, _step); \
        /*vC1 = AE_SRAI32(vC1, 2);    */           \
        vF1 = AE_MOVF24X2_FROMINT32X2(vC1);    \
        AE_S32X2F24_##_T##P(vF1, p_x1, _step); \
        /*vC3 = AE_SRAI32(vC3, 2); */              \
        vF2 = AE_MOVF24X2_FROMINT32X2(vC3);    \
        AE_S32X2F24_##_T##P(vF2, p_x2, _step); \
      }

#define FFT_BUTTERFLY_R2_shift1(idx)          \
    {                                  \
    vA1 = AE_L32X2F24_I(p_x0, 8);      \
    vA2 = AE_L32X2F24_I(p_x0, 16);     \
    vA3 = AE_L32X2F24_I(p_x0, 24);     \
    AE_L32X2F24_IP(vF2, p_x0, 32);     \
    vA0 = (vF2);                       \
    vB0 = AE_ADD32S(vA0, vA1);         \
    vB2 = AE_SUB32S(vA0, vA1);         \
    vB1 = AE_ADD32S(vA2, vA3);         \
    vB3 = AE_SUB32S(vA2, vA3);         \
    vB0 = AE_SRAI32(vB0, 1);            \
    vB2 = AE_SRAI32(vB2, 1);            \
    vB1 = AE_SRAI32(vB1, 1);            \
    vB3 = AE_SRAI32(vB3, 1);            \
    vF1 = AE_MOVF24X2_FROMINT32X2(vB0); \
    AE_S32X2F24_X(vF1, p_y0, idx);     \
    vF2 = AE_MOVF24X2_FROMINT32X2(vB1); \
    AE_S32X2F24_X(vF2, p_y1, idx);     \
    vF1 = AE_MOVF24X2_FROMINT32X2(vB2); \
    AE_S32X2F24_X(vF1, p_y2, idx);     \
    vF2 = AE_MOVF24X2_FROMINT32X2(vB3); \
    AE_S32X2F24_X(vF2, p_y3, idx);     \
    }


static void fft_cplx_24x24_lbut2_s3(int32_t *y, int32_t *x, int N)
{
    ae_int32x2  vA0, vA1, vA2, vA3, vB0, vB1, vB2, vB3;
    ae_f24x2    vF1, vF2;

    int32_t i, i0, i1, i2, i3, ai;
    ae_f24x2 * restrict p_y0 = (ae_f24x2 *)(y);
    ae_f24x2 * restrict p_y1 = (p_y0 + (N >> 2));
    ae_f24x2 * restrict p_y2 = (p_y1 + (N >> 2));
    ae_f24x2 * restrict p_y3 = (p_y2 + (N >> 2));
    ae_f24x2 * restrict p_x0 = (ae_f24x2 *)(x);

    i = NSA(N) + 1;
    ai = ((int32_t)0x1) << i;
    i0 = 0;
    for (i = 0; i < (N >> 4); i++)
    {
        i1 = AE_ADDBRBA32(i0, ai);
        i2 = AE_ADDBRBA32(i1, ai);
        i3 = AE_ADDBRBA32(i2, ai);
        //----------------------------------------------------------------
        // Read eight inputs, and perform radix4 decomposition            
        //----------------------------------------------------------------
        FFT_BUTTERFLY_R2_shift1(i0);
        FFT_BUTTERFLY_R2_shift1(i1);
        FFT_BUTTERFLY_R2_shift1(i2);
        FFT_BUTTERFLY_R2_shift1(i3);
        i0 = AE_ADDBRBA32(i3, ai);
    }

}


/*
   scaled fft with reordering
   N=16 - size of FFT
   NOTE: y is input and output, x - temporary
*/
void fft16_24x24(int32_t *y, int32_t *x, const int32_t *ptwd)
{
        ae_f24x2 * p_x0, * p_x1, * restrict p_y0;
  const ae_f24x2 * restrict p_twd;
  ae_int24x2  vA0, vA1, vA2, vA3, vB0, vB1, vB2, vB3;
  ae_f24x2    vF0, vF1, vF2, vF3, vT1, vT2, vT3;
  ae_f32x2    vFR;

  const int   step_right  = (int)sizeof(ae_f24x2);
  const int   step_left   = -step_right;
  const int   step_down   = 4*step_right;
  const int   step_nextc  = step_right - 3*step_down;
  const int   step_hupleft= 7*step_left;
  const int   step_upleft = 15*step_left;

  p_x0  = (ae_f24x2 *)(x);
  p_x1  = (ae_f24x2 *)(x+2*15);
  p_y0  = (ae_f24x2 *)(y);
  p_twd = (const ae_f24x2 *)(ptwd)+3;

  /*** reordering ***/
  AE_L32X2F24_XP(vF0, p_y0, step_right);
  AE_L32X2F24_XP(vF1, p_y0, step_right);
  AE_L32X2F24_XP(vF2, p_y0, step_right);
  AE_L32X2F24_XP(vF3, p_y0, step_right);
  vA0 = (vF0);
  vA1 = (vF1);
  vA2 = (vF2);
  vA3 = (vF3);

  vB0 = AE_SELP24_HH(vA0, vA1);
  vB3 = AE_SELP24_LL(vA1, vA0);
  AE_L32X2F24_XP(vF0, p_y0, step_right);
  AE_L32X2F24_XP(vF1, p_y0, step_right);
  vA0 = (vF0);
  vA1 = (vF1);
  vB1 = AE_SELP24_HH(vA2, vA3);
  vB2 = AE_SELP24_LL(vA3, vA2);
  AE_L32X2F24_XP(vF2, p_y0, step_right);
  AE_L32X2F24_XP(vF3, p_y0, step_right);
  vA2 = (vF2);
  vA3 = (vF3);
  vF0 = (vB0);
  vF3 = (vB3);
  vF1 = (vB1);
  vF2 = (vB2);
  AE_S32X2F24_XP(vF0, p_x0, step_right);
  AE_S32X2F24_XP(vF3, p_x1, step_left);
  AE_S32X2F24_XP(vF1, p_x0, step_right);
  AE_S32X2F24_XP(vF2, p_x1, step_left);

  vB0 = AE_SELP24_HH(vA0, vA1);
  vB3 = AE_SELP24_LL(vA1, vA0);
  AE_L32X2F24_XP(vF0, p_y0, step_right);
  AE_L32X2F24_XP(vF1, p_y0, step_right);
  vA0 = (vF0);
  vA1 = (vF1);
  vB1 = AE_SELP24_HH(vA2, vA3);
  vB2 = AE_SELP24_LL(vA3, vA2);
  AE_L32X2F24_XP(vF2, p_y0, step_right);
  AE_L32X2F24_XP(vF3, p_y0, step_right);
  vA2 = (vF2);
  vA3 = (vF3);
  vF0 = (vB0);
  vF3 = (vB3);
  vF1 = (vB1);
  vF2 = (vB2);
  AE_S32X2F24_XP(vF0, p_x0, step_right);
  AE_S32X2F24_XP(vF3, p_x1, step_left);
  AE_S32X2F24_XP(vF1, p_x0, step_right);
  AE_S32X2F24_XP(vF2, p_x1, step_left);

  vB0 = AE_SELP24_HH(vA0, vA1);
  vB3 = AE_SELP24_LL(vA1, vA0);
  AE_L32X2F24_XP(vF0, p_y0, step_right);
  AE_L32X2F24_XP(vF1, p_y0, step_right);
  vA0 = (vF0);
  vA1 = (vF1);
  vB1 = AE_SELP24_HH(vA2, vA3);
  vB2 = AE_SELP24_LL(vA3, vA2);
  AE_L32X2F24_XP(vF2, p_y0, step_right);
  AE_L32X2F24_XP(vF3, p_y0, step_upleft);
  vA2 = (vF2);
  vA3 = (vF3);
  vF0 = (vB0);
  vF3 = (vB3);
  vF1 = (vB1);
  vF2 = (vB2);
  AE_S32X2F24_XP(vF0, p_x0, step_right);
  AE_S32X2F24_XP(vF3, p_x1, step_left);
  AE_S32X2F24_XP(vF1, p_x0, step_right);
  AE_S32X2F24_XP(vF2, p_x1, step_left);

  vB0 = AE_SELP24_HH(vA0, vA1);
  vB3 = AE_SELP24_LL(vA1, vA0);
  vB1 = AE_SELP24_HH(vA2, vA3);
  vB2 = AE_SELP24_LL(vA3, vA2);
  vF0 = (vB0);
  vF3 = (vB3);
  vF1 = (vB1);
  vF2 = (vB2);
  AE_S32X2F24_XP(vF0, p_x0, step_right);
  AE_S32X2F24_XP(vF3, p_x1, step_left);
  AE_S32X2F24_XP(vF1, p_x0, step_hupleft);
  AE_S32X2F24_I(vF2, p_x1, 0);

  /*** first radix4 stage ***/
  /* Elements 0, 4, 8, 12 */
  vF0 = AE_L32X2F24_I(p_x1, (0-8)*(int)sizeof(ae_f24x2));
  vF1 = AE_L32X2F24_I(p_x1, (4-8)*(int)sizeof(ae_f24x2));
  vF2 = AE_L32X2F24_I(p_x1, (8-8)*(int)sizeof(ae_f24x2));
  vF3 = AE_L32X2F24_I(p_x1, (12-8)*(int)sizeof(ae_f24x2));
  vA0 = (vF0);
  vA1 = (vF1);
  vA2 = (vF2);
  vA3 = (vF3);

  vA0 = AE_SRAIP24(vA0, 3);
  vA1 = AE_SRAIP24(vA1, 3);
  vA2 = AE_SRAIP24(vA2, 3);
  vA3 = AE_SRAIP24(vA3, 3);

  vB0 = AE_ADDSP24S(vA0, vA2);
  vB2 = AE_SUBSP24S(vA0, vA2);
  vB1 = AE_ADDSP24S(vA1, vA3);
  vA2 = AE_SUBSP24S(vA3, vA1);
  vB3 = AE_SUBSP24S(vA1, vA3);
  vB3 = AE_SELP24_LH(vA2, vB3);

  vA0 = AE_ADDSP24S(vB0, vB1);
  vA1 = AE_SUBSP24S(vB0, vB1);
  vA2 = AE_ADDSP24S(vB2, vB3);
  vA3 = AE_SUBSP24S(vB2, vB3);

  vB0 = AE_SRAIP24(vA0, 2);
  vB1 = AE_SRAIP24(vA1, 2);
  vB2 = AE_SRAIP24(vA2, 2);
  vB3 = AE_SRAIP24(vA3, 2);

  vF0 = (vB0);
  vF1 = (vB1);
  vF2 = (vB2);
  vF3 = (vB3);

  AE_S32X2F24_I(vF0, p_x1, (0-8)*(int)sizeof(ae_f24x2));
  AE_S32X2F24_I(vF3, p_x1, (4-8)*(int)sizeof(ae_f24x2));
  AE_S32X2F24_I(vF1, p_x1, (8-8)*(int)sizeof(ae_f24x2));
  AE_S32X2F24_I(vF2, p_x1, (12-8)*(int)sizeof(ae_f24x2));

  /* Elements 1, 5, 9, 13 */
  vF0 = AE_L32X2F24_I(p_x1, (1-8)*(int)sizeof(ae_f24x2));
  vF1 = AE_L32X2F24_I(p_x1, (5-8)*(int)sizeof(ae_f24x2));
  vF2 = AE_L32X2F24_I(p_x1, (9-8)*(int)sizeof(ae_f24x2));
  vF3 = AE_L32X2F24_I(p_x1, (13-8)*(int)sizeof(ae_f24x2));
  vA0 = (vF0);
  vA1 = (vF1);
  vA2 = (vF2);
  vA3 = (vF3);

  vA0 = AE_SRAIP24(vA0, 3);
  vA1 = AE_SRAIP24(vA1, 3);
  vA2 = AE_SRAIP24(vA2, 3);
  vA3 = AE_SRAIP24(vA3, 3);

  vB0 = AE_ADDSP24S(vA0, vA2);
  vB2 = AE_SUBSP24S(vA0, vA2);
  vB1 = AE_ADDSP24S(vA1, vA3);
  vA2 = AE_SUBSP24S(vA3, vA1);
  vB3 = AE_SUBSP24S(vA1, vA3);
  vB3 = AE_SELP24_LH(vA2, vB3);

  vA0 = AE_ADDSP24S(vB0, vB1);
  vA1 = AE_SUBSP24S(vB0, vB1);
  vA2 = AE_ADDSP24S(vB2, vB3);
  vA3 = AE_SUBSP24S(vB2, vB3);

  AE_L32X2F24_IP(vT1, p_twd, 8);
  AE_L32X2F24_IP(vT2, p_twd, 8);
  AE_L32X2F24_IP(vT3, p_twd, 8);

  vF3 = (vA3);
  vFR = AE_MULFC24RA(vF3, vT1);
  vA3 = AE_MOVINT24X2_FROMF32X2(vFR);

  vF1 = (vA1);
  vFR = AE_MULFC24RA(vF1, vT2);
  vA1 = AE_MOVINT24X2_FROMF32X2(vFR);

  vF2 = (vA2);
  vFR = AE_MULFC24RA(vF2, vT3);
  vA2 = AE_MOVINT24X2_FROMF32X2(vFR);

  vB0 = AE_SRAIP24(vA0, 2);
  vB1 = AE_SRAIP24(vA1, 2);
  vB2 = AE_SRAIP24(vA2, 2);
  vB3 = AE_SRAIP24(vA3, 2);

  vF0 = (vB0);
  vF1 = (vB1);
  vF2 = (vB2);
  vF3 = (vB3);

  AE_S32X2F24_I(vF0, p_x1, (1-8)*(int)sizeof(ae_f24x2));
  AE_S32X2F24_I(vF3, p_x1, (5-8)*(int)sizeof(ae_f24x2));
  AE_S32X2F24_I(vF1, p_x1, (9-8)*(int)sizeof(ae_f24x2));
  AE_S32X2F24_I(vF2, p_x1, (13-8)*(int)sizeof(ae_f24x2));

  /* Elements 2, 6, 10, 14 */
  vF0 = AE_L32X2F24_I(p_x1, (2-8)*(int)sizeof(ae_f24x2));
  vF1 = AE_L32X2F24_I(p_x1, (6-8)*(int)sizeof(ae_f24x2));
  vF2 = AE_L32X2F24_I(p_x1, (10-8)*(int)sizeof(ae_f24x2));
  vF3 = AE_L32X2F24_I(p_x1, (14-8)*(int)sizeof(ae_f24x2));
  vA0 = (vF0);
  vA1 = (vF1);
  vA2 = (vF2);
  vA3 = (vF3);

  vA0 = AE_SRAIP24(vA0, 3);
  vA1 = AE_SRAIP24(vA1, 3);
  vA2 = AE_SRAIP24(vA2, 3);
  vA3 = AE_SRAIP24(vA3, 3);

  vB0 = AE_ADDSP24S(vA0, vA2);
  vB2 = AE_SUBSP24S(vA0, vA2);
  vB1 = AE_ADDSP24S(vA1, vA3);
  vA2 = AE_SUBSP24S(vA3, vA1);
  vB3 = AE_SUBSP24S(vA1, vA3);
  vB3 = AE_SELP24_LH(vA2, vB3);

  vA0 = AE_ADDSP24S(vB0, vB1);
  vA1 = AE_SUBSP24S(vB0, vB1);
  vA2 = AE_ADDSP24S(vB2, vB3);
  vA3 = AE_SUBSP24S(vB2, vB3);

  AE_L32X2F24_IP(vT1, p_twd, 8);
  AE_L32X2F24_IP(vT2, p_twd, 8);
  AE_L32X2F24_IP(vT3, p_twd, 8);

  vF3 = (vA3);
  vFR = AE_MULFC24RA(vF3, vT1);
  vA3 = AE_MOVINT24X2_FROMF32X2(vFR);

  vF1 = (vA1);
  vFR = AE_MULFC24RA(vF1, vT2);
  vA1 = AE_MOVINT24X2_FROMF32X2(vFR);

  vF2 = (vA2);
  vFR = AE_MULFC24RA(vF2, vT3);
  vA2 = AE_MOVINT24X2_FROMF32X2(vFR);

  vB0 = AE_SRAIP24(vA0, 2);
  vB1 = AE_SRAIP24(vA1, 2);
  vB2 = AE_SRAIP24(vA2, 2);
  vB3 = AE_SRAIP24(vA3, 2);

  vF0 = (vB0);
  vF1 = (vB1);
  vF2 = (vB2);
  vF3 = (vB3);

  AE_S32X2F24_I(vF0, p_x1, (2-8)*(int)sizeof(ae_f24x2));
  AE_S32X2F24_I(vF3, p_x1, (6-8)*(int)sizeof(ae_f24x2));
  AE_S32X2F24_I(vF1, p_x1, (10-8)*(int)sizeof(ae_f24x2));
  AE_S32X2F24_I(vF2, p_x1, (14-8)*(int)sizeof(ae_f24x2));

  /* Elements 3, 7, 11, 15 */
  vF0 = AE_L32X2F24_I(p_x1, (3-8)*(int)sizeof(ae_f24x2));
  vF1 = AE_L32X2F24_I(p_x1, (7-8)*(int)sizeof(ae_f24x2));
  vF2 = AE_L32X2F24_I(p_x1, (11-8)*(int)sizeof(ae_f24x2));
  vF3 = AE_L32X2F24_I(p_x1, (15-8)*(int)sizeof(ae_f24x2));
  vA0 = (vF0);
  vA1 = (vF1);
  vA2 = (vF2);
  vA3 = (vF3);

  vA0 = AE_SRAIP24(vA0, 3);
  vA1 = AE_SRAIP24(vA1, 3);
  vA2 = AE_SRAIP24(vA2, 3);
  vA3 = AE_SRAIP24(vA3, 3);

  vB0 = AE_ADDSP24S(vA0, vA2);
  vB2 = AE_SUBSP24S(vA0, vA2);
  vB1 = AE_ADDSP24S(vA1, vA3);
  vA2 = AE_SUBSP24S(vA3, vA1);
  vB3 = AE_SUBSP24S(vA1, vA3);
  vB3 = AE_SELP24_LH(vA2, vB3);

  vA0 = AE_ADDSP24S(vB0, vB1);
  vA1 = AE_SUBSP24S(vB0, vB1);
  vA2 = AE_ADDSP24S(vB2, vB3);
  vA3 = AE_SUBSP24S(vB2, vB3);

  AE_L32X2F24_IP(vT1, p_twd, 8);
  AE_L32X2F24_IP(vT2, p_twd, 8);
  AE_L32X2F24_IP(vT3, p_twd, 8);

  vF3 = (vA3);
  vFR = AE_MULFC24RA(vF3, vT1);
  vA3 = AE_MOVINT24X2_FROMF32X2(vFR);

  vF1 = (vA1);
  vFR = AE_MULFC24RA(vF1, vT2);
  vA1 = AE_MOVINT24X2_FROMF32X2(vFR);

  vF2 = (vA2);
  vFR = AE_MULFC24RA(vF2, vT3);
  vA2 = AE_MOVINT24X2_FROMF32X2(vFR);

  vB0 = AE_SRAIP24(vA0, 2);
  vB1 = AE_SRAIP24(vA1, 2);
  vB2 = AE_SRAIP24(vA2, 2);
  vB3 = AE_SRAIP24(vA3, 2);

  vF0 = (vB0);
  vF1 = (vB1);
  vF2 = (vB2);
  vF3 = (vB3);

  AE_S32X2F24_I(vF0, p_x1, (3-8)*(int)sizeof(ae_f24x2));
  AE_S32X2F24_I(vF3, p_x1, (7-8)*(int)sizeof(ae_f24x2));
  AE_S32X2F24_I(vF1, p_x1, (11-8)*(int)sizeof(ae_f24x2));
  AE_S32X2F24_I(vF2, p_x1, (15-8)*(int)sizeof(ae_f24x2));

  /*** last radix4 stage ***/
  AE_L32X2F24_XP(vF0, p_x0, step_right);
  AE_L32X2F24_XP(vF1, p_x0, step_right);
  AE_L32X2F24_XP(vF2, p_x0, step_right);
  AE_L32X2F24_XP(vF3, p_x0, step_right);
  vA0 = (vF0);
  vA1 = (vF1);
  vA2 = (vF2);
  vA3 = (vF3);

  vB0 = AE_ADDSP24S(vA0, vA2);
  vB2 = AE_SUBSP24S(vA0, vA2);
  vB1 = AE_ADDSP24S(vA1, vA3);
  vA2 = AE_SUBSP24S(vA3, vA1);
  vB3 = AE_SUBSP24S(vA1, vA3);
  vB3 = AE_SELP24_LH(vA2, vB3);

  vA0 = AE_ADDSP24S(vB0, vB1);
  vA2 = AE_SUBSP24S(vB0, vB1);
  vA1 = AE_SUBSP24S(vB2, vB3);
  vA3 = AE_ADDSP24S(vB2, vB3);

  vF0 = (vA0);
  vF1 = (vA1);
  vF2 = (vA2);
  vF3 = (vA3);
  AE_S32X2F24_XP(vF0, p_y0, step_down);
  AE_S32X2F24_XP(vF1, p_y0, step_down);
  AE_S32X2F24_XP(vF2, p_y0, step_down);
  AE_S32X2F24_XP(vF3, p_y0, step_nextc);

  AE_L32X2F24_XP(vF0, p_x0, step_right);
  AE_L32X2F24_XP(vF1, p_x0, step_right);
  AE_L32X2F24_XP(vF2, p_x0, step_right);
  AE_L32X2F24_XP(vF3, p_x0, step_right);
  vA0 = (vF0);
  vA1 = (vF1);
  vA2 = (vF2);
  vA3 = (vF3);

  vB0 = AE_ADDSP24S(vA0, vA2);
  vB2 = AE_SUBSP24S(vA0, vA2);
  vB1 = AE_ADDSP24S(vA1, vA3);
  vA2 = AE_SUBSP24S(vA3, vA1);
  vB3 = AE_SUBSP24S(vA1, vA3);
  vB3 = AE_SELP24_LH(vA2, vB3);

  vA0 = AE_ADDSP24S(vB0, vB1);
  vA2 = AE_SUBSP24S(vB0, vB1);
  vA1 = AE_SUBSP24S(vB2, vB3);
  vA3 = AE_ADDSP24S(vB2, vB3);

  vF0 = (vA0);
  vF1 = (vA1);
  vF2 = (vA2);
  vF3 = (vA3);
  AE_S32X2F24_XP(vF0, p_y0, step_down);
  AE_S32X2F24_XP(vF1, p_y0, step_down);
  AE_S32X2F24_XP(vF2, p_y0, step_down);
  AE_S32X2F24_XP(vF3, p_y0, step_nextc);

  AE_L32X2F24_XP(vF0, p_x0, step_right);
  AE_L32X2F24_XP(vF1, p_x0, step_right);
  AE_L32X2F24_XP(vF2, p_x0, step_right);
  AE_L32X2F24_XP(vF3, p_x0, step_right);
  vA0 = (vF0);
  vA1 = (vF1);
  vA2 = (vF2);
  vA3 = (vF3);

  vB0 = AE_ADDSP24S(vA0, vA2);
  vB2 = AE_SUBSP24S(vA0, vA2);
  vB1 = AE_ADDSP24S(vA1, vA3);
  vA2 = AE_SUBSP24S(vA3, vA1);
  vB3 = AE_SUBSP24S(vA1, vA3);
  vB3 = AE_SELP24_LH(vA2, vB3);

  vA0 = AE_ADDSP24S(vB0, vB1);
  vA2 = AE_SUBSP24S(vB0, vB1);
  vA1 = AE_SUBSP24S(vB2, vB3);
  vA3 = AE_ADDSP24S(vB2, vB3);

  vF0 = (vA0);
  vF1 = (vA1);
  vF2 = (vA2);
  vF3 = (vA3);
  AE_S32X2F24_XP(vF0, p_y0, step_down);
  AE_S32X2F24_XP(vF1, p_y0, step_down);
  AE_S32X2F24_XP(vF2, p_y0, step_down);
  AE_S32X2F24_XP(vF3, p_y0, step_nextc);

  AE_L32X2F24_XP(vF0, p_x0, step_right);
  AE_L32X2F24_XP(vF1, p_x0, step_right);
  AE_L32X2F24_XP(vF2, p_x0, step_right);
  AE_L32X2F24_XP(vF3, p_x0, step_right);
  vA0 = (vF0);
  vA1 = (vF1);
  vA2 = (vF2);
  vA3 = (vF3);

  vB0 = AE_ADDSP24S(vA0, vA2);
  vB2 = AE_SUBSP24S(vA0, vA2);
  vB1 = AE_ADDSP24S(vA1, vA3);
  vA2 = AE_SUBSP24S(vA3, vA1);
  vB3 = AE_SUBSP24S(vA1, vA3);
  vB3 = AE_SELP24_LH(vA2, vB3);

  vA0 = AE_ADDSP24S(vB0, vB1);
  vA2 = AE_SUBSP24S(vB0, vB1);
  vA1 = AE_SUBSP24S(vB2, vB3);
  vA3 = AE_ADDSP24S(vB2, vB3);

  vF0 = (vA0);
  vF1 = (vA1);
  vF2 = (vA2);
  vF3 = (vA3);
  AE_S32X2F24_XP(vF0, p_y0, step_down);
  AE_S32X2F24_XP(vF1, p_y0, step_down);
  AE_S32X2F24_XP(vF2, p_y0, step_down);
  AE_S32X2F24_I(vF3, p_y0, 0);
}

/*
   scaled fft with reordering
   N=32 - size of FFT
   NOTE: y is input and output, x - temporary
*/
void fft32_24x24(int32_t *y, int32_t *x, const int32_t *ptwd)
{
    const int N = 32;
          ae_int32x2 * px0, * px1, * restrict py0;
          ae_f24x2 * restrict p_x0,
                   * restrict p_x1,
                   * restrict p_x2,
                   * restrict p_x3;
    const ae_f24x2 * restrict p_twd;
    ae_int32x2  vA0, vA1, vA2, vA3, vB0, vB1, vB2, vB3;
    ae_int32x2  vC0, vC1, vC2, vC3;
    ae_f24x2    vF1, vF2, vT1, vT2, vT3;
    ae_f24x2    vT1_1, vT2_1, vT3_1;
    ae_f24x2    vT1_2, vT2_2, vT3_2;
    ae_f32x2    vF0;

    int lc, k;
    int step_circ, stride;

    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);

    px0  = (ae_int32x2 *)(x);
    px1  = (ae_int32x2 *)(x+2*N-2);
    py0  = (ae_int32x2 *)(y);

    // reordering
    __Pragma("loop_count min=2");
    for (k=0; k<(N>>2); k++) 
    { 
        AE_L32X2_IP(vA0, py0, sizeof(ae_int32x2));
        AE_L32X2_IP(vA1, py0, sizeof(ae_int32x2));
        AE_L32X2_IP(vA2, py0, sizeof(ae_int32x2));
        AE_L32X2_IP(vA3, py0, sizeof(ae_int32x2));
        vB0 = AE_SEL32_HH(vA0, vA1);
        vB3 = AE_SEL32_LL(vA1, vA0);
        vB1 = AE_SEL32_HH(vA2, vA3);
        vB2 = AE_SEL32_LL(vA3, vA2);
        AE_S32X2_IP(vB0, px0,       sizeof(ae_int32x2));
        AE_S32X2_XP(vB3, px1, -(int)sizeof(ae_int32x2));
        AE_S32X2_IP(vB1, px0,       sizeof(ae_int32x2));
        AE_S32X2_XP(vB2, px1, -(int)sizeof(ae_int32x2));
    }
#if 0
    fft_cplx24x24(y, x, cfft24_32, 3);
#else
    __Pragma("no_reorder");
    p_twd = (const ae_f24x2 *)ptwd;

    WUR_AE_CBEGIN0( (uintptr_t)p_twd );
    WUR_AE_CEND0  ( (uintptr_t)(p_twd)+6*N );
    step_circ=24;
    stride = N; // The stride is quartered with every iteration of the outer loop.

    { // first iteration for pre-scaling

      // Set up pointers to access "N/4", "N/2", "3N/4" complex point or  
      // "N/2", "N", "3N/2" half word                                    
      //-----------------------------------------------------------------
      p_x0 = (ae_f24x2 *)x;
      p_x1 = p_x0 + (stride>>2);
      p_x2 = p_x1 + (stride>>2);
      p_x3 = p_x2 + (stride>>2);

      vA0 = AE_L32X2F24_I(p_x0, 0);
      vA1 = AE_L32X2F24_I(p_x1, 0);
      vA2 = AE_L32X2F24_I(p_x2, 0);
      vA3 = AE_L32X2F24_I(p_x3, 0);
      /*
      vA0 = AE_SRAI32(vA0, 3);
      vA1 = AE_SRAI32(vA1, 3);
      vA2 = AE_SRAI32(vA2, 3);
      vA3 = AE_SRAI32(vA3, 3);
      */
      vB0 = AE_ADD32S(vA0, vA2);
      vB2 = AE_SUB32S(vA0, vA2);
      vB1 = AE_ADD32S(vA1, vA3);
      vB3 = AE_SUB32S(vA1, vA3);

      vT2 = AE_L32X2F24_I(p_twd, 8);

      vA1 = AE_SUB32S(vB0, vB1);
      vA1 = AE_SRAI32R(vA1, 3);
      for(lc=0; lc<((N>>2)-1); lc++)
      {
        FFT_BUTTERFLY_S3_0(I, 2*4,XC,step_circ)
      };
      FFT_BUTTERFLY_S3_0_LAST(I, 2*4,XC,step_circ)

      stride >>= 2;
      step_circ<<=2;
    }

    {
      //-----------------------------------------------------------------
      // Set up pointers to access "N/4", "N/2", "3N/4" complex point or  
      // "N/2", "N", "3N/2" half word                                    
      //-----------------------------------------------------------------
      p_x0 = (ae_f24x2 *)x;
      p_x1 = p_x0 + (stride>>2);
      p_x2 = p_x1 + (stride>>2);
      p_x3 = p_x2 + (stride>>2);

      vA0 = AE_L32X2F24_I(p_x0, 0);
      vA1 = AE_L32X2F24_I(p_x1, 0);
      vA2 = AE_L32X2F24_I(p_x2, 0);
      vA3 = AE_L32X2F24_I(p_x3, 0);

      vB0 = AE_ADD32S(vA0, vA2);
      vB2 = AE_SUB32S(vA0, vA2);
      vB1 = AE_ADD32S(vA1, vA3);
      vB3 = AE_SUB32S(vA1, vA3);

      vT2_1 = AE_L32X2F24_I(p_twd, 8);          
      vT3_1 = AE_L32X2F24_I(p_twd, 16);         
      AE_L32X2F24_XC(vT1_1, p_twd,step_circ); 
      vT2_2 = AE_L32X2F24_I(p_twd, 8);          
      vT3_2 = AE_L32X2F24_I(p_twd, 16);         
      AE_L32X2F24_XC(vT1_2, p_twd,step_circ); 

      for(lc=0; lc<((N>>4)-1); lc++) 
      {
        vT1 = vT1_1; vT2 = vT2_1; vT3 = vT3_1;
        FFT_BUTTERFLY_S3_T3F(I, 2*4)  
        vT1 = vT1_2; vT2 = vT2_2; vT3 = vT3_2;
        FFT_BUTTERFLY_S3_T3F(I, 14*4)      
        
        vT1 = vT1_1; vT2 = vT2_1; vT3 = vT3_1;
        FFT_BUTTERFLY_S3_T3F(I, 2*4)  
        vT1 = vT1_2; vT2 = vT2_2; vT3 = vT3_2;
        FFT_BUTTERFLY_S3_T3F(I, 14*4)      
      };
      {
        vT1 = vT1_1; vT2 = vT2_1; vT3 = vT3_1;
        FFT_BUTTERFLY_S3_T3F(I, 2*4)  
        vT1 = vT1_2; vT2 = vT2_2; vT3 = vT3_2;
        FFT_BUTTERFLY_S3_T3F(I, 14*4)      
  
        vT1 = vT1_1; vT2 = vT2_1; vT3 = vT3_1;
        FFT_BUTTERFLY_S3_T3F(I, 2*4)
        vT1 = vT1_2; vT2 = vT2_2; vT3 = vT3_2;
        FFT_BUTTERFLY_S3_T3F_LAST(I, 14*4)
      }

      stride >>= 2;
    }

    fft_cplx_24x24_lbut2_s3(y, x, N);
#endif
}
