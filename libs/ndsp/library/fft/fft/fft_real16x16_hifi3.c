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
    Real/halfcomplex FFT 32x16 with scaling option 3
	Integrit, 2006-2017
*/

#include "NatureDSP_Signal_fft.h"
#include "fft_real_twiddles.h"
#include "fft_cplx_twiddles.h"


#if (XCHAL_HAVE_HIFI3Z)

inline_ void MULCx2(ae_int16x4 *x, const ae_int16x4 *t)
{
    ae_f16x4 f;
    f = AE_MULFC16RAS(AE_MOVF16X4_FROMINT16X4(*x), AE_MOVF16X4_FROMINT16X4(*t));
    *x = AE_MOVINT16X4_FROMF16X4(f);
}

static void splitPart_x2(int16_t * x, 
                         int N, 
                         int shift /*Additional scaling  0  - used with ScalingOpt=3, 
                                                         1  - used with ScalingOpt=2 */)
{
    int i, step;
    ae_int32x2  vA0, vA1, vB0, vB1, vC0;
    xtbool4 mov_even = (int)0x5;
    ae_int16x4 * restrict px0 = (ae_int16x4 *)&x[2];
    ae_int16x4 * restrict px1 = (ae_int16x4 *)&x[N - 2 * (1 + 1)];
    ae_int16x4 * restrict py0 = (ae_int16x4 *)&x[2];
    ae_int16x4 * restrict py1 = (ae_int16x4 *)&x[N - 2 * (1 + 1)];
    ae_int32x2 * restrict pw;
    ae_valign v0 = AE_LA64_PP(px0);
    ae_valign vy0 = AE_ZALIGN64();

    NASSERT_ALIGN8(x);
    ASSERT(shift==0 || shift==1); 

    WUR_AE_SAR(2|shift); // set shift = 1 for AE_ADDANDSUBRNG16RAS_S2 and scaling 0 or 1 for AE_ADDANDSUBRNG16RAS_S1
    // split part: remainder
    {
        ae_q56s tmp = AE_CVTQ48A32S(N);
        step = 1 << (AE_NSAQ56S(tmp) - (38 - MAX_RFFT_PWR));
        step *= sizeof(int32_t);
    }
    pw = (ae_int32x2*)(step + (uintptr_t)twiddleSplit);
    // load data and prepare pointers for pre-increment
    // first and last samples
    vA0 = AE_L16X2M_I((ae_p16x2s *)x, 0);
    vA0 = AE_SRAA32(vA0, shift);

    vA1 = AE_SEL32_LH(vA0, vA0);
    vB0 = AE_ADD32S(vA0, vA1);
    vB1 = AE_SUB32S(vA0, vA1);
    vA1 = AE_MOVI(0);
    vB0 = AE_SEL32_HH(vB0, vA1);
    vB1 = AE_SEL32_HH(vB1, vA1);
    AE_S16X2M_I(vB0, (ae_p16x2s *)x, 0);
    AE_S16X2M_I(vB1, (ae_p16x2s *)(x + N), 0);


     /* 7 cycles per pipeline stage in steady state with unroll=1 */
    i = (N - 2) >> 3; // for (i = 1; i < (N >> 2) - 1; i += 2)
    do
    {   

        ae_int16x4 x0, x1, s, d, t0, t1, tw;
        ae_int32x2 tmp0, tmp1;
        AE_LA16X4_IP(x0, v0, px0);
        AE_L16X4_XP(x1, px1, -(int)sizeof(*px1));
        x1 = AE_SEL16_5432(x1, x1);

        s = x0; 
        d = x1;
        AE_ADDANDSUBRNG16RAS_S1(s, d);  

        t0 = s;
        AE_MOVT16X4(t0, d, mov_even);
        t1 = d;
        AE_MOVT16X4(t1, s, mov_even);

        AE_L32_XP(tmp0, castxcc(const ae_int32, pw), step);
        AE_L32_XP(tmp1, castxcc(const ae_int32, pw), step);

        tw = AE_SEL16_5432(
            AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(tmp0)), AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(tmp1)));
        MULCx2(&t1, &tw);

        s = t0; 
        d = t1; 
        AE_ADDANDSUBRNG16RAS_S2(s, d); 

        t1 = AE_MUL16JS(s);
        t1 = AE_SHORTSWAP(t1);
        AE_S16X4_XP(t1, py1, -8);
        AE_SA16X4_IP(d, vy0, py0);
    } while (--i);

    AE_SA64POS_FP(vy0, py0);
    {

        ae_int16x4 x0, x1, s, d, t0, t1, tw;
        ae_int32x2 tmp0, tmp1;

        AE_LA16X4_IP(x0, v0, px0);
        AE_L16X4_XP(x1, px1, -(int)sizeof(*px1));
        x1 = AE_SEL16_5432(x1, x1);

        s = x0;
        d = x1;
        AE_ADDANDSUBRNG16RAS_S1(s, d); 

        t0 = s;
        AE_MOVT16X4(t0, d, mov_even);
        t1 = d;
        AE_MOVT16X4(t1, s, mov_even);

        AE_L32_XP(tmp0, castxcc(const ae_int32, pw), step);
        AE_L32_XP(tmp1, castxcc(const ae_int32, pw), step);

        tw = AE_SEL16_5432(
            AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(tmp0)), AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(tmp1)));
        MULCx2(&t1, &tw);

        t0 = AE_SRAI16(t0, 1);
        t1 = AE_SRAI16(t1, 1);
        s = AE_ADD16S(t0, t1);
        d = AE_SUB16S(t0, t1);

        ((int16_t*)py0)[0]/*x[N/2-2]*/ = AE_MOVAD16_3(d);
        ((int16_t*)py0)[1]/*x[N / 2 - 1]*/ = AE_MOVAD16_2(d);

        ((int16_t*)py1)[2]/*x[N/2+2]*/ = AE_MOVAD16_3(s);
        s = AE_NEG16S(s);
        ((int16_t*)py1)[3]/*x[N / 2 + 3]*/ = AE_MOVAD16_2(s);
    }
    // middle sample
    vA0 = AE_L16X2M_I((ae_p16x2s*)&x[N / 2], 0); //x[N/2]
    vA0 = AE_SRAA32(vA0, shift);
    vB0 = AE_NEG32S(vA0);
    vC0 = AE_SEL32_HL(vA0, vB0);
    AE_S16X2M_I(vC0, (ae_p16x2s*)&x[N / 2], 0);
}
#else //#if (XCHAL_HAVE_HIFI3Z)
static void splitPart_x2(int16_t * x, int N)
{
  int i, step;

  ae_p16x2s * restrict p_x0, * restrict p_x1;
  const ae_p16x2s * restrict p_twd;

  ae_int32x2  vA0, vA1, vB0, vB1, vC0, vC1, vR;
  ae_int16x4  vT;
  ae_f32x2    vF0, vF1;
  ae_f16x4    vF2;

  NASSERT_ALIGN8(x);

  // split part: remainder
  {
    ae_q56s tmp = AE_CVTQ48A32S(N);
    step = 1<<(AE_NSAQ56S(tmp)-(38-MAX_RFFT_PWR));
    step *= sizeof(int32_t);
  }

  p_twd = (const ae_p16x2s *)twiddleSplit;
  p_x0 = (ae_p16x2s *)x;
  p_x1 = (ae_p16x2s *)(x+N);

  // load data and prepare pointers for pre-increment
  // first and last samples
  vA0 = AE_L16X2M_I(p_x0, 0);
  vA1 = AE_SEL32_LH(vA0, vA0);
  vB0 = AE_ADD32S(vA0, vA1);
  vB1 = AE_SUB32S(vA0, vA1);
  vA1 = AE_MOVI(0);
  vB0 = AE_SEL32_HH(vB0, vA1);
  vB1 = AE_SEL32_HH(vB1, vA1);
  AE_S16X2M_I(vB0, p_x0, 0);
  AE_S16X2M_I(vB1, p_x1, 0);

  AE_L16X2M_IU(vA0, p_x0, 4);
  AE_L16X2M_IU(vA1, p_x1, -4);

  vR = AE_MOVI(1);
  vR = AE_SLAI32(vR, 8);
  /*   Hifi3>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
    15 cycles per pipeline stage in steady state with unroll=2
   
*/
  for (i = 1; i < N>>2; i++) 
  {
    // load twiddle
    AE_L16X2M_XU(vB1, p_twd, step);
    vB1 = AE_SRAI32(vB1, 8);
    vT = AE_CVT16X4(vB1, vB1);

    // ADD/SUBB
    vB0 = AE_ADD32S(vA0, vA1);
    vB1 = AE_SUB32S(vA0, vA1);

    vA0 = AE_SEL32_HL(vB1, vB0);
    vB1 = AE_SEL32_HL(vB0, vB1);

    // do rotation
    vF1 = (vA0);
    vF2 = (vT);
    vF0 = AE_MULFC32X16RAS_H(vF1, vF2);
    vB0 = (vF0);

    // load next data
    AE_L16X2M_IU(vA0, p_x0, 4);
    AE_L16X2M_IU(vA1, p_x1, -4);

    // SUM/DIFF
    vC0 = AE_SUB32S(vB1, vB0);
    vC1 = AE_ADD32S(vB1, vB0);
    vB1 = AE_NEG32S(vC1);
    vC1 = AE_SEL32_HL(vC1, vB1);

    vC0 = AE_ADD32S(vC0, vR);
    vC1 = AE_ADD32S(vC1, vR);
    vC0 = AE_SRAI32(vC0, 1);
    vC1 = AE_SRAI32(vC1, 1);

    AE_S16X2M_I(vC0, p_x0, -4);
    AE_S16X2M_I(vC1, p_x1, 4);
  } 

  // middle sample
  vB0 = AE_NEG32S(vA0);
  vC0 = AE_SEL32_HL(vA0, vB0);
  AE_S16X2M_I(vC0, p_x0, 0);
}

/*
    Split converter with additional scaling = 1
*/
static void splitPart_x2_shift1(int16_t * x, int N)
{
    int i, step;

    ae_p16x2s * restrict p_x0, *restrict p_x1;
    const ae_p16x2s * restrict p_twd;

    ae_int32x2  vA0, vA1, vB0, vB1, vC0, vC1, vR;
    ae_int16x4  vT;
    ae_f32x2    vF0, vF1;
    ae_f16x4    vF2;

    NASSERT_ALIGN8(x);

    // split part: remainder
    {
        ae_q56s tmp = AE_CVTQ48A32S(N);
        step = 1 << (AE_NSAQ56S(tmp) - (38 - MAX_RFFT_PWR));
        step *= sizeof(int32_t);
    }

    p_twd = (const ae_p16x2s *)twiddleSplit;
    p_x0 = (ae_p16x2s *)x;
    p_x1 = (ae_p16x2s *)(x + N);

    // load data and prepare pointers for pre-increment
    // first and last samples
    vA0 = AE_L16X2M_I(p_x0, 0);
    vA0 = AE_SRAI32(vA0, 1);

    vA1 = AE_SEL32_LH(vA0, vA0);
    vB0 = AE_ADD32S(vA0, vA1);
    vB1 = AE_SUB32S(vA0, vA1);


    vA1 = AE_MOVI(0);
    vB0 = AE_SEL32_HH(vB0, vA1);
    vB1 = AE_SEL32_HH(vB1, vA1);
    AE_S16X2M_I(vB0, p_x0, 0);
    AE_S16X2M_I(vB1, p_x1, 0);

    AE_L16X2M_IU(vA0, p_x0, 4);
    AE_L16X2M_IU(vA1, p_x1, -4);
    vA0 = AE_SRAI32(vA0, 1);
    vA1 = AE_SRAI32(vA1, 1);

    vR = AE_MOVI(1);
    vR = AE_SLAI32(vR, 8);
    /*   Hifi3>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    15 cycles per pipeline stage in steady state with unroll=2

    */
    for (i = 1; i < N >> 2; i++)
    {
        // load twiddle
        AE_L16X2M_XU(vB1, p_twd, step);
        vB1 = AE_SRAI32(vB1, 8);
        vT = AE_CVT16X4(vB1, vB1);

        // ADD/SUBB
        vB0 = AE_ADD32S(vA0, vA1);
        vB1 = AE_SUB32S(vA0, vA1);

        vA0 = AE_SEL32_HL(vB1, vB0);
        vB1 = AE_SEL32_HL(vB0, vB1);

        // do rotation
        vF1 = (vA0);
        vF2 = (vT);
        vF0 = AE_MULFC32X16RAS_H(vF1, vF2);
        vB0 = (vF0);

        // load next data
        AE_L16X2M_IU(vA0, p_x0, 4);
        AE_L16X2M_IU(vA1, p_x1, -4);
        vA0 = AE_SRAI32(vA0, 1);
        vA1 = AE_SRAI32(vA1, 1);
        // SUM/DIFF
        vC0 = AE_SUB32S(vB1, vB0);
        vC1 = AE_ADD32S(vB1, vB0);
        vB1 = AE_NEG32S(vC1);
        vC1 = AE_SEL32_HL(vC1, vB1);

        vC0 = AE_ADD32S(vC0, vR);
        vC1 = AE_ADD32S(vC1, vR);
        vC0 = AE_SRAI32(vC0, 1);
        vC1 = AE_SRAI32(vC1, 1);

        AE_S16X2M_I(vC0, p_x0, -4);
        AE_S16X2M_I(vC1, p_x1, 4);
    }

    // middle sample
    vB0 = AE_NEG32S(vA0);
   
    vC0 = AE_SEL32_HL(vA0, vB0);    
    AE_S16X2M_I(vC0, p_x0, 0);
} //splitPart_x2_shift1
#endif //#if (XCHAL_HAVE_HIFI3Z)

/*-------------------------------------------------------------------------
  FFT on Real Data
  These functions make FFT on real data forming half of spectrum
      Scaling: 
      +-------------------+----------------------------------------+
      |      Function     |           Scaling options              |
      +-------------------+----------------------------------------+
      |  fft_real16x16    |  2 - 16-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  fft_real32x32    |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  fft_real32x16    |  3 - fixed scaling before each stage   | 
      |  fft_real24x24    |  0 - no scaling                        | 
      |                   |  1 - 24-bit scaling                    |
      |                   |  2 - 32-bit scaling on the first stage |
      |                   |  and 24-bit scaling later              |
      |                   |  3 - fixed scaling before each stage   |
      +-------------------+----------------------------------------+
  NOTES:
  1. Bit-reversal reordering is done here. 
  2. FFT runs in-place so INPUT DATA WILL APPEAR DAMAGED after the call.
  3. Real data FFT function calls fft_cplx() to apply complex FFT of size
     N/2 to input data and then transforms the resulting spectrum.
  4. 32x32 FFTs support mixed radix transforms 
  5. N - FFT size

  Precision:
  32x32  32-bit input/outputs, 32-bit twiddles
  24x24  24-bit input/outputs, 24-bit twiddles
  32x16  32-bit input/outputs, 16-bit twiddles
  16x16  16-bit input/outputs, 16-bit twiddles

  Input:
  x[N]          input signal
  scalingOpt    scaling option (see table above)
  Output:
  y[(N/2+1)*2]  output spectrum (positive side)

  Returned value: total number of right shifts occurred during scaling 
                  procedure

  Restrictions:
  x,y           should not overlap
  x,y           aligned on a 8-bytes boundary
-------------------------------------------------------------------------*/
int fft_real16x16( int16_t* y,int16_t* x,fft_handle_t h,int scalingOpt)
{
  int scale;
  int N;

  NASSERT_ALIGN8(x);
  NASSERT_ALIGN8(y);
  NASSERT(scalingOpt == 3 || scalingOpt == 2);

  N=(((const tFftDescr*)h)->N)<<1;

  scale = fft_cplx16x16(y, x, h, scalingOpt);    /* fft of half-size              */

#if (XCHAL_HAVE_HIFI3Z)
    if (scalingOpt == 3)
    {
        splitPart_x2(y, N, 0);
    }
    else
    {
        splitPart_x2(y, N, 1);
        scale += 1;
    }
#else //#if (XCHAL_HAVE_HIFI3Z)
    if (scalingOpt == 3)
    {
        splitPart_x2(y, N);
    }
    else
    {
        splitPart_x2_shift1(y, N);
        scale += 1;
    }

#endif

  return scale;
}
