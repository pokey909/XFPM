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
#include "NatureDSP_Signal_vector.h"
#include "common.h"


#if (XCHAL_HAVE_HIFI3Z)

inline_ void MULCx2(ae_int16x4 *x, const ae_int16x4 *t)
{
    ae_f16x4 f;
    f = AE_MULFC16RAS(AE_MOVF16X4_FROMINT16X4(*x), AE_MOVF16X4_FROMINT16X4(*t));
    *x = AE_MOVINT16X4_FROMF16X4(f);
}




/*
    In-place inverse split part of FFT:
    x[N+2]  input (N+2 samples)/output(N samples)
    N       size of FFT
    Internal scaling is 2-lshift 

*/
static void isplitPart16x16_ie_lshift(complex_fract16 *_X, int N, const complex_fract16 *tw, int twstep, int lshift)
{
    int i, step;
    int16_t *x = (int16_t*)_X;
    ae_p16x2s * restrict p_x0, *restrict p_x1;

    ae_int32x2  vA0, vA1, vB0, vB1, vC0, vR;


    NASSERT_ALIGN8(_X);

    p_x0 = (ae_p16x2s *)x;
    p_x1 = (ae_p16x2s *)(x + N);

    vR = AE_MOVI(1);
    vR = AE_SLAI32(vR, 8);

    // first point
    vA0 = AE_L16X2M_I(p_x0, 0);
    vA1 = AE_L16X2M_I(p_x1, 0);
    vA0 = AE_SLAA32(vA0, lshift-1);
    vA1 = AE_SLAA32(vA1, lshift-1);

    vB0 = AE_ADD32S(vA0, vA1);
    vB1 = AE_SUB32S(vA0, vA1);
    vB0 = AE_SEL32_HH(vB0, vB1);
    vA1 = AE_MOVI(1);
    vA1 = AE_SLAI32(vA1, 8);
    vB1 = AE_MOVI(0);
    vB0 = AE_ADD32(vB0, vA1);
    vB0 = AE_SRAI32(vB0, 1);
    AE_S16X2M_I(vB0, p_x0, 0);
    AE_S16X2M_I(vB1, p_x1, 0);

    AE_L16X2M_IU(vA0, p_x0, 4);
    AE_L16X2M_IU(vA1, p_x1, -4);
    vA0 = AE_SLAA32(vA0, lshift - 1);
    vA1 = AE_SLAA32(vA1, lshift - 1);

    ae_int16x4 * restrict px0 = (ae_int16x4 *)&x[2];
    ae_int16x4 * restrict px1 = (ae_int16x4 *)&x[N - 2 * (1 + 1)];
    ae_int16x4 * restrict py0 = (ae_int16x4 *)&x[2];
    ae_int16x4 * restrict py1 = (ae_int16x4 *)&x[N - 2 * (1 + 1)];
    ae_valign v0 = AE_LA64_PP(px0);
    ae_valign vy0 = AE_ZALIGN64();
    xtbool4 mov_even = (int)0x5;

    ae_int32x2 * restrict pw = (ae_int32x2*)(twstep + tw);
    step = twstep*sizeof(complex_fract16); 

    // set shift = 2 for AE_ADDANDSUBRNG16RAS_S2 and shift=0 for AE_ADDANDSUBRNG16RAS_S0
    WUR_AE_SAR(4);

    i = (N - 2) >> 3;
    /*  8 cycles per pipeline stage in steady state with unroll=1 */
    do//  for (i = 1; i < (N >> 2)-1; i+=2)
    {
        ae_int16x4 s, d, t0, t1, tw;
        ae_int32x2 tmp0, tmp1;

        AE_LA16X4_IP(s, v0, px0);
        AE_L16X4_XP(d, px1, -(int)sizeof(*px1));
        s = AE_SLAA16S(s, lshift); 
        d = AE_SLAA16S(d, lshift);

        d = AE_SEL16_5432(d, d);

        AE_ADDANDSUBRNG16RAS_S2(s, d);

        /*
        t0 = {s1_im, d1_re, s0_im, d0_re }
        t1 = {d1_im, s1_re, d0_im, s0_re }
        */
        t0 = s;
        t1 = d;
        AE_MOVT16X4(t0, d, mov_even);
        AE_MOVT16X4(t1, s, mov_even);
        s = t0;
        d = t1;
        /*
        x = xre + j * xim
        (-j*x)' = xim + j*xre
        -(-j*(-j*x)' * tw)' = -( (-j*x)' * (-j*tw)  )' = (j*x) *(-j*tw)'
        */
        /* Re and Im part are swapped */
        AE_L32_XP(tmp0, castxcc(const ae_int32, pw), step);
        AE_L32_XP(tmp1, castxcc(const ae_int32, pw), step);
        /* Pack 2 twiddles */
        tw = AE_SEL16_5432(AE_MOVINT16X4_FROMINT32X2(tmp0),
            AE_MOVINT16X4_FROMINT32X2(tmp1));
        tw = AE_CONJ16S(tw); 

        d = AE_MUL16JS(d);
        MULCx2(&d, &tw);

        AE_ADDANDSUBRNG16RAS_S1(s, d);

        d = AE_MUL16JS(d);
        d = AE_SHORTSWAP(d);

        AE_SA16X4_IP(s, vy0, py0);
        AE_S16X4_XP(d, py1, -8);
    } while (--i);

    AE_SA64POS_FP(vy0, py0);

    {
        ae_int16x4 s, d, t0, t1, tw;
        ae_int32x2 tmp0, tmp1;

        AE_LA16X4_IP(s, v0, px0);
        AE_L16X4_XP(d, px1, -(int)sizeof(*px1));
        s = AE_SLAA16S(s, lshift);
        d = AE_SLAA16S(d, lshift);

        d = AE_SEL16_5432(d, d);

        AE_ADDANDSUBRNG16RAS_S2(s, d);

        /*
        t0 = {s1_im, d1_re, s0_im, d0_re }
        t1 = {d1_im, s1_re, d0_im, s0_re }
        */
        t0 = s;
        t1 = d;
        AE_MOVT16X4(t0, d, mov_even);
        AE_MOVT16X4(t1, s, mov_even);
        s = t0;
        d = t1;

        /* Re and Im part are swapped */
        AE_L32_XP(tmp0, castxcc(const ae_int32, pw), step);
        AE_L32_XP(tmp1, castxcc(const ae_int32, pw), step);
        /* Pack 2 twiddles */
        tw = AE_SEL16_5432(AE_MOVINT16X4_FROMINT32X2(tmp0),
            AE_MOVINT16X4_FROMINT32X2(tmp1));
        tw = AE_CONJ16S(tw);

        d = AE_MUL16JS(d);
        MULCx2(&d, &tw);

        AE_ADDANDSUBRNG16RAS_S1(s, d);

        d = AE_MUL16JS(d);
        d = AE_SHORTSWAP(d);

        ((int16_t*)py0)[0]/*x[N/2-2]*/ = AE_MOVAD16_3(s);
        ((int16_t*)py0)[1]/*x[N / 2 - 1]*/ = AE_MOVAD16_2(s);

        ((int16_t*)py1)[2]/*x[N/2+2]*/ = AE_MOVAD16_1(d);
        ((int16_t*)py1)[3]/*x[N / 2 + 3]*/ = AE_MOVAD16_0(d);

    }
    // middle sample
    vA0 = AE_L16X2M_I((ae_p16x2s *)(x + N / 2), 0);
    vA0 = AE_SLAA32(vA0, lshift - 1);

    vB0 = AE_NEG32S(vA0);
    vC0 = AE_SEL32_HL(vA0, vB0);
    AE_S16X2M_I(vC0, (ae_p16x2s *)(x + N / 2), 0);
}
#else //#if (XCHAL_HAVE_HIFI3Z)
#define _CONJ32(_x) {_x = AE_SEL32_HL(_x, AE_NEG32S(_x) ); }
#define _INVRE32(_x) {_x = AE_SEL32_HL( AE_NEG32S(_x), _x ); }

/*
    In-place inverse split part of FFT:
    x[N+2]  input (N+2 samples)/output(N samples)
    N       size of FFT
    Internal scaling is 2-lshift
*/
static void isplitPart16x16_ie_lshift(complex_fract16 *_X, int N, const complex_fract16 *tw, int twstep, int lshift)
{
    int i, step;
    int16_t *X = (int16_t*)_X; 

    ae_p16x2s * restrict p_x0, *restrict p_x1;
    const ae_p16x2s * restrict p_twd;

    ae_int32x2  vA0, vA1, vB0, vB1, vC0, vC1, vR;
    ae_int16x4  vT;
    ae_f32x2    vF0, vF1;
    ae_f16x4    vF2;
    NASSERT_ALIGN8(X);

    step = twstep*sizeof(complex_fract16);
    p_twd = (const ae_p16x2s *)tw;

    p_x0 = (ae_p16x2s *)X;
    p_x1 = (ae_p16x2s *)(X + N);

    vR = AE_MOVI(1);
    vR = AE_SLAI32(vR, 9);

    // first point
    vA0 = AE_L16X2M_I(p_x0, 0);
    vA1 = AE_L16X2M_I(p_x1, 0);

    vA0 = AE_SLAA32(vA0, lshift); 
    vA1 = AE_SLAA32(vA1, lshift);

    vB0 = AE_ADD32S(vA0, vA1);
    vB1 = AE_SUB32S(vA0, vA1);
    vB0 = AE_SEL32_HH(vB0, vB1);
    vA1 = AE_MOVI(1);
    vA1 = AE_SLAI32(vA1, 9);
    vB1 = AE_MOVI(0);
    vB0 = AE_ADD32(vB0, vA1);
    vB0 = AE_SRAI32(vB0, 2);
    AE_S16X2M_I(vB0, p_x0, 0);
    AE_S16X2M_I(vB1, p_x1, 0);

    AE_L16X2M_IU(vA0, p_x0, 4);
    AE_L16X2M_IU(vA1, p_x1, -4);

    vA0 = AE_SLAA32(vA0, lshift);
    vA1 = AE_SLAA32(vA1, lshift);

    for (i = 1; i < (N >> 2); i++)
    {
        // load twiddle
        AE_L16X2M_XU(vB1, p_twd, step);
        _INVRE32(vB1);
        vB1 = AE_SRAI32(vB1, 8);
        vT = AE_CVT16X4(vB1, vB1);

        // ADD/SUBB
        vB0 = AE_ADD32S(vA0, vA1);
        vB1 = AE_SUB32S(vA0, vA1);

        vA0 = AE_SEL32_LH(vB0, vB1);
        vB1 = AE_SEL32_HL(vB0, vB1);

        // do rotation
        vF1 = (vA0);
        vF2 = (vT);
        vF0 = AE_MULFC32X16RAS_H(vF1, vF2);
        vB0 = (vF0);
        vA1 = AE_NEG32S(vB0);
        vB0 = AE_SEL32_LH(vA1, vA1);

        // load next data
        AE_L16X2M_IU(vA0, p_x0, 4);
        AE_L16X2M_IU(vA1, p_x1, -4);
        vA0 = AE_SLAA32(vA0, lshift);
        vA1 = AE_SLAA32(vA1, lshift);
        // SUM/DIFF
        vC0 = AE_ADD32S(vB1, vB0);
        vC1 = AE_SUB32S(vB1, vB0);
        vB1 = AE_NEG32S(vC1);
        vC1 = AE_SEL32_HL(vC1, vB1);

        vC0 = AE_ADD32S(vC0, vR);
        vC1 = AE_ADD32S(vC1, vR);

        vC0 = AE_SRAI32(vC0, 2);
        vC1 = AE_SRAI32(vC1, 2);

        AE_S16X2M_I(vC0, p_x0, -4);
        AE_S16X2M_I(vC1, p_x1, 4);
    }

    // middle sample
    vB0 = AE_NEG32S(vA0);
    vC0 = AE_SEL32_HL(vA0, vB0);
    vC0 = AE_SRAI32(vC0, 1);
    AE_S16X2M_I(vC0, p_x0, 0);
} //isplitPart_x2
#endif //#if (XCHAL_HAVE_HIFI3Z)

/*-------------------------------------------------------------------------
  Inverse FFT on Real Data with Optimized Memory Usage
  These functions make inverse FFT on real data from half of spectrum with
  optimized memory usage.
  Scaling: 
      +-----------------------+--------------------------------------+
      |      Function         |           Scaling options            |
      +-----------------------+--------------------------------------+
      | ifft_real16x16_ie     |  2 - 16-bit dynamic scaling          |
      | ifft_real32x16_ie     |  3 - fixed scaling before each stage |
      | ifft_real24x24_ie     |  3 - fixed scaling before each stage |
      | ifft_real24x24_ie_24p |  3 - fixed scaling before each stage |
      | ifft_real32x16_ie_24p |  1 - 24-bit scaling                  |
      | ifft_real32x32_ie     |  2 - 32-bit dynamic scaling          |    
      |                       |  3 - fixed scaling before each stage |   
      +-----------------------+--------------------------------------+
  NOTES:
  1. Bit-reversing reordering is done here.
  2. INPUT DATA MAY APPEAR DAMAGED after the call.
  3. FFT functions may use input and output buffers for temporal storage
     of intermediate 32-bit data, so FFT functions with 24-bit packed
     I/O (Nx3-byte data) require that the buffers are large enough to 
     keep Nx4-byte data.
  4. FFT of size N may be supplied with constant data (twiddle factors) 
     of a larger-sized FFT = N*twdstep.

  Precision:
  16x16_ie      16-bit input/outputs, 16-bit data, 16-bit twiddles
  24x24_ie      24-bit input/outputs, 24-bit data, 24-bit twiddles
  32x16_ie      32-bit input/outputs, 32-bit data, 16-bit twiddles
  32x32_ie      32-bit input/outputs, 32-bit data, 32-bit twiddles
  24x24_ie_24p  24-bit packed input/outputs, 24-bit data, 24-bit twiddles
  32x16_ie_24p  24-bit packed input/outputs, 32-bit data, 16-bit twiddles
  f_ie          floating point

  Input:
  x                     input spectrum (positive side). Real and imaginary
                        data are interleaved and real data goes first:
  --------------+----------+-----------------+----------------
  Function      |   Size   |  Allocated Size |       type    |
  --------------+----------+-----------------+----------------
  16x16_ie      |   N/2+1  |      N/2+1      |complex_fract16|
  24x24_ie      |   N/2+1  |      N/2+1      |complex_fract32|
  32x16_ie      |   N/2+1  |      N/2+1      |complex_fract32|
  32x32_ie      |   N/2+1  |      N/2+1      |complex_fract32|
  24x24_ie_24p  |   3*(N+2)|      4*N+8      |       uint8_t |
  32x16_ie_24p  |   3*(N+2)|      4*N+8      |       uint8_t |
  f_ie          |   N/2+1  |      N/2+1      | complex_float |
  --------------+----------+-----------------+----------------

  twd[2*N*twdstep*3/4]  twiddle factor table of a complex-valued FFT
                        of size N*twdstep
  N                     FFT size
  twdstep               twiddle step
  scalingOpt            scaling option (see table above), not applicable 
                        to the floating point function
  Output:
  y                     output spectrum. Real and imaginary data are 
                        interleaved and real data goes first:
  --------------+----------+-----------------+-----------
  Function      |   Size   |  Allocated Size |  type    |
  --------------+----------+-----------------+-----------
  16x16_ie      |     N    |      N          |  int16_t |
  24x24_ie      |     N    |      N          |   f24    |
  32x16_ie      |     N    |      N          |  int32_t |
  32x32_ie      |     N    |      N          |  int32_t |
  24x24_ie_24p  |    3*N   |      4*N+8      |  uint8_t |
  32x16_ie_24p  |    3*N   |      4*N+8      |  uint8_t |
  f_ie          |      N   |      N          | float32_t|
  --------------+----------+-----------------+-----------

  Returned value: total number of right shifts occurred during scaling
  procedure

  Restrictions:
  x,y   should not overlap
  x,y   aligned on 8-bytes boundary
-------------------------------------------------------------------------*/

int ifft_real16x16_ie(int16_t* y, complex_fract16* x, const complex_fract16* twd, int twdstep, int N, int scalingOpt)
{

    int scale, bexp; 
    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT_ALIGN8(twd);
    NASSERT(x != (complex_fract16*)y);
    NASSERT(N == 256 || N == 512 || N == 1024);
    NASSERT(scalingOpt == 2);

    bexp = vec_bexp16((int16_t*)x, N+2) - 16; 
    scale = 2-bexp; 
    isplitPart16x16_ie_lshift(x, N, twd, twdstep, bexp);

    scale += ifft_cplx16x16_ie((complex_fract16*)y, x, twd, twdstep*2, N/2, scalingOpt);
    
    return scale;
}
