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
Reference C code
Integrit, 2006-2017
*/
#include "NatureDSP_Signal_fft.h"
#include "NatureDSP_Signal_vector.h"
#include "common.h"

#define SWAP_PTR(_x, _y) {complex_fract16 *tmp = _x; _x = _y ; _y = tmp; } 
#define DFT4XI2(x0, x1, x2, x3)/* output x0, x3, x1, x2*/\
{\
    ae_int16x4 t1, t2, t3; \
    AE_ADDANDSUBRNG16RAS_S1(x0, x2); \
    AE_ADDANDSUBRNG16RAS_S1(x1, x3); \
    \
    x3 = AE_MOVINT16X4_FROMF16X4(AE_MUL16JS(AE_MOVF16X4_FROMINT16X4(x3))); \
    \
    AE_ADDANDSUBRNG16RAS_S2(x0, x1); \
    AE_ADDANDSUBRNG16RAS_S2(x2, x3); \
    \
    t1 = x3; \
    t2 = x1; \
    t3 = x2; \
    x1 = t1; \
    x2 = t2; \
    x3 = t3; \
}

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
#if  (XCHAL_HAVE_HIFI3Z)
inline_ void MULCx2(ae_int16x4 *x, const ae_int16x4 *t)
{
    ae_f16x4 f;
    f = AE_MULFC16RAS(AE_MOVF16X4_FROMINT16X4(*x), AE_MOVF16X4_FROMINT16X4(*t));
    *x = AE_MOVINT16X4_FROMF16X4(f);
}

/*
First stage radix 4
*/
static int stage_first_DFT4_16x16_ie(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    const int R = 4; // stage radix
    int i;
    const int stride = N / R;
    ae_int16x4 * restrict _py;
    ae_int16x4 * restrict _px;
    ae_int32  * restrict ptw1;
    ae_int32  * restrict ptw2;
    ae_int32  * restrict ptw3;

    int shift, tw_inc, s;
    const int min_shift = 3;

    if (*bexp < min_shift)
    {
        shift = min_shift - *bexp;
    }
    else
    {
        shift = 0;
    }

    switch (shift)
    {
    case 0: WUR_AE_SAR(0);   break;
    case 1: WUR_AE_SAR(2);   break;
    case 2: WUR_AE_SAR(4);   break;
    case 3: WUR_AE_SAR(5);   break;
    }

    tw_inc = tw_step;
    s = *bexp;

    shift = 3 - s;
    _py = (ae_int16x4 *)y;
    _px = (ae_int16x4 *)x;

    ptw1 = (ae_int32  *)tw;
    ptw2 = (ae_int32  *)(N / 4 * tw_step * sizeof(complex_fract16)+(uintptr_t)tw);
    ptw3 = (ae_int32  *)(2 * N / 4 * tw_step * sizeof(complex_fract16)+(uintptr_t)tw);

    WUR_AE_SAR(5); // Set scaling = 3 
    __Pragma("loop_count min=2 factor=2");
    // First stage: 13 cycles per pipeline stage in steady state with unroll=1
    for (i = 0; i < N / R / 2; i++)
    {
        ae_int16x4 _t1, _t2, _t3;
        ae_int32x2 t10, t20, t30;
        ae_int32x2 _t11, _t21, _t31;
        ae_int16x4 _x0, _x1, _x2, _x3;

        AE_L32_XP(t10, ptw1, tw_inc * sizeof(complex_fract16));
        AE_L32_XP(t20, ptw2, tw_inc * sizeof(complex_fract16));
        AE_L32_XP(t30, ptw3, tw_inc * sizeof(complex_fract16));

        AE_L32_XP(_t11, ptw1, tw_inc * sizeof(complex_fract16));
        AE_L32_XP(_t21, ptw2, tw_inc * sizeof(complex_fract16));
        AE_L32_XP(_t31, ptw3, tw_inc * sizeof(complex_fract16));

        t10 = AE_SEL32_LL(t10, _t11);
        t20 = AE_SEL32_LL(t20, _t21);
        t30 = AE_SEL32_LL(t30, _t31);

        AE_L16X4_XP(_x0, _px, stride*sizeof(complex_fract16));
        AE_L16X4_XP(_x1, _px, stride*sizeof(complex_fract16));
        AE_L16X4_XP(_x2, _px, stride*sizeof(complex_fract16));
        AE_L16X4_XP(_x3, _px, sizeof(_x3)-3 * stride*sizeof(complex_fract16));

        // Normalize input data
        _x0 = AE_SLAA16S(_x0, s);
        _x1 = AE_SLAA16S(_x1, s);
        _x2 = AE_SLAA16S(_x2, s);
        _x3 = AE_SLAA16S(_x3, s);

        DFT4XI2(_x0, _x1, _x2, _x3);

        _t1 = AE_MOVINT16X4_FROMF32X2(t10);
        _t2 = AE_MOVINT16X4_FROMF32X2(t20);
        _t3 = AE_MOVINT16X4_FROMF32X2(t30);

        MULCx2(&_x1, &_t1);
        MULCx2(&_x2, &_t2);
        MULCx2(&_x3, &_t3);

        AE_S16X4RNG_IP(AE_SEL16_7632(_x0, _x1), _py, sizeof(*_py));
        AE_S16X4RNG_IP(AE_SEL16_7632(_x2, _x3), _py, sizeof(*_py));
        AE_S16X4RNG_IP(AE_SEL16_5410(_x0, _x1), _py, sizeof(*_py));
        AE_S16X4RNG_IP(AE_SEL16_5410(_x2, _x3), _py, sizeof(*_py));
    } //for (i = 0; i < N / R / _v; i++)

    {
        /* Bits 5, 4, and 3 respectively of AE_SAR are set if bits 14, 13 and 12
        of any quarter of d is different than their respective sign bit.  */
        int sar = RUR_AE_SAR();
        *bexp = NSA((sar << (30 - 5)) | 1);
    }

    *v *= R;
    return shift;
} //stage_first_DFT4_16x16_ie
#else //#if  (XCHAL_HAVE_HIFI3Z)
static int stage_first_DFT4_16x16_ie(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    const int R = 4; // stage radix
    const int stride = N / R;
    const int _v = v[0];
    int i, s;
    ae_int16x4 * restrict _py;
    ae_int16x4 * restrict _px;
    ae_int32x2  t10, t20, t30, t11, t21, t31;
    ae_int32x2  x10, x20, x30, x11, x21, x31;

    ae_p16x2s * restrict  ptw16x2;
    ae_int16x4 acc16 = AE_MOVINT16X4_FROMF32X2(AE_MOVI(0));
    const int tw_inc0 = N / 4 * tw_step/_v * sizeof(complex_fract16);
    int shift;

    s = *bexp;

    ptw16x2 = (ae_p16x2s *)(/*-tw_inc +*/(uintptr_t)tw);

    shift = 3 - s;
    _py = (ae_int16x4 *)y;
    _px = (ae_int16x4 *)x;


    WUR_AE_SAR(5); // Set scaling = 3 
    __Pragma("loop_count min=2 factor=2");
    // First stage: 11 cycles per pipeline stage in steady state with unroll=1
    for (i = 0; i < N / R / 2; i++)
    {
        ae_int16x4 _x0, _x1, _x2, _x3;

        t10 = AE_L16X2M_X(ptw16x2, 0        );
        t20 = AE_L16X2M_X(ptw16x2, tw_inc0  );
        t30 = AE_L16X2M_X(ptw16x2, 2*tw_inc0);
        ptw16x2 += tw_step; 
        t11 = AE_L16X2M_X(ptw16x2, 0        );
        t21 = AE_L16X2M_X(ptw16x2, tw_inc0  );
        t31 = AE_L16X2M_X(ptw16x2, 2*tw_inc0);
        ptw16x2 += tw_step;

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

        AE_L16X4_XP(_x0, _px, stride*sizeof(complex_fract16));
        AE_L16X4_XP(_x1, _px, stride*sizeof(complex_fract16));
        AE_L16X4_XP(_x2, _px, stride*sizeof(complex_fract16));
        AE_L16X4_XP(_x3, _px, sizeof(_x3)-3 * stride*sizeof(complex_fract16));

        DFT4XI2_HIFI3(_x0, _x1, _x2, _x3, shift);

        x10 = AE_MULFC32X16RAS_H(t10, _x1);
        x11 = AE_MULFC32X16RAS_L(t11, _x1);
        x20 = AE_MULFC32X16RAS_H(t20, _x2);
        x21 = AE_MULFC32X16RAS_L(t21, _x2);
        x30 = AE_MULFC32X16RAS_H(t30, _x3);
        x31 = AE_MULFC32X16RAS_L(t31, _x3);

        _x1 = AE_ROUND16X4F32SASYM(x10, x11);
        _x2 = AE_ROUND16X4F32SASYM(x20, x21);
        _x3 = AE_ROUND16X4F32SASYM(x30, x31);

        AE_S16X4_IP(AE_SEL16_7632(_x0, _x1), _py, sizeof(*_py));
        AE_S16X4_IP(AE_SEL16_7632(_x2, _x3), _py, sizeof(*_py));
        AE_S16X4_IP(AE_SEL16_5410(_x0, _x1), _py, sizeof(*_py));
        AE_S16X4_IP(AE_SEL16_5410(_x2, _x3), _py, sizeof(*_py));

        _x0 = AE_ABS16S(_x0);
        _x1 = AE_ABS16S(_x1);
        _x2 = AE_ABS16S(_x2);
        _x3 = AE_ABS16S(_x3);

        acc16 = AE_OR16(acc16, _x0);
        acc16 = AE_OR16(acc16, _x1);
        acc16 = AE_OR16(acc16, _x2);
        acc16 = AE_OR16(acc16, _x3);

    } //for (i = 0; i < N / R / _v; i++)
    ae_int16x4 tmp0 = AE_SEL16_5432(acc16, acc16);  //AE_INTSWAP
    acc16 = AE_OR16(acc16, tmp0);
    acc16 = AE_OR16(acc16, AE_SHORTSWAP(acc16));
    int  tmpi = AE_MOVAD16_0(acc16) >> (15 - 3);
    WUR_AE_SAR(tmpi << 3);

    {
        /* Bits 5, 4, and 3 respectively of AE_SAR are set if bits 14, 13 and 12
        of any quarter of d is different than their respective sign bit.  */
        int sar = RUR_AE_SAR();
        *bexp = NSA((sar << (30 - 5)) | 1);
    }

    *v *= R;  
    return shift;
} //stage_first_DFT4_16x16_ie, hifi3
#endif //#if  (XCHAL_HAVE_HIFI3Z)


#if  (XCHAL_HAVE_HIFI3Z)
static int stage_last_DFT4_16x16_ie(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    int i;
    const int R = 4; // stage radix
    const int min_shift = 2;
    const int stride = N / R;
    ae_int16x4 * restrict _py;
    ae_int16x4 * restrict _px;
    int shift;

    if (*bexp < min_shift)
    {
        shift = min_shift - *bexp;
    }
    else
    {
        shift = 0;
    }

    switch (shift)
    {
        case 0: WUR_AE_SAR(0);   break;
        case 1: WUR_AE_SAR(2);   break;
        case 2: WUR_AE_SAR(4);   break;
        case 3: WUR_AE_SAR(5);   break;
    }
    
    _py = (ae_int16x4 *)y;
    _px = (ae_int16x4 *)x;

    // Last phase, without twiddles
    __Pragma("loop_count min=2 factor=2");
    for (i = 0; i < N / 8; i++)
    { /* 10 cycles per pipeline stage in steady state with unroll=2 */
        ae_int16x4 _x0, _x1, _x2, _x3;
        AE_L16X4_XP(_x0, _px, stride*sizeof(complex_fract16));
        AE_L16X4_XP(_x1, _px, stride*sizeof(complex_fract16));
        AE_L16X4_XP(_x2, _px, stride*sizeof(complex_fract16));
        AE_L16X4_XP(_x3, _px, sizeof(_x3)-3 * stride*sizeof(complex_fract16));

        DFT4XI2(_x0, _x1, _x2, _x3);

        AE_S16X4RNG_XP(_x0, _py, stride*sizeof(complex_fract16));
        AE_S16X4RNG_XP(_x1, _py, stride*sizeof(complex_fract16));
        AE_S16X4RNG_XP(_x2, _py, stride*sizeof(complex_fract16));
        AE_S16X4RNG_XP(_x3, _py, sizeof(_x3)-3 * stride*sizeof(complex_fract16));
    }
    return shift ;
 } //stage_last_DFT4_16x16_ie, hifi3z
#else //#if  (XCHAL_HAVE_HIFI3Z)

static int stage_last_DFT4_16x16_ie(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    const int R = 4; // stage radix
    const int stride = N / R;
    const int _v = v[0];
    const int min_shift = 2;
    int i;
    ae_int16x4 * restrict _py;
    ae_int16x4 * restrict _px;
    int shift;

    if (*bexp < min_shift)
    {
        shift = min_shift - *bexp;
    }
    else
    {
        shift = 0;
    }

    if (_v == stride)
    {
        _py = (ae_int16x4 *)y;
        _px = (ae_int16x4 *)x;
        // Last phase, without twiddles
        __Pragma("loop_count min=2 factor=2");
        for (i = 0; i < N / 8; i++)
        {
            ae_int16x4 _x0, _x1, _x2, _x3;
            AE_L16X4_XP(_x0, _px, stride*sizeof(complex_fract16));
            AE_L16X4_XP(_x1, _px, stride*sizeof(complex_fract16));
            AE_L16X4_XP(_x2, _px, stride*sizeof(complex_fract16));
            AE_L16X4_XP(_x3, _px, sizeof(_x3)-3 * stride*sizeof(complex_fract16));

            DFT4XI2_HIFI3(_x0, _x1, _x2, _x3, shift);

            AE_S16X4_XP(_x0, _py, stride*sizeof(complex_fract16));
            AE_S16X4_XP(_x1, _py, stride*sizeof(complex_fract16));
            AE_S16X4_XP(_x2, _py, stride*sizeof(complex_fract16));
            AE_S16X4_XP(_x3, _py, sizeof(_x3)-3 * stride*sizeof(complex_fract16));
        }
    }
    return shift ;
} //stage_last_DFT4_16x16_ie, hifi3z

#endif //#if  (XCHAL_HAVE_HIFI3Z)

extern int stage_inner_DFT4_16x16_ie(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp);


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
int fft_cplx16x16_ie(complex_fract16* y, complex_fract16* x, const complex_fract16* twd, int twdstep, int N, int scalingOpt)
{
    int bexp, shift = 0;
    int v = 1;
    complex_fract16 *pdest = y;
    int log2N = 30 - NSA(N);
    ae_int16x4 * restrict px;
    ae_int16x4 * restrict py;

    NASSERT_ALIGN8(twd);
    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT(x != y);
    NASSERT(scalingOpt == 2);
    NASSERT(N == 128 || N == 256 || N == 512 || N == 1024);
    {
        int i;
        ae_int16x4 acc = AE_MOVINT16X4_FROMINT32X2(AE_MOVI(0)), tmp;

        __Pragma("loop_count min=4 factor=4");
        px = (ae_int16x4*)x;
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
           
    shift += stage_first_DFT4_16x16_ie((const int16_t*)twd, (int16_t*)x, (int16_t*)y, N, &v, twdstep, &bexp);
    SWAP_PTR(x, y);
    log2N -= 2;
    twdstep *= 4;

    while (log2N > 2)
    {
        shift += stage_inner_DFT4_16x16_ie((const int16_t*)twd, (int16_t*)x, (int16_t*)y, N, &v, twdstep, &bexp);
        SWAP_PTR(x, y);
        log2N -= 2;
        twdstep *= 4;
    }

    if (y != pdest)
    {
        /* Execute the last stage inplace */
        y = x;
    }

    /* Last stage */      
    if (log2N & 1)
    {
     //   shift += fft_stageS2_DFT2_16x16(NULL, (int16_t*)x, (int16_t*)y, N, &v, 0, &bexp);
        const int stride = N / 2;
        int i; 
        px = (ae_int16x4 *)x;
        py = (ae_int16x4 *)y;
        if (bexp == 0)
        {
            shift++;    
             __Pragma("loop_count min=2 factor=2");
            for (i = 0; i < N / 4; i++)
            {
                ae_int16x4 x0, x1, tmp0;
                /* 3 cycles per pipeline stage in steady state with unroll=1 */
                AE_L16X4_XP(x0, px, stride*sizeof(complex_fract16));
                AE_L16X4_XP(x1, px, sizeof(*px) - stride*sizeof(complex_fract16));

                x0 = AE_SRAI16(x0, 1);
                x1 = AE_SRAI16(x1, 1);

                tmp0 = AE_ADD16S(x0, x1);
                x1 = AE_SUB16S(x0, x1); //AE_ADDANDSUBRNG16RAS_S1(x0, x1); 

                AE_S16X4_XP(tmp0, py, stride*sizeof(complex_fract16));
                AE_S16X4_XP(x1, py, sizeof(*px) - stride*sizeof(complex_fract16));
            }
        }
        else
        {
            __Pragma("loop_count min=2 factor=2");
            for (i = 0; i < N / 4; i++)
            {
                ae_int16x4 x0, x1, tmp0;
                /* 4 cycles per pipeline stage in steady state with unroll=2*/
                AE_L16X4_XP(x0, px, stride*sizeof(complex_fract16));
                AE_L16X4_XP(x1, px, sizeof(*px) - stride*sizeof(complex_fract16));

                tmp0 = AE_ADD16S(x0, x1);
                x1 = AE_SUB16S(x0, x1); 

                AE_S16X4_XP(tmp0, py, stride*sizeof(complex_fract16));
                AE_S16X4_XP(x1, py, sizeof(*px) - stride*sizeof(complex_fract16));
            }
        }
    }
    else
    {
        shift += stage_last_DFT4_16x16_ie(NULL, (int16_t*)x, (int16_t*)y, N, &v, 0, &bexp);
    }
    return shift;
} /* fft_cplx16x16_ie*/ 


