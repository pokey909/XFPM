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
Internal stages of fft_cplx16x16_ie, ifft_cplx16x16_ie.
*/

int stage_inner_DFT4_16x16_ie(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    const int R = 4; // stage radix
    int i;
    const int stride = N / R;
    ae_int16x4 * restrict _py;
    ae_int16x4 * restrict _px;
    ae_int32  * restrict ptw1;

    const int _v = v[0];


    int shift;
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

    ptw1 = (ae_int32  *)tw;


    ae_int16x4 _t1, _t2, _t3;
    ae_int32x2 t10, t20, t30;

    _py = (ae_int16x4 *)y;
    _px = (ae_int16x4 *)x;


    ASSERT(shift >= 0 && shift <= 3);

#if 1
    //  Two loops
    int j;
     
    if (_v == 4)
    {
        int tw_inc = N / 4 * tw_step / _v * sizeof(complex_fract16); 

        __Pragma("loop_count min=2 factor=2");
        for (i = 0; i < N / R / _v; i++)
        {
           

            AE_L32_XP(t10, ptw1, tw_inc);
            AE_L32_XP(t20, ptw1, tw_inc);
            AE_L32_XP(t30, ptw1, tw_step * sizeof(complex_fract16) - 2 * tw_inc);

            _t1 = AE_MOVINT16X4_FROMF32X2(t10);
            _t2 = AE_MOVINT16X4_FROMF32X2(t20);
            _t3 = AE_MOVINT16X4_FROMF32X2(t30);
            /* 14 cycles per pipeline stage in steady state with unroll=2 ?!! */
            ASSERT(_v >= 2);
           
            //for (j = 0; j < _v; j += 2)
            {
                ae_int16x4 _x0, _x1, _x2, _x3;
                AE_L16X4_XP(_x0, _px, stride*sizeof(complex_fract16));
                AE_L16X4_XP(_x1, _px, stride*sizeof(complex_fract16));
                AE_L16X4_XP(_x2, _px, stride*sizeof(complex_fract16));
                AE_L16X4_XP(_x3, _px, sizeof(_x3)-3 * stride*sizeof(complex_fract16));

                DFT4XI2(_x0, _x1, _x2, _x3);

                MULCx2(&_x1, &_t1);
                MULCx2(&_x2, &_t2);
                MULCx2(&_x3, &_t3);

                AE_S16X4RNG_XP(_x0, _py, _v * sizeof(complex_fract16));
                AE_S16X4RNG_XP(_x1, _py, _v * sizeof(complex_fract16));
                AE_S16X4RNG_XP(_x2, _py, _v * sizeof(complex_fract16));
                AE_S16X4RNG_XP(_x3, _py, sizeof(ae_int16x4) - (3 * _v) * sizeof(complex_fract16) );
               
            } 
            {
                ae_int16x4 _x0, _x1, _x2, _x3;
                AE_L16X4_XP(_x0, _px, stride*sizeof(complex_fract16));
                AE_L16X4_XP(_x1, _px, stride*sizeof(complex_fract16));
                AE_L16X4_XP(_x2, _px, stride*sizeof(complex_fract16));
                AE_L16X4_XP(_x3, _px, sizeof(_x3)-3 * stride*sizeof(complex_fract16));

                DFT4XI2(_x0, _x1, _x2, _x3);

                MULCx2(&_x1, &_t1);
                MULCx2(&_x2, &_t2);
                MULCx2(&_x3, &_t3);

                AE_S16X4RNG_XP(_x0, _py, _v * sizeof(complex_fract16));
                AE_S16X4RNG_XP(_x1, _py, _v * sizeof(complex_fract16));
                AE_S16X4RNG_XP(_x2, _py, _v * sizeof(complex_fract16));
                AE_S16X4RNG_XP(_x3, _py, sizeof(ae_int16x4) );
            }
        } //for (i = 0; i < N / R / _v; i++)
    }
    else //if _v==4
    {
        int tw_inc = N / 4 * tw_step / _v * sizeof(complex_fract16);
        __Pragma("loop_count min=2 factor=2");
        for (i = 0; i < N / R / _v; i++)
        {/* ~ 17 cycles in outer loop */


            AE_L32_XP(t10, ptw1, tw_inc);
            AE_L32_XP(t20, ptw1, tw_inc);
            AE_L32_XP(t30, ptw1, tw_step * sizeof(complex_fract16)-2 * tw_inc);

            _t1 = AE_MOVINT16X4_FROMF32X2(t10);
            _t2 = AE_MOVINT16X4_FROMF32X2(t20);
            _t3 = AE_MOVINT16X4_FROMF32X2(t30);
            
            ASSERT(_v >= 2);
            __Pragma("loop_count min=2 factor=2");
            for (j = 0; j < _v; j += 2)
            {  /* 7 cycles per pipeline stage in steady state with unroll=1 */
                ae_int16x4 _x0, _x1, _x2, _x3;
                AE_L16X4_XP(_x0, _px, stride*sizeof(complex_fract16));
                AE_L16X4_XP(_x1, _px, stride*sizeof(complex_fract16));
                AE_L16X4_XP(_x2, _px, stride*sizeof(complex_fract16));
                AE_L16X4_XP(_x3, _px, sizeof(_x3)-3 * stride*sizeof(complex_fract16));

                DFT4XI2(_x0, _x1, _x2, _x3);

                MULCx2(&_x1, &_t1);
                MULCx2(&_x2, &_t2);
                MULCx2(&_x3, &_t3);

                AE_S16X4RNG_XP(_x0, _py, _v * sizeof(complex_fract16));
                AE_S16X4RNG_XP(_x1, _py, _v * sizeof(complex_fract16));
                AE_S16X4RNG_XP(_x2, _py, _v * sizeof(complex_fract16));
                AE_S16X4RNG_XP(_x3, _py, sizeof(_x3) - (3 * _v) * sizeof(complex_fract16));                
            } //for (j = 0; j < _v; j++)
            _py = (ae_int16x4*)(3 *_v * sizeof(complex_fract16) + (uintptr_t)_py);
        } //for (i = 0; i < N / R / _v; i++)
    }
#else
    uint32_t  flag = 0;
    uint32_t  dflag = (_v > 1) ? 0x80000000 >> (30 - NSA(_v) - 2) :
                                 0x00000000;
    const int tw_inc0 = N / 4 * tw_step / _v * sizeof(complex_fract16);

    ae_int32x2 start_incs = AE_MOVDA32X2((sizeof(ae_int16x4)-3 * v[0] * sizeof(complex_fract16)),
        -2 * tw_inc0);


    _py = (ae_int16x4*)y;
    ASSERT(_v >= 2);
    //  Merged loops
    //   10 cycles per pipeline stage in steady state with unroll=1
    __Pragma("loop_count min=2 factor=2");
    for (i = 0; i < N / 8; i++)
    {
        ae_int16x4 _x0, _x1, _x2, _x3;

        int py_inc;// = (sizeof(ae_int16x4) - 3 * v[0] * sizeof(complex_fract16));
        int tw_inc;// = -2*tw_inc0;

        py_inc = AE_MOVAD32_H(start_incs);
        tw_inc = AE_MOVAD32_L(start_incs);

        flag += dflag;
        /*
        if (flag == 0)
        {
        ASSERT(j == (_v/2 - 1));
        py_inc = sizeof(ae_int16x4)* (3 * v[0] / 2 + 1);
        tw_inc = tw_step * sizeof(complex_fract16);
        }
        */
        XT_MOVEQZ(py_inc, sizeof(ae_int16x4), flag);
        XT_MOVEQZ(tw_inc, tw_step * sizeof(complex_fract16)-2 * tw_inc0, flag);


        AE_L32_XP(t10, ptw1, tw_inc0);
        AE_L32_XP(t20, ptw1, tw_inc0);
        AE_L32_XP(t30, ptw1, tw_inc);

        _t1 = AE_MOVINT16X4_FROMF32X2(t10);
        _t2 = AE_MOVINT16X4_FROMF32X2(t20);
        _t3 = AE_MOVINT16X4_FROMF32X2(t30);

        AE_L16X4_XP(_x0, _px, stride*sizeof(complex_fract16));
        AE_L16X4_XP(_x1, _px, stride*sizeof(complex_fract16));
        AE_L16X4_XP(_x2, _px, stride*sizeof(complex_fract16));
        AE_L16X4_XP(_x3, _px, sizeof(_x3)-3 * stride*sizeof(complex_fract16));

        DFT4XI2(_x0, _x1, _x2, _x3);

        MULCx2(&_x1, &_t1);
        MULCx2(&_x2, &_t2);
        MULCx2(&_x3, &_t3);

        AE_S16X4RNG_XP(_x0, _py, _v * sizeof(complex_fract16));
        AE_S16X4RNG_XP(_x1, _py, _v * sizeof(complex_fract16));
        AE_S16X4RNG_XP(_x2, _py, _v * sizeof(complex_fract16));
        AE_S16X4RNG_XP(_x3, _py, py_inc);
    } //for (i = 0; i < N / R / _v; i++)
#endif
    {
        /* Bits 5, 4, and 3 respectively of AE_SAR are set if bits 14, 13 and 12
        of any quarter of d is different than their respective sign bit.  */
        int sar = RUR_AE_SAR();
        *bexp = NSA((sar << (30 - 5)) | 1);
    }

    *v *= R;
    return shift;
} //stage_inner_DFT4_16x16_ie, hifi3z
#else //#if (XCHAL_HAVE_HIFI3Z)

/*
    Internal stages of fft_cplx16x16_ie, ifft_cplx16x16_ie.
*/
int stage_inner_DFT4_16x16_ie(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    const int R = 4; // stage radix
    const int min_shift = 3;
    int i;
    const int stride = N / R;
    ae_int16x4 * restrict _py;
    ae_int16x4 * restrict _px;
    ae_p16x2s * restrict  ptw16x2;
    const int _v = v[0];
    ae_int16x4 acc16 = AE_MOVINT16X4_FROMF32X2(AE_MOVI(0));
    const int tw_inc0 = N / 4 * tw_step / _v * sizeof(complex_fract16);
    int shift;

    if (*bexp < min_shift)
    {
        shift = min_shift - *bexp;
    }
    else
    {
        shift = 0;
    }

    ptw16x2 = (ae_p16x2s *)(/*-tw_inc +*/(uintptr_t)tw);

    _py = (ae_int16x4 *)y;
    _px = (ae_int16x4 *)x;

    ASSERT(shift >= 0 && shift <= 3);

    //  Two loops
    int j;
    __Pragma("loop_count min=2");
    for (i = 0; i < N / R / _v; i++)
    {
        ae_int32x2    t10, t20, t30;
        _py = (ae_int16x4*)(4 * _v * i  * sizeof(complex_fract16)+(uintptr_t)y);

        t10 = AE_L16X2M_X(ptw16x2, 0);
        t20 = AE_L16X2M_X(ptw16x2, tw_inc0);
        t30 = AE_L16X2M_X(ptw16x2, 2 * tw_inc0);
        ptw16x2 += tw_step;

        t10 = AE_SLAI32(t10, 8);
        t20 = AE_SLAI32(t20, 8);
        t30 = AE_SLAI32(t30, 8);

        t10 = AE_SEL32_LH(t10, t10);
        t20 = AE_SEL32_LH(t20, t20);
        t30 = AE_SEL32_LH(t30, t30);

        /* 14 cycles per pipeline stage in steady state with unroll=2 ?!! */
        ASSERT(_v >= 2);
        __Pragma("loop_count min=2 factor=2");
        for (j = 0; j < _v; j += 2)
        {
            ae_int16x4 _x0, _x1, _x2, _x3;
            ae_int32x2  x10, x20, x30, x11, x21, x31;

            AE_L16X4_XP(_x0, _px, stride*sizeof(complex_fract16));
            AE_L16X4_XP(_x1, _px, stride*sizeof(complex_fract16));
            AE_L16X4_XP(_x2, _px, stride*sizeof(complex_fract16));
            AE_L16X4_XP(_x3, _px, sizeof(_x3)-3 * stride*sizeof(complex_fract16));

            DFT4XI2_HIFI3(_x0, _x1, _x2, _x3, shift);

            x10 = AE_MULFC32X16RAS_H(t10, _x1);
            x11 = AE_MULFC32X16RAS_L(t10, _x1);
            x20 = AE_MULFC32X16RAS_H(t20, _x2);
            x21 = AE_MULFC32X16RAS_L(t20, _x2);
            x30 = AE_MULFC32X16RAS_H(t30, _x3);
            x31 = AE_MULFC32X16RAS_L(t30, _x3);

            _x1 = AE_ROUND16X4F32SASYM(x10, x11);
            _x2 = AE_ROUND16X4F32SASYM(x20, x21);
            _x3 = AE_ROUND16X4F32SASYM(x30, x31);

            AE_S16X4_X(_x0, _py, 0);
            AE_S16X4_X(_x1, _py, (1 * _v) * sizeof(complex_fract16));
            AE_S16X4_X(_x2, _py, (2 * _v) * sizeof(complex_fract16));
            AE_S16X4_X(_x3, _py, (3 * _v) * sizeof(complex_fract16));

            _x0 = AE_ABS16S(_x0);
            _x1 = AE_ABS16S(_x1);
            _x2 = AE_ABS16S(_x2);
            _x3 = AE_ABS16S(_x3);

            acc16 = AE_OR16(acc16, _x0);
            acc16 = AE_OR16(acc16, _x1);
            acc16 = AE_OR16(acc16, _x2);
            acc16 = AE_OR16(acc16, _x3);
            _py++;
        } //for (j = 0; j < _v; j++)
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
} //stage_inner_DFT4_16x16_ie, hifi3
#endif //#if (XCHAL_HAVE_HIFI3Z)
