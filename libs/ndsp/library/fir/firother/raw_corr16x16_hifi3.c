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
  NatureDSP Signal Processing Library. FIR part
    helper for correlation/convolution
    C code optimized for HiFi3
*/
/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "raw_corr16x16.h"

/*-----------------------------------------------------
    raw correlation:
    Input:
    x[N+M-1] padded with extra 3 zeroes
    y[M]
    Output:
    r[N]
    restriction:
    M should be a multiple of 4 and >0
-----------------------------------------------------*/
void raw_corr16x16(int16_t * r, const int16_t * restrict x, const int16_t * restrict y, int N, int M)
#if !(defined(AE_MULZAAAAQ16))
{
        /* code with Quad MAC support */
    const ae_int16x4 * pX0;
    const ae_int16x4 * pX1;
    const ae_int16x4 * pX2;
    const ae_int16x4 * pX3;
    const ae_int16x4 * S0;
    const ae_int16x4 * S0_1;
    const ae_int16x4 * pY;

    ae_int16x4 * restrict pR;

    ae_int64   q0, q1, q2, q3;
    ae_int64   q4, q5, q6, q7;
    ae_int16x4 x0, x1, x2, x3;
    ae_int16x4 y0;
    ae_f32x2   t0, t1;

    ae_valign ar1, ar2, ar3, ar;

    int n, m;

    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(y, 8);
    NASSERT(M > 0 && M % 2 == 0);

    /* compute N&~3 correlation results */
    ar = AE_ZALIGN64();

    pX0 = (const ae_int16x4 *)x;
    pR = (ae_int16x4 *)r;

    for (n = 0; n < (N&~7); n += 8)
    {
        pX1 = (const ae_int16x4 *)((int16_t *) pX0 + 1);
        pX2 = (const ae_int16x4 *)((int16_t *) pX0 + 2);
        pX3 = (const ae_int16x4 *)((int16_t *) pX0 + 3);

        ar1 = AE_LA64_PP(pX1);
        ar2 = AE_LA64_PP(pX2);
        ar3 = AE_LA64_PP(pX3);

        AE_L16X4_IP(x0, pX0, +8);
        AE_LA16X4_IP(x1, ar1, pX1);
        AE_LA16X4_IP(x2, ar2, pX2);
        AE_LA16X4_IP(x3, ar3, pX3);

        S0 = pX0;

        pY = (const ae_int16x4 *)y;
        AE_L16X4_IP(y0, pY, 8);
        {
            ae_f32x2 yy0, yy1;
            yy0 = AE_CVT32X2F16_32(y0);
            yy1 = AE_CVT32X2F16_10(y0);
            q0 = AE_MULZAAD32X16_H3_L2(yy0, x0);
            AE_MULAAD32X16_H1_L0  (q0, yy1, x0);
            q1 = AE_MULZAAD32X16_H3_L2(yy0, x1);
            AE_MULAAD32X16_H1_L0  (q1, yy1, x1);
            q2 = AE_MULZAAD32X16_H3_L2(yy0, x2);
            AE_MULAAD32X16_H1_L0  (q2, yy1, x2);
            q3 = AE_MULZAAD32X16_H3_L2(yy0, x3);
            AE_MULAAD32X16_H1_L0  (q3, yy1, x3);
            
            AE_L16X4_IP(x0, pX0, +8);
            S0_1 = pX0;
            AE_LA16X4_IP(x1, ar1, pX1);
            AE_LA16X4_IP(x2, ar2, pX2);
            AE_LA16X4_IP(x3, ar3, pX3);
            
            q4 = AE_MULZAAD32X16_H3_L2(yy0, x0);
            AE_MULAAD32X16_H1_L0  (q4, yy1, x0);
            q5 = AE_MULZAAD32X16_H3_L2(yy0, x1);
            AE_MULAAD32X16_H1_L0  (q5, yy1, x1);
            q6 = AE_MULZAAD32X16_H3_L2(yy0, x2);
            AE_MULAAD32X16_H1_L0  (q6, yy1, x2);
            q7 = AE_MULZAAD32X16_H3_L2(yy0, x3);
            AE_MULAAD32X16_H1_L0  (q7, yy1, x3);
        }

        for (m = 0; m < (M >> 2) - 1; m++)
        {
            {
                ae_f32x2 yy0, yy1;
                AE_L16X4_IP(x0, S0, +8);
                AE_L16X4_IP(y0, pY, 8);
                yy0 = AE_CVT32X2F16_32(y0);
                yy1 = AE_CVT32X2F16_10(y0);
                
                AE_MULAAD32X16_H3_L2(q0, yy0, x0);
                AE_MULAAD32X16_H1_L0(q0, yy1, x0);
                AE_MULAAD32X16_H3_L2(q1, yy0, x1);
                AE_MULAAD32X16_H1_L0(q1, yy1, x1);
                AE_MULAAD32X16_H3_L2(q2, yy0, x2);
                AE_MULAAD32X16_H1_L0(q2, yy1, x2);
                AE_MULAAD32X16_H3_L2(q3, yy0, x3);
                AE_MULAAD32X16_H1_L0(q3, yy1, x3);
                
                AE_L16X4_IP(x0, S0_1, 8);
                AE_LA16X4_IP(x1, ar1, pX1);
                AE_LA16X4_IP(x2, ar2, pX2);
                AE_LA16X4_IP(x3, ar3, pX3);
                
                AE_MULAAD32X16_H3_L2(q4, yy0, x0);
                AE_MULAAD32X16_H1_L0(q4, yy1, x0);
                AE_MULAAD32X16_H3_L2(q5, yy0, x1);
                AE_MULAAD32X16_H1_L0(q5, yy1, x1);
                AE_MULAAD32X16_H3_L2(q6, yy0, x2);
                AE_MULAAD32X16_H1_L0(q6, yy1, x2);
                AE_MULAAD32X16_H3_L2(q7, yy0, x3);
                AE_MULAAD32X16_H1_L0(q7, yy1, x3);
            }
        }
        t0 = AE_TRUNCA32X2F64S(q0, q1, 17);
        t1 = AE_TRUNCA32X2F64S(q2, q3, 17);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);

        t0 = AE_TRUNCA32X2F64S(q4, q5, 17);
        t1 = AE_TRUNCA32X2F64S(q6, q7, 17);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
    }
    /* process last 1...7 samples */
    N &= 7;
    if (N & 4)
    {
        pX1 = (const ae_int16x4 *)((int16_t *)pX0 + 1);
        pX2 = (const ae_int16x4 *)((int16_t *)pX0 + 2);
        pX3 = (const ae_int16x4 *)((int16_t *)pX0 + 3);

        ar1 = AE_LA64_PP(pX1);
        ar2 = AE_LA64_PP(pX2);
        ar3 = AE_LA64_PP(pX3);

        AE_L16X4_IP(x0, pX0, +8);
        AE_LA16X4_IP(x1, ar1, pX1);
        AE_LA16X4_IP(x2, ar2, pX2);
        AE_LA16X4_IP(x3, ar3, pX3);

        S0 = pX0;

        pY = (const ae_int16x4 *)y;
        AE_L16X4_IP(y0, pY, 8);
        {
            ae_f32x2 yy0, yy1;
            yy0 = AE_CVT32X2F16_32(y0);
            yy1 = AE_CVT32X2F16_10(y0);
            q0 = AE_MULZAAD32X16_H3_L2(yy0, x0);
            AE_MULAAD32X16_H1_L0(q0, yy1, x0);
            q1 = AE_MULZAAD32X16_H3_L2(yy0, x1);
            AE_MULAAD32X16_H1_L0(q1, yy1, x1);
            q2 = AE_MULZAAD32X16_H3_L2(yy0, x2);
            AE_MULAAD32X16_H1_L0(q2, yy1, x2);
            q3 = AE_MULZAAD32X16_H3_L2(yy0, x3);
            AE_MULAAD32X16_H1_L0(q3, yy1, x3);
        }

        for (m = 0; m < (M >> 2) - 1; m++)
        {
            AE_L16X4_IP(x0, S0, +8);
            AE_LA16X4_IP(x1, ar1, pX1);
            AE_LA16X4_IP(x2, ar2, pX2);
            AE_LA16X4_IP(x3, ar3, pX3);

            AE_L16X4_IP(y0, pY, 8);
            {
                ae_f32x2 yy0, yy1;
                yy0 = AE_CVT32X2F16_32(y0);
                yy1 = AE_CVT32X2F16_10(y0);
                AE_MULAAD32X16_H3_L2(q0, yy0, x0);
                AE_MULAAD32X16_H1_L0(q0, yy1, x0);
                AE_MULAAD32X16_H3_L2(q1, yy0, x1);
                AE_MULAAD32X16_H1_L0(q1, yy1, x1);
                AE_MULAAD32X16_H3_L2(q2, yy0, x2);
                AE_MULAAD32X16_H1_L0(q2, yy1, x2);
                AE_MULAAD32X16_H3_L2(q3, yy0, x3);
                AE_MULAAD32X16_H1_L0(q3, yy1, x3);
            }
        }
        t0 = AE_TRUNCA32X2F64S(q0, q1, 17);
        t1 = AE_TRUNCA32X2F64S(q2, q3, 17);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
    }
    AE_SA64POS_FP(ar, pR);

    N &= 3;
    if (N)
    {
        pX1 = (const ae_int16x4 *)((int16_t *)pX0 + 1);
        pX2 = (const ae_int16x4 *)((int16_t *)pX0 + 2);

        ar1 = AE_LA64_PP(pX1);
        ar2 = AE_LA64_PP(pX2);

        AE_L16X4_IP(x0, pX0, +8);
        AE_LA16X4_IP(x1, ar1, pX1);
        AE_LA16X4_IP(x2, ar2, pX2);

        pY = (const ae_int16x4 *)y;
        AE_L16X4_IP(y0, pY, 8);
        {
            ae_f32x2 yy0, yy1;
            yy0 = AE_CVT32X2F16_32(y0);
            yy1 = AE_CVT32X2F16_10(y0);
            q0 = AE_MULZAAD32X16_H3_L2(yy0, x0);
            AE_MULAAD32X16_H1_L0(q0, yy1, x0);
            q1 = AE_MULZAAD32X16_H3_L2(yy0, x1);
            AE_MULAAD32X16_H1_L0(q1, yy1, x1);
            q2 = AE_MULZAAD32X16_H3_L2(yy0, x2);
            AE_MULAAD32X16_H1_L0(q2, yy1, x2);
        }

        for (m = 0; m < (M >> 2) - 1; m++)
        {
            AE_L16X4_IP(x0, pX0, +8);
            AE_LA16X4_IP(x1, ar1, pX1);
            AE_LA16X4_IP(x2, ar2, pX2);

            AE_L16X4_IP(y0, pY, 8);
            {
                ae_f32x2 yy0, yy1;
                yy0 = AE_CVT32X2F16_32(y0);
                yy1 = AE_CVT32X2F16_10(y0);
                AE_MULAAD32X16_H3_L2(q0, yy0, x0);
                AE_MULAAD32X16_H1_L0(q0, yy1, x0);
                AE_MULAAD32X16_H3_L2(q1, yy0, x1);
                AE_MULAAD32X16_H1_L0(q1, yy1, x1);
                AE_MULAAD32X16_H3_L2(q2, yy0, x2);
                AE_MULAAD32X16_H1_L0(q2, yy1, x2);
            }
        }
        t0 = AE_TRUNCA32X2F64S(q1, q0, 17);
        t1 = AE_TRUNCA32X2F64S(q2, q2, 17);
        x0 = AE_ROUND16X4F32SASYM(t1, t0);

        AE_S16_0_IP(x0, castxcc(ae_int16, pR), +2);
        if (N > 1) { x0 = AE_SEL16_4321(x0, x0); AE_S16_0_IP(x0, castxcc(ae_int16, pR), +2); }
        if (N > 2) { x0 = AE_SEL16_4321(x0, x0); AE_S16_0_IP(x0, castxcc(ae_int16, pR), +2); }
    }
} /* raw_corr16x16() */
#else
{
    /* code with Quad MAC support */
    const ae_int16x4 * pX0;
    const ae_int16x4 * pX1;
    const ae_int16x4 * pX2;
    const ae_int16x4 * pX3;
    const ae_int16x4 * S0;
    const ae_int16x4 * S0_1;
    const ae_int16x4 * pY;

    ae_int16x4 * restrict pR;

    ae_int64   q0, q1, q2, q3;
    ae_int64   q4, q5, q6, q7;
    ae_int16x4 x0, x1, x2, x3;
    ae_int16x4 y0;
    ae_f32x2   t0, t1;

    ae_valign ar1, ar2, ar3, ar;

    int n, m;

    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(y, 8);
    NASSERT(M > 0 && M % 2 == 0);

    /* compute N&~3 correlation results */
    ar = AE_ZALIGN64();

    pX0 = (const ae_int16x4 *)x;
    pR = (ae_int16x4 *)r;

    for (n = 0; n < (N&~7); n += 8)
    {

        AE_L16X4_IP(x0, pX0, +8);
        S0 = pX0;
        AE_L16X4_IP(x1, pX0, +8);
        S0_1 = pX0;
        AE_L16X4_IP(x2, S0_1, +8);

        pY = (const ae_int16x4 *)y;
        AE_L16X4_IP(y0, pY, 8);

        AE_MULFQ16X2_FIR_3(q0, q1, x0, x1, y0);
        AE_MULFQ16X2_FIR_1(q2, q3, x0, x1, y0);
        AE_MULFQ16X2_FIR_3(q4, q5, x1, x2, y0);
        AE_MULFQ16X2_FIR_1(q6, q7, x1, x2, y0);

        x0 = x1; x1=x2;

        for (m = 0; m < (M >> 2) - 1; m++)
        {
            AE_L16X4_IP(x2, S0_1, +8);
            AE_L16X4_IP(y0, pY, 8);

            AE_MULAFQ16X2_FIR_3(q0, q1, x0, x1, y0);
            AE_MULAFQ16X2_FIR_1(q2, q3, x0, x1, y0);
            AE_MULAFQ16X2_FIR_3(q4, q5, x1, x2, y0);
            AE_MULAFQ16X2_FIR_1(q6, q7, x1, x2, y0);
            x0 = x1; x1 = x2;
        }
        t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
        t1 = AE_TRUNCA32X2F64S(q2, q3, 32);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);

        t0 = AE_TRUNCA32X2F64S(q4, q5, 32);
        t1 = AE_TRUNCA32X2F64S(q6, q7, 32);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
    }
    /* process last 1...7 samples */
    N &= 7;
    if (N & 4)
    {
        pX1 = (const ae_int16x4 *)((int16_t *)pX0 + 1);
        pX2 = (const ae_int16x4 *)((int16_t *)pX0 + 2);
        pX3 = (const ae_int16x4 *)((int16_t *)pX0 + 3);

        ar1 = AE_LA64_PP(pX1);
        ar2 = AE_LA64_PP(pX2);
        ar3 = AE_LA64_PP(pX3);

        AE_L16X4_IP (x0, pX0, +8);
        AE_LA16X4_IP(x1, ar1, pX1);
        AE_LA16X4_IP(x2, ar2, pX2);
        AE_LA16X4_IP(x3, ar3, pX3);

        S0 = pX0;

        pY = (const ae_int16x4 *)y;
        AE_L16X4_IP(y0, pY, 8);

        q0 = AE_MULZAAAAQ16(x0, y0);
        q1 = AE_MULZAAAAQ16(x1, y0);
        q2 = AE_MULZAAAAQ16(x2, y0);
        q3 = AE_MULZAAAAQ16(x3, y0);

        for (m = 0; m < (M >> 2) - 1; m++)
        {
            AE_L16X4_IP(x0, S0, +8);
            AE_LA16X4_IP(x1, ar1, pX1);
            AE_LA16X4_IP(x2, ar2, pX2);
            AE_LA16X4_IP(x3, ar3, pX3);

            AE_L16X4_IP(y0, pY, 8);

            AE_MULAAAAQ16(q0, x0, y0);
            AE_MULAAAAQ16(q1, x1, y0);
            AE_MULAAAAQ16(q2, x2, y0);
            AE_MULAAAAQ16(q3, x3, y0);
        }
        t0 = AE_TRUNCA32X2F64S(q0, q1, 33);
        t1 = AE_TRUNCA32X2F64S(q2, q3, 33);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
    }
    AE_SA64POS_FP(ar, pR);

    N &= 3;
    if (N)
    {
        pX1 = (const ae_int16x4 *)((int16_t *)pX0 + 1);
        pX2 = (const ae_int16x4 *)((int16_t *)pX0 + 2);

        ar1 = AE_LA64_PP(pX1);
        ar2 = AE_LA64_PP(pX2);

        AE_L16X4_IP(x0, pX0, +8);
        AE_LA16X4_IP(x1, ar1, pX1);
        AE_LA16X4_IP(x2, ar2, pX2);

        pY = (const ae_int16x4 *)y;
        AE_L16X4_IP(y0, pY, 8);

        q0 = AE_MULZAAAAQ16(x0, y0);
        q1 = AE_MULZAAAAQ16(x1, y0);
        q2 = AE_MULZAAAAQ16(x2, y0);

        for (m = 0; m < (M >> 2) - 1; m++)
        {
            AE_L16X4_IP(x0, pX0, +8);
            AE_LA16X4_IP(x1, ar1, pX1);
            AE_LA16X4_IP(x2, ar2, pX2);

            AE_L16X4_IP(y0, pY, 8);

            AE_MULAAAAQ16(q0, x0, y0);
            AE_MULAAAAQ16(q1, x1, y0);
            AE_MULAAAAQ16(q2, x2, y0);
        }

        t0 = AE_TRUNCA32X2F64S(q1, q0, 33);
        t1 = AE_TRUNCA32X2F64S(q2, q2, 33);
        x0 = AE_ROUND16X4F32SASYM(t1, t0);

        AE_S16_0_IP(x0, castxcc(ae_int16, pR), +2);
        if (N > 1) 
        { 
          x0 = AE_SEL16_4321(x0, x0); 
          AE_S16_0_IP(x0, castxcc(ae_int16, pR), +2); 
        }
        if (N > 2) 
        {
          x0 = AE_SEL16_4321(x0, x0);
          AE_S16_0_IP(x0, castxcc(ae_int16, pR), +2); 
        }
    }
} /* raw_corr16x16() */
#endif
