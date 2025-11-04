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
    Real data circular convolution, 16x16-bit
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Circular Convolution
  Performs circular convolution between vectors x (of length N) and y (of 
  length M)  resulting in vector r of length N.

  Precision: 
  16x16     16x16-bit data, 16-bit outputs
  24x24     24x24-bit data, 24-bit outputs
  32x16     32x16-bit data, 32-bit outputs 
  32x32     32x32-bit data, 32-bit outputs
  f         floating point

  Input:
  x[N]      input data, Q15, Q31 or floating point
  y[M]      input data, Q15, Q31 or floating point
  N         length of x
  M         length of y
  Output:
  r[N]      output data, Q15, Q31 or floating point

  Restriction:
  x,y,r     should not overlap
  x,y,r     aligned on an 8-bytes boundary
  N,M       multiples of 4 and >0
-------------------------------------------------------------------------*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

void fir_convol16x16( int16_t * restrict r,
                const int16_t * restrict x,
                const int16_t * restrict y,
                int N, int M )
#if !(defined(AE_MULZAAAAQ16) )
{
    /* code with Quad MAC option */
    //
    // Circular convolution algorithm:
    //
    //   r[n] = sum( x[mod(n-m,N)]*y[m] )
    //        m=0..M-1
    //
    //   where n = 0..N-1
    //
    const int16_t    * xn;
    const ae_int16x4 * pX0;
    const ae_int16x4 * pX1;
    const ae_int16x4 * pX2;
    const ae_int16x4 * pX3;
    const ae_int16x4 * pX3_1;
    const ae_int16x4 * pY;

    ae_int16x4 * restrict pR;

    ae_int64   q0, q1, q2, q3;
    ae_int64   q4, q5, q6, q7;
    ae_int16x4 x0, x1, x2, x3;
    ae_int16x4 y0;
    ae_f32x2   t0, t1;

    ae_valign ar0, ar1, ar2;

    int n, m;

    NASSERT(r);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(r, 8);
    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(y, 8);
    NASSERT((M>0) && ((M % 4) == 0));
    NASSERT((N>0) && ((N % 4) == 0));

    //
    // Setup pointers and circular addressing for array x[N].
    //
    xn = x;
    pR = (ae_int16x4 *)r;

    WUR_AE_CBEGIN0((uintptr_t)(x + 0));
    WUR_AE_CEND0  ((uintptr_t)(x + N));

    for (n = 0; n < (N&~7); n += 8)
    {
        pX0 = (const ae_int16x4 *)(xn + 4);
        pX1 = (const ae_int16x4 *)(xn + 5);
        pX2 = (const ae_int16x4 *)(xn + 6);
        pX3 = (const ae_int16x4 *)(xn + 4);

        xn += 8;

        AE_LA16X4NEG_PC(ar0, pX0);
        AE_LA16X4NEG_PC(ar1, pX1);
        AE_LA16X4NEG_PC(ar2, pX2);
        
        AE_LA16X4_RIC(x0, ar0, pX0);
        AE_LA16X4_RIC(x1, ar1, pX1);
        AE_LA16X4_RIC(x2, ar2, pX2);
        AE_L16X4_RIC (x3, pX3);

        pY = (const ae_int16x4 *)y;
        AE_L16X4_IP(y0, pY, 8);
        {
            ae_f32x2 yy0, yy1;
            yy0 = AE_CVT32X2F16_32(y0);
            yy1 = AE_CVT32X2F16_10(y0);
            q4 = AE_MULZAAD32X16_H3_L2(yy0, x0);
            AE_MULAAD32X16_H1_L0  (q4, yy1, x0);
            q5 = AE_MULZAAD32X16_H3_L2(yy0, x1);
            AE_MULAAD32X16_H1_L0  (q5, yy1, x1);
            q6 = AE_MULZAAD32X16_H3_L2(yy0, x2);
            AE_MULAAD32X16_H1_L0  (q6, yy1, x2);
            q7 = AE_MULZAAD32X16_H3_L2(yy0, x3);
            AE_MULAAD32X16_H1_L0  (q7, yy1, x3);
        }

        pX3_1 = pX3;

        AE_LA16X4_RIC(x0, ar0, pX0);
        AE_LA16X4_RIC(x1, ar1, pX1);
        AE_LA16X4_RIC(x2, ar2, pX2);
        AE_L16X4_RIC (x3, pX3_1);
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
        }

        for (m = 0; m < (M >> 2) - 1; m++)
        {
            ae_f32x2 yy0, yy1;
            AE_L16X4_RIC(x3, pX3);
            AE_L16X4_IP(y0, pY, 8);
            yy0 = AE_CVT32X2F16_32(y0);
            yy1 = AE_CVT32X2F16_10(y0);
            AE_MULAAD32X16_H3_L2(q4, yy0, x0);
            AE_MULAAD32X16_H1_L0(q4, yy1, x0);
            AE_MULAAD32X16_H3_L2(q5, yy0, x1);
            AE_MULAAD32X16_H1_L0(q5, yy1, x1);
            AE_MULAAD32X16_H3_L2(q6, yy0, x2);
            AE_MULAAD32X16_H1_L0(q6, yy1, x2);
            AE_MULAAD32X16_H3_L2(q7, yy0, x3);
            AE_MULAAD32X16_H1_L0(q7, yy1, x3);

            AE_L16X4_RIC(x3, pX3_1);
            AE_LA16X4_RIC(x0, ar0, pX0);
            AE_LA16X4_RIC(x1, ar1, pX1);
            AE_LA16X4_RIC(x2, ar2, pX2);

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
        t0 = AE_TRUNCA32X2F64S(q0, q1, 17);
        t1 = AE_TRUNCA32X2F64S(q2, q3, 17);
        AE_S16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), pR, 8);

        t0 = AE_TRUNCA32X2F64S(q4, q5, 17);
        t1 = AE_TRUNCA32X2F64S(q6, q7, 17);
        AE_S16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), pR, 8);
    }
    
    if (N & 4)
    {
        pX0 = (const ae_int16x4 *)(xn + 0);
        pX1 = (const ae_int16x4 *)(xn + 1);
        pX2 = (const ae_int16x4 *)(xn + 2);
        pX3 = (const ae_int16x4 *)(xn);

        AE_LA16X4NEG_PC(ar0, pX0);
        AE_LA16X4NEG_PC(ar1, pX1);
        AE_LA16X4NEG_PC(ar2, pX2);

        AE_L16X4_RIC (x3, pX3);
        AE_LA16X4_RIC(x0, ar0, pX0);
        AE_LA16X4_RIC(x1, ar1, pX1);
        AE_LA16X4_RIC(x2, ar2, pX2);

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
        }

        for (m = 0; m < (M >> 2) - 1; m++)
        {
            AE_L16X4_RIC(x3, pX3);
            AE_LA16X4_RIC(x0, ar0, pX0);
            AE_LA16X4_RIC(x1, ar1, pX1);
            AE_LA16X4_RIC(x2, ar2, pX2);
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
        AE_S16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), pR, 8);
    }
}
#else
{
    /* code with Quad MAC option */
    //
    // Circular convolution algorithm:
    //
    //   r[n] = sum( x[mod(n-m,N)]*y[m] )
    //        m=0..M-1
    //
    //   where n = 0..N-1
    //
    const int16_t    * xn;
    const ae_int16x4 * pX;
    const ae_int16x4 * pY;

    ae_int16x4 * restrict pR;

    ae_int64   q0, q1, q2, q3;
    ae_int64   q4, q5, q6, q7;
    ae_int16x4 x0, x1, x2;
    ae_int16x4 y0;
    ae_f32x2   t0, t1;

    int n, m;

    NASSERT(r);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(r, 8);
    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(y, 8);
    NASSERT((M>0) && ((M % 4) == 0));
    NASSERT((N>0) && ((N % 4) == 0));

    //
    // Setup pointers and circular addressing for array x[N].
    //
    xn = x;
    pR = (ae_int16x4 *)r;

    WUR_AE_CBEGIN0((uintptr_t)(x + 0));
    WUR_AE_CEND0  ((uintptr_t)(x + N));

    for (n = 0; n < (N&~7); n += 8)
    {
        pX = (const ae_int16x4 *)(xn + 4);

        xn += 8;

        AE_L16X4_RIC (x2, pX);
        AE_L16X4_RIC (x1, pX);
        AE_L16X4_RIC (x0, pX);

        pY = (const ae_int16x4 *)y;
        AE_L16X4_IP(y0, pY, 8);

        AE_MULFQ16X2_FIR_3(q7, q6, x2, x1, y0);
        AE_MULFQ16X2_FIR_1(q5, q4, x2, x1, y0);
        AE_MULFQ16X2_FIR_3(q3, q2, x1, x0, y0);
        AE_MULFQ16X2_FIR_1(q1, q0, x1, x0, y0);

        x2 = x1; x1 = x0;

        for (m = 0; m < (M >> 2) - 1; m++)
        {
            AE_L16X4_RIC(x0, pX);
            AE_L16X4_IP(y0, pY, 8);

            AE_MULAFQ16X2_FIR_3(q7, q6, x2, x1, y0);
            AE_MULAFQ16X2_FIR_1(q5, q4, x2, x1, y0);
            AE_MULAFQ16X2_FIR_3(q3, q2, x1, x0, y0);
            AE_MULAFQ16X2_FIR_1(q1, q0, x1, x0, y0);
            x2 = x1; x1 = x0;
        }

        t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
        t1 = AE_TRUNCA32X2F64S(q2, q3, 32);
        AE_S16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), pR, 8);

        t0 = AE_TRUNCA32X2F64S(q4, q5, 32);
        t1 = AE_TRUNCA32X2F64S(q6, q7, 32);
        AE_S16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), pR, 8);
    }
    
    if (N & 4)
    {
        pX = (const ae_int16x4 *)(xn);

        AE_L16X4_RIC (x2, pX);
        AE_L16X4_RIC (x1, pX);

        pY = (const ae_int16x4 *)y;
        AE_L16X4_IP(y0, pY, 8);

        AE_MULFQ16X2_FIR_3(q3, q2, x2, x1, y0);
        AE_MULFQ16X2_FIR_1(q1, q0, x2, x1, y0);
        x2 = x1;

        for (m = 0; m < (M >> 2) - 1; m++)
        {
            AE_L16X4_RIC(x1, pX);
            AE_L16X4_IP (y0, pY, 8);

            AE_MULAFQ16X2_FIR_3(q3, q2, x2, x1, y0);
            AE_MULAFQ16X2_FIR_1(q1, q0, x2, x1, y0);
            x2 = x1;
        }
        t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
        t1 = AE_TRUNCA32X2F64S(q2, q3, 32);
        AE_S16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), pR, 8);
    }
} /* fir_convol16x16() */
#endif
