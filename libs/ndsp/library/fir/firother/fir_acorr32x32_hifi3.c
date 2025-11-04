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
    Real data circular auto-correlation, 32x32-bit
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Circular Autocorrelation 
  Estimates the auto-correlation of vector x. Returns autocorrelation of 
  length N.

  Precision: 
  16x16     16-bit data, 16-bit outputs
  24x24     24-bit data, 24-bit outputs
  32x32     32-bit data, 32-bit outputs
  f         floating point

  Input:
  x[N]      input data, Q15, Q31 or floating point
  N         length of x
  Output:
  r[N]      output data, Q15, Q31 or floating point

  Restrictions:
  x,r       should not overlap
  x,r       aligned on an 8-bytes boundary
  N         multiple of 4 and >0
-------------------------------------------------------------------------*/
/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

void fir_acorr32x32( int32_t * restrict r,
               const int32_t * restrict x,
               int N )
{
    //
    // Circular autocorrelation algorithm:
    //
    //   r[n] = sum( x[mod(n+m,N)]*x[m] )
    //        m=0..N-1
    //
    //   where n = 0..N-1
    //

    const ae_f32x2 * X;
    const ae_f32x2 * S;
    const ae_f32x2 * Y;
    ae_f32x2 * restrict R;

    ae_f64     q0, q1, q2, q3;
    ae_int32x2 t0, t1, t2;
    ae_f32x2 x0, x1, x2;
    ae_f32x2 y0;
    ae_f32x2 y0_low, y0_high;

    int n, m;

    NASSERT(r);
    NASSERT(x);
    NASSERT_ALIGN(r, 8);
    NASSERT_ALIGN(x, 8);
    NASSERT((N > 0) && ((N % 4) == 0));

    //
    // Setup pointers and circular addressing for array x[N].
    //

    X = (const ae_f32x2 *)x;
    R = (ae_f32x2 *)r;

    WUR_AE_CBEGIN0((uintptr_t)(x + 0));
    WUR_AE_CEND0  ((uintptr_t)(x + N));

    for (n = 0; n < (N >> 2); n++)
    {
        AE_L32X2_XC(t0, X, +8); x0 = (t0);
        AE_L32X2_XC(t1, X, +8); x1 = (t1);
        S = X;
        AE_L32X2_XC(t2, S, +8); x2 = (t2);

        Y = (const ae_f32x2 *)x;

        AE_L32X2_IP(t0, Y, +8); y0 = (t0);

        // Q16.47 <- Q16.47 + [ Q31*Q31 - 15 ] w/ symmetric rounding
#if (defined(AE_MULAAFD32R_HH_LL))
        q0 = 0;
        AE_MULAAFD32R_HH_LL(q0, x0, y0);
#else
        q0 = AE_MULF32R_HH(x0, y0);
        AE_MULAF32R_LL(q0, x0, y0);
#endif
        
        // Q16.47 <- Q16.47 + [ Q31*Q31 - 15 ] w/ symmetric rounding
        q1 = AE_MULF32R_LH(x0, y0);
        AE_MULAF32R_LH(q1, y0, x1);

        // Q16.47 <- Q16.47 + [ Q31*Q31 - 15 ] w/ symmetric rounding
#if (defined(AE_MULAAFD32R_HH_LL))
        q2 = 0;
        AE_MULAAFD32R_HH_LL(q2, x1, y0);
#else
        q2 = AE_MULF32R_HH(x1, y0);
        AE_MULAF32R_LL(q2, x1, y0);
#endif
        // Q16.47 <- Q16.47 + [ Q31*Q31 - 15 ] w/ symmetric rounding
        q3 = AE_MULF32R_LH(x1, y0);
        AE_MULAF32R_LH(q3, y0, x2);

        x0 = x1;
        x1 = x2;
        AE_L32X2_XC(t2, S, +8);
        x2 = (t2);

#if 0
        __Pragma("loop_count min=2")
        for (m = 0; m < (N >> 1) - 1; m++)
        {
            AE_L32X2_IP(t0, Y, +8);
            y0 = (t0);

            // Q16.47 <- Q16.47 + [ Q31*Q31 - 15 ] w/ symmetric rounding
#if (defined(AE_MULAAFD32R_HH_LL))
            AE_MULAAFD32R_HH_LL(q0, x0, y0);
#else
            AE_MULAF32R_HH(q0, x0, y0);
            AE_MULAF32R_LL(q0, x0, y0);
#endif

            // Q16.47 <- Q16.47 + [ Q31*Q31 - 15 ] w/ symmetric rounding
            AE_MULAF32R_LH(q1, x0, y0);
            AE_MULAF32R_LH(q1, y0, x1);

            // Q16.47 <- Q16.47 + [ Q31*Q31 - 15 ] w/ symmetric rounding
#if (defined(AE_MULAAFD32R_HH_LL))
            AE_MULAAFD32R_HH_LL(q2, x1, y0);
#else
            AE_MULAF32R_HH(q2, x1, y0);
            AE_MULAF32R_LL(q2, x1, y0);
#endif

            // Q16.47 <- Q16.47 + [ Q31*Q31 - 15 ] w/ symmetric rounding
            AE_MULAF32R_LH(q3, x1, y0);
            AE_MULAF32R_LH(q3, y0, x2);

            x0 = x1;
            x1 = x2;
            AE_L32X2_XC(t2, S, +8);
            x2 = (t2);
        }
#else
        __Pragma("loop_count min=1")
        for (m = 0; m < (N >> 1) - 1; m++)
        {
            //AE_L32X2_IP(t0, Y, +8);
            //y0 = (t0);
            AE_L32_IP(t0, castxcc(ae_int32,Y), +4);
            y0_high = (t0);
            AE_L32_IP(t0, castxcc(ae_int32,Y), +4);
            y0_low = (t0);

            // Q16.47 <- Q16.47 + [ Q31*Q31 - 15 ] w/ symmetric rounding
            AE_MULAF32R_HH(q0, x0, y0_high);
            AE_MULAF32R_LL(q0, x0, y0_low);

            // Q16.47 <- Q16.47 + [ Q31*Q31 - 15 ] w/ symmetric rounding
            AE_MULAF32R_LL(q1, x0, y0_high);
            AE_MULAF32R_LH(q1, y0_low, x1);

            // Q16.47 <- Q16.47 + [ Q31*Q31 - 15 ] w/ symmetric rounding
            AE_MULAF32R_HH(q2, x1, y0_high);
            AE_MULAF32R_LL(q2, x1, y0_low);

            // Q16.47 <- Q16.47 + [ Q31*Q31 - 15 ] w/ symmetric rounding
            AE_MULAF32R_LL(q3, x1, y0_high);
            AE_MULAF32R_LH(q3, y0_low, x2);

            x0 = x1;
            x1 = x2;
            AE_L32X2_XC(t2, S, +8);
            x2 = (t2);
        }
#endif

        // Convert and save 4 outputs.
        // 2xQ31 <- 2xQ16.47 - 16 w/ rounding and saturation.
        AE_S32RA64S_IP(q0, castxcc(ae_f32, R), +4);
        AE_S32RA64S_IP(q1, castxcc(ae_f32, R), +4);
        AE_S32RA64S_IP(q2, castxcc(ae_f32, R), +4);
        AE_S32RA64S_IP(q3, castxcc(ae_f32, R), +4);
    }
} /* fir_acorr32x32() */
