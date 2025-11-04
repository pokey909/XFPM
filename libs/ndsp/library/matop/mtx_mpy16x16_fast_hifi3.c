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
  NatureDSP Signal Processing Library. Matrix Operations
    Matrix Multiply
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/
/* Signal Processing Library API. */
#include "NatureDSP_Signal_matop.h"
#include "NatureDSP_types.h"
#include "common.h"

#ifdef COMPILER_MSVC
#include <malloc.h>
#else
#include <alloca.h>
#endif

/*-------------------------------------------------------------------------
  Matrix Multiply
  These functions compute the expression z = 2^lsh * x * y for the matrices 
  x and y. The columnar dimension of x must match the row dimension of y. 
  The resulting matrix has the same number of rows as x and the same number 
  of columns as y.
  NOTE: lsh factor is not relevant for floating point routines.

  Functions require scratch memory for storing intermediate data. This 
  scratch memory area should be aligned on 8 byte boundary and its size is 
  calculated by macros SCRATCH_MTX_MPY24X24(M,N,P), 
  SCRATCH_MTX_MPY32X32(M,N,P), SCRATCH_MTX_MPY16X16(M,N,P).

  Two versions of functions available: regular version (mtx_mpy32x32,
  mtx_mpy24x24, mtx_mpy16x16, mtx_mpyf) with arbitrary arguments and faster
  version (mtx_mpy32x32_fast, mtx_mpy24x24_fast, mtx_mpy16x16_fast, 
  mtx_mpyf_fast) that apply some restrictions.

  Precision:
  32x32 32-bit inputs, 32-bit output
  24x24 24-bit inputs, 24-bit output
  16x16 16-bit inputs, 16-bit output
  f     floating point

  Input:
  x[M*N]      input matrix x, Q15, Q31 or floating point
  y[N*P]      input matrix y, Q15, Q31 or floating point
  M           number of rows in matrix x and z
  N           number of columns in matrix x and number of rows in matrix y
  P           number of columns in matrices y and z
  lsh         left shift applied to the result (applied to the fixed-
              point functions only) 
  Output:
  z[M*P]      output matrix z, Q15, Q31 or floating point 
  Scratch:
  pScr        size in bytes defined by macros SCRATCH_MTX_MPY32X32,
              SCRATCH_MTX_MPY24X24, SCRATCH_MTX_MPY16X16

  Restrictions:
  For regular routines (mtx_mpy32x32, mtx_mpy24x24, mtx_mpy16x16, mtx_mpyf):
  x,y,z   should not overlap

  For faster routines (mtx_mpy32x32_fast, mtx_mpy24x24_fast, 
                       mtx_mpy16x16_fast, mtx_mpyf_fast):
  x,y,z   should not overlap
  x,y,z   aligned on 8-byte boundary
  M,N,P   multiplies of 4
  lsh     should be in range:
          -31...31 for mtx_mpy32x32, mtx_mpy32x32_fast,
                   mtx_mpy24x24, mtx_mpy24x24_fast
          -15...15 for mtx_mpy16x16, mtx_mpy16x16_fast  

-------------------------------------------------------------------------*/
void mtx_mpy16x16_fast (  
                     int16_t* restrict z,
               const int16_t* restrict x,
               const int16_t* restrict y,
               int M, int N, int P, int lsh )

{
#ifndef AE_MULAAAAQ16
    int m, n, p;
          ae_int16x4 * restrict pz;
    const ae_int16x4 * restrict py;
    const ae_p16x2s * restrict px;
    const ae_p16x2s * restrict px0;
    const ae_p16x2s * restrict px1;
    ae_int16x4 t1, t2;

    NASSERT(N % 4 == 0);
    NASSERT(M % 4 == 0);
    NASSERT(P % 4 == 0);
    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT_ALIGN8(z);
    NASSERT(lsh >= -15 && lsh <= 15);
    if ((M <= 0) || (P <= 0)) return;
    if (N <= 0)
    {
        for (m = 0; m < M * P; m++) z[m] = 0;
        return;
    }

    for (p = 0; p < P; p += 4)
    {
        pz = (ae_int16x4 *)(z);
        px = (const ae_p16x2s *)(x - 2);
        for (m = 0; m < M; m += 2)
        {
            ae_int64    C0, C1, C2, C3, C4, C5, C6, C7;
            ae_int32x2  vah, vbh, vch, vdh;

            C0 = C1 = C2 = C3 = C4 = C5 = C6 = C7 = AE_ZERO64();

            px0 = (const ae_p16x2s *)px;
            px1 = (const ae_p16x2s *)XT_ADDX2(N, (uintptr_t)px0);
            py = (const ae_int16x4 *)(y);
            for (n = 0; n < N; n += 2)
            {
                ae_int32x2 x0, x1;
                ae_int16x4 y1, y2;

                AE_L16X2M_IU(x0, px0, 4);
                AE_L16X2M_IU(x1, px1, 4);

                AE_L16X4_XP(t1, py, P*sizeof(int16_t));
                AE_L16X4_XP(t2, py, P*sizeof(int16_t));
                y1 = AE_SEL16_7362(t1, t2);
                y2 = AE_SEL16_5410(t2, t1);
                y2 = AE_SEL16_5146(y2, t2);

                AE_MULAAD32X16_H3_L2(C0, x0, y1);
                AE_MULAAD32X16_H1_L0(C1, x0, y1);
                AE_MULAAD32X16_H3_L2(C2, x0, y2);
                AE_MULAAD32X16_H1_L0(C3, x0, y2);
                AE_MULAAD32X16_H3_L2(C4, x1, y1);
                AE_MULAAD32X16_H1_L0(C5, x1, y1);
                AE_MULAAD32X16_H3_L2(C6, x1, y2);
                AE_MULAAD32X16_H1_L0(C7, x1, y2);
            }
            vah = AE_TRUNCA32X2F64S(C0, C1, lsh + 25);
            vbh = AE_TRUNCA32X2F64S(C2, C3, lsh + 25);
            vch = AE_TRUNCA32X2F64S(C4, C5, lsh + 25);
            vdh = AE_TRUNCA32X2F64S(C6, C7, lsh + 25);
            t1 = AE_ROUND16X4F32SASYM(vah, vbh);
            t2 = AE_ROUND16X4F32SASYM(vch, vdh);
            AE_S16X4_XP(t1, pz, P*sizeof(int16_t));
            AE_S16X4_XP(t2, pz, P*sizeof(int16_t));

            px = (const ae_p16x2s *)XT_ADDX4(N, (uintptr_t)px);
        }
        y += 4;
        z += 4;
    }
#else
    int m, n, p;
    const ae_int16x4 * restrict px;
    const ae_int16x4 * restrict px0;
    const ae_int16x4 * restrict px1;
    const ae_int16x4 * restrict py;
          ae_int16x4 * restrict pz;
    ae_int16x4 t0, t1, t2, t3;
    ae_int16x4 x0, x1;
    ae_int16x4 y0, y1, y2, y3;
    int16_t * restrict scr;

    NASSERT(lsh >= -15 && lsh <= 15);
    NASSERT(N % 4 == 0);
    NASSERT(M % 4 == 0);
    NASSERT(P % 4 == 0);
    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT_ALIGN8(z);
    if ((M <= 0) || (P <= 0)) return;
    if (N <= 0)
    {
        for (m = 0; m < M * P; m++) z[m] = 0;
        return;
    }

    scr = (int16_t *)alloca(N * 4 * sizeof(int16_t));
    NASSERT_ALIGN8(scr);

    for (p = 0; p < P; p += 4)
    {
        py = (const ae_int16x4 *)y;
        pz = (ae_int16x4 *)scr;
        for (n = 0; n < N; n += 4)
        {
            AE_L16X4_XP(y0, py, P*sizeof(int16_t));
            AE_L16X4_XP(y1, py, P*sizeof(int16_t));
            AE_L16X4_XP(y2, py, P*sizeof(int16_t));
            AE_L16X4_XP(y3, py, P*sizeof(int16_t));

            t1 = AE_SEL16_5410(y0, y1);
            t0 = AE_SEL16_7632(y0, y1);
            t3 = AE_SEL16_5410(y2, y3);
            t2 = AE_SEL16_7632(y2, y3);
            y1 = AE_SEL16_6420(t0, t2);
            y0 = AE_SEL16I(t0, t2, 7);
            y3 = AE_SEL16_6420(t1, t3);
            y2 = AE_SEL16I(t1, t3, 7);

            AE_S16X4_IP(y0, pz, sizeof(ae_int16x4));
            AE_S16X4_IP(y1, pz, sizeof(ae_int16x4));
            AE_S16X4_IP(y2, pz, sizeof(ae_int16x4));
            AE_S16X4_IP(y3, pz, sizeof(ae_int16x4));
        }

        px = (const ae_int16x4 *)x;
        pz = (ae_int16x4 *)z;
        for (m = 0; m < M; m += 2)
        {
            ae_int32x2 a0, a1, a2, a3;
            ae_int64 B0, B1, B2, B3, B4, B5, B6, B7;
            py = (const ae_int16x4 *)scr;
            px0 = (const ae_int16x4 *)px;
            px1 = (const ae_int16x4 *)XT_ADDX2(N, (uintptr_t)px);

            B0 = B1 = B2 = B3 = B4 = B5 = B6 = B7 = AE_ZERO64();
            for (n = 0; n < N; n += 4)
            {
                AE_L16X4_IP(y0, py, sizeof(ae_int16x4));
                AE_L16X4_IP(y1, py, sizeof(ae_int16x4));
                AE_L16X4_IP(y2, py, sizeof(ae_int16x4));
                AE_L16X4_IP(y3, py, sizeof(ae_int16x4));

                AE_L16X4_XP(x0, px0, sizeof(ae_int16x4));
                AE_L16X4_XP(x1, px1, sizeof(ae_int16x4));

                AE_MULAAAAQ16(B0, x0, y0);
                AE_MULAAAAQ16(B1, x0, y1);
                AE_MULAAAAQ16(B2, x0, y2);
                AE_MULAAAAQ16(B3, x0, y3);
                AE_MULAAAAQ16(B4, x1, y0);
                AE_MULAAAAQ16(B5, x1, y1);
                AE_MULAAAAQ16(B6, x1, y2);
                AE_MULAAAAQ16(B7, x1, y3);
            }
            a0 = AE_TRUNCA32X2F64S(B0, B1, lsh + 33);
            a1 = AE_TRUNCA32X2F64S(B2, B3, lsh + 33);
            a2 = AE_TRUNCA32X2F64S(B4, B5, lsh + 33);
            a3 = AE_TRUNCA32X2F64S(B6, B7, lsh + 33);
            t0 = AE_ROUND16X4F32SASYM(a0, a1);
            t1 = AE_ROUND16X4F32SASYM(a2, a3);
            AE_S16X4_XP(t0, pz, P*sizeof(int16_t));
            AE_S16X4_XP(t1, pz, P*sizeof(int16_t));

            px = (const ae_int16x4 *)XT_ADDX4(N, (uintptr_t)px);
        }
        z += 4;
        y += 4;
    }
#endif
} /* mtx_mpy16x16_fast() */

