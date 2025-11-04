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
void mtx_mpy24x24_fast (  
                     f24* restrict z,
               const f24* restrict x,
               const f24* restrict y,
               int M, int N, int P, int lsh )

{
    const ae_f24x2 * restrict px0;
    const ae_f24x2 * restrict px1;
    const ae_f24x2 * restrict px2;
    const ae_f24x2 * restrict px3;
    const ae_f24x2 * restrict py;
          ae_f24x2 * restrict pz;

    ae_f24x2 vx0, vx1, vx2, vx3,
             vy0, vy1, vt0, vt1,
             vz0, vz1, vz2, vz3;
    ae_f64   ACC00, ACC01, ACC10, ACC11,
             ACC20, ACC21, ACC30, ACC31;
    int n, m, p;

    NASSERT(N % 4 == 0);
    NASSERT(M % 4 == 0);
    NASSERT(P % 4 == 0);
    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT_ALIGN8(z);
    NASSERT(lsh >= -31 && lsh <= 31);
    if ((M <= 0) || (P <= 0)) return;
    if (N <= 0)
    {
        for (m = 0; m < M * P; m++) z[m] = 0;
        return;
    }


    for (p = 0; p < P; p += 2)
    {
        pz = (ae_f24x2 *)z;
        px3 = (const ae_f24x2 *)(x);

        for (m = 0; m < M; m += 4)
        {
            px0 = (const ae_f24x2 *)(px3);
            px1 = (const ae_f24x2 *)((f24 *)px0 + N);
            px2 = (const ae_f24x2 *)((f24 *)px1 + N);
            px3 = (const ae_f24x2 *)((f24 *)px2 + N);
            py = (const ae_f24x2 *)y;

            ACC00 = ACC01 = ACC10 = ACC11 =
                ACC20 = ACC21 = ACC30 = ACC31 = AE_ZERO64();

            /* Innermost loop: compute 2 values for 4 rows */
            __Pragma("loop_count factor=2");
            for (n = 0; n < (N >> 1); n++)
            {
                /* load data from 'x' */
                AE_L32X2F24_IP(vx0, px0, sizeof(ae_f24x2));
                AE_L32X2F24_IP(vx1, px1, sizeof(ae_f24x2));
                AE_L32X2F24_IP(vx2, px2, sizeof(ae_f24x2));
                AE_L32X2F24_IP(vx3, px3, sizeof(ae_f24x2));
                /* load data from 'y' */
                AE_L32X2F24_XP(vt0, py, P*sizeof(f24));
                AE_L32X2F24_XP(vt1, py, P*sizeof(f24));
                vy0 = AE_SEL24_HH(vt0, vt1);
                vy1 = AE_SEL24_LL(vt0, vt1);
                /* perform multiplications */
                AE_MULAAFD24_HH_LL(ACC00, vx0, vy0);
                AE_MULAAFD24_HH_LL(ACC01, vx0, vy1);
                AE_MULAAFD24_HH_LL(ACC10, vx1, vy0);
                AE_MULAAFD24_HH_LL(ACC11, vx1, vy1);
                AE_MULAAFD24_HH_LL(ACC20, vx2, vy0);
                AE_MULAAFD24_HH_LL(ACC21, vx2, vy1);
                AE_MULAAFD24_HH_LL(ACC30, vx3, vy0);
                AE_MULAAFD24_HH_LL(ACC31, vx3, vy1);
            }
            /* format values */
            ACC00 = AE_SLAA64(ACC00, lsh);
            ACC01 = AE_SLAA64(ACC01, lsh);
            ACC10 = AE_SLAA64(ACC10, lsh);
            ACC11 = AE_SLAA64(ACC11, lsh);
            ACC20 = AE_SLAA64(ACC20, lsh);
            ACC21 = AE_SLAA64(ACC21, lsh);
            ACC30 = AE_SLAA64(ACC30, lsh);
            ACC31 = AE_SLAA64(ACC31, lsh);
            vz0 = AE_ROUND24X2F48SASYM(ACC00, ACC01);
            vz1 = AE_ROUND24X2F48SASYM(ACC10, ACC11);
            vz2 = AE_ROUND24X2F48SASYM(ACC20, ACC21);
            vz3 = AE_ROUND24X2F48SASYM(ACC30, ACC31);
            /* save values */
            AE_S32X2F24_XP(vz0, pz, P*sizeof(f24));
            AE_S32X2F24_XP(vz1, pz, P*sizeof(f24));
            AE_S32X2F24_XP(vz2, pz, P*sizeof(f24));
            AE_S32X2F24_XP(vz3, pz, P*sizeof(f24));

        }
        y += 2;
        z += 2;
    }
} /* mtx_mpy24x24_fast() */
