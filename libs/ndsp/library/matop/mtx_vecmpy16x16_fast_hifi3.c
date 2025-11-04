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
    Matrix by Vector Multiply
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/
/* Signal Processing Library API. */
#include "NatureDSP_Signal_matop.h"
#include "NatureDSP_types.h"
#include "common.h"

/*-------------------------------------------------------------------------
  Matrix by Vector Multiply
  These functions compute the expression z = 2^lsh * x * y for the matrices 
  x and vector y. 
  NOTE: lsh factor is not relevant for floating point routines.

  Two versions of functions available: regular version (mtx_vecmpy32x32, 
  mtx_vecmpy24x24, mtx_vecmpy16x16,mtx_vecmpyf) with arbitrary arguments 
  and faster version (mtx_vecmpy32x32_fast, mtx_vecmpy24x24_fast, 
  mtx_vecmpy16x16_fast, mtx_vecmpyf_fast) that apply some restrictions.

  Precision: 
  32x32 32-bit input, 32-bit output
  24x24 24-bit input, 24-bit output
  16x16 16-bit input, 16-bit output
  f     floating point

  Input:
  x[M*N] input matrix,Q31,Q15 or floating point
  y[N]   input vector,Q31,Q15 or floating point
  M      number of rows in matrix x
  N      number of columns in matrix x
  lsh    additional left shift(applied to the fixed-
         point functions only) 
  Output:
  z[M]   output vector,Q31,Q15 or floating point

  Restriction:
  For regular routines (mtx_vecmpy32x32, mtx_vecmpy24x24,
                        mtx_vecmpy16x16, mtx_vecmpyf)
  x,y,z should not overlap

  For faster routines (mtx_vecmpy32x32_fast, mtx_vecmpy24x24_fast,
                       mtx_vecmpy16x16_fast, mtx_vecmpyf_fast)
  x,y,z   should not overlap
  x,y     aligned on 8-byte boundary
  N, M    multiples of 4
  lsh     should be in range:
          -31...31 for mtx_vecmpy32x32, mtx_vecmpy32x32_fast,
                   mtx_vecmpy24x24, mtx_vecmpy24x24_fast
          -15...15 for mtx_vecmpy16x16, mtx_vecmpy16x16_fast  
-------------------------------------------------------------------------*/
void mtx_vecmpy16x16_fast (  int16_t* restrict z,
                       const int16_t* restrict x,
                       const int16_t* restrict y,
                       int M, int N, int lsh)
{
#ifndef AE_MULAAAAQ16
  const ae_int16x4 * restrict px0;
  const ae_int16x4 * restrict px1;
  const ae_int16x4 * restrict px2;
  const ae_int16x4 * restrict px3;
  const ae_int16x4 * restrict py;
        ae_int16   * restrict pz;
  int m,n;
  ae_valign al_z;

  NASSERT(x);
  NASSERT(y);
  NASSERT(z);
  NASSERT_ALIGN8(x);
  NASSERT_ALIGN8(y);
  ASSERT((N&3)==0);
  ASSERT((M&3)==0);
  NASSERT(lsh >= -15 && lsh <= 15);
  pz  = (      ae_int16   *)(z);
  px3 = (const ae_int16x4 *)(x);
  al_z = AE_ZALIGN64();

  if (N<=0)
  {
      ae_int16x4 zero;
      zero = AE_ZERO16();
      for (m = 0; m < (M>>2); m++)
      {
          AE_SA16X4_IP(zero, al_z, castxcc(ae_int16x4,pz));
      }
      AE_SA64POS_FP(al_z, pz);
      return;
  }

  /* Compute multiplication by 4 values per iteration */
  for(m=0; m<M; m+=4)
  {
    ae_int32x2 va0h, va1h;
    ae_int16x4 out0;
    ae_int64   ACC0, ACC1, ACC2, ACC3;
    
    px0 = (const ae_int16x4 *)(px3);
    px1 = (const ae_int16x4 *)((int16_t *)px0+N);
    px2 = (const ae_int16x4 *)((int16_t *)px1+N);
    px3 = (const ae_int16x4 *)((int16_t *)px2+N);
    py  = (const ae_int16x4 *)y;
    ACC0 = ACC1 = ACC2 = ACC3 = AE_ZERO64();

    __Pragma("loop_count min=1")
    for(n=0; n<N>>2; n++)
    {
      ae_int16x4 vxh0,vxh1,vxh2,vxh3,vyh;
      ae_int32x2 y0,y1;
      /* load input data */
      AE_L16X4_IP(vxh0, px0, 8);
      AE_L16X4_IP(vxh1, px1, 8);
      AE_L16X4_IP(vxh2, px2, 8);
      AE_L16X4_IP(vxh3, px3, 8);
      AE_L16X4_IP(vyh , py , 8);
      /* convert data from 'y' to 32 bit */
      y0 = AE_SEXT32X2D16_32(vyh);
      y1 = AE_SEXT32X2D16_10(vyh);
      /* perform multiplications */
      AE_MULAAD32X16_H3_L2(ACC0, y0, vxh0);
      AE_MULAAD32X16_H1_L0(ACC0, y1, vxh0);
      AE_MULAAD32X16_H3_L2(ACC1, y0, vxh1);
      AE_MULAAD32X16_H1_L0(ACC1, y1, vxh1);
      AE_MULAAD32X16_H3_L2(ACC2, y0, vxh2);
      AE_MULAAD32X16_H1_L0(ACC2, y1, vxh2);
      AE_MULAAD32X16_H3_L2(ACC3, y0, vxh3);
      AE_MULAAD32X16_H1_L0(ACC3, y1, vxh3);
    }
    /* truncate accumulators to 32-bit data w/ saturation */
    va1h = AE_TRUNCA32X2F64S(ACC2, ACC3, 33+lsh);
    va0h = AE_TRUNCA32X2F64S(ACC0, ACC1, 33+lsh);
    /* convert to 16-bit data */
    out0 = AE_ROUND16X4F32SASYM(va0h, va1h);
    /* Store values */
    AE_SA16X4_IP(out0, al_z, castxcc(ae_int16x4,pz));
  }
  AE_SA64POS_FP(al_z, pz);
#else
    int m, n;
    const ae_int16x4 * restrict px = (const ae_int16x4 *)x;
    const ae_int16x4 * restrict py = (const ae_int16x4 *)y;
          ae_int16x4 * restrict pZ = (      ae_int16x4 *)z;
    ae_valign az;

    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT((N & 3) == 0);
    NASSERT((M & 3) == 0);
    NASSERT(lsh >= -15 && lsh <= 15);
    az = AE_ZALIGN64();
    for (m = 0; m < M; m += 4)
    {
        ae_int32x2 a0, a1;
        ae_int64 B0, B1, B2, B3;
        ae_int16x4 t;

        py = (const ae_int16x4 *)y;
        B0 = B1 = B2 = B3 = AE_ZERO64();
        for (n = 0; n < (N >> 2); n++)
        {
            ae_int16x4 x0, x1, x2, x3, y0;
            x3 = AE_L16X4_X(px, 3 * N * 2);
            x2 = AE_L16X4_X(px, 2 * N * 2);
            x1 = AE_L16X4_X(px, N * 2);
            AE_L16X4_IP(x0, px, 8);
            AE_L16X4_IP(y0, py, 8);
            AE_MULAAAAQ16(B0, x0, y0);
            AE_MULAAAAQ16(B1, x1, y0);
            AE_MULAAAAQ16(B2, x2, y0);
            AE_MULAAAAQ16(B3, x3, y0);
        }
        a0 = AE_TRUNCA32X2F64S(B0, B1, lsh + 33);
        a1 = AE_TRUNCA32X2F64S(B2, B3, lsh + 33);
        t = AE_ROUND16X4F32SASYM(a0, a1);
        AE_SA16X4_IP(t, az, pZ);
        px += 3 * N / 4;
    }
    AE_SA64POS_FP(az, pZ);
#endif
} /* mtx_vecmpy16x16_fast() */
