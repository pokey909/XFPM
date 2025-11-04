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
void mtx_vecmpy16x16 (  int16_t* restrict z,
                  const int16_t* restrict x,
                  const int16_t* restrict y,
                  int M, int N, int lsh)
{
#ifndef AE_MULAAAAQ16
  const ae_int16x4 * restrict px0;
  const ae_int16x4 * restrict px1;
  const ae_int16x4 * restrict py;
        ae_int16   * restrict pz;
  ae_int32x2 va0h,t1,t2;
  ae_int32x2 y1,y2;
  ae_int16x4 out0,out1,vx0,vx1;
  ae_int16x4 vy0,vy1,vy2;
  ae_int64 ACC0, ACC1;
  ae_valign x0_align,x1_align,y_align;
  int m,n;
  
  NASSERT(x);
  NASSERT(y);
  NASSERT(z);
  NASSERT(lsh >= -15 && lsh <= 15);
  if (N<=0 || M<=0)/* exceptional situation */
  {
    for (m=0; m<M; m++) z[m]=0;
    return;
  }
  
  pz = (ae_int16 *)z;

  /* Compute multiplication by 2 values per iteration */
  for(m=0; m<(M&~1); m+=2)
  {
    px0 = (const ae_int16x4 *)(x+m*N);
    px1 = (const ae_int16x4 *)((int16_t *)px0+N);
    py  = (const ae_int16x4 *)(y);
    x0_align = AE_LA64_PP(px0);
    x1_align = AE_LA64_PP(px1);
    y_align  = AE_LA64_PP(py);
    ACC0 = ACC1 = AE_ZERO64();

    for(n=0; n<(N>>2); n++)
    {
      AE_LA16X4_IP(vx0, x0_align, px0);
      AE_LA16X4_IP(vx1, x1_align, px1);
      AE_LA16X4_IP(vy0, y_align , py );
      y1 = AE_SEXT32X2D16_32(vy0);
      y2 = AE_SEXT32X2D16_10(vy0);
      AE_MULAAD32X16_H3_L2(ACC0,y1,vx0);
      AE_MULAAD32X16_H1_L0(ACC0,y2,vx0);
      AE_MULAAD32X16_H3_L2(ACC1,y1,vx1);
      AE_MULAAD32X16_H1_L0(ACC1,y2,vx1);
    }
    
    switch(N&3)
    {
      case 1:
        AE_L16_IP(vy0, castxcc(ae_int16,py ), 2);
        AE_L16_IP(vx0, castxcc(ae_int16,px0), 2);
        AE_L16_IP(vx1, castxcc(ae_int16,px1), 2);
        y1 = AE_SEXT32X2D16_32(vy0);
        AE_MULA32X16_L0(ACC0,y1,vx0);
        AE_MULA32X16_L0(ACC1,y1,vx1);
        break;
      case 2:
        AE_L16_IP(vy0, castxcc(ae_int16,py), 2);
        AE_L16_IP(vy1, castxcc(ae_int16,py), 2);
        AE_LA16X4_IP(vx0, x0_align, px0);
        AE_LA16X4_IP(vx1, x1_align, px1);
        y1 = AE_SEXT32X2D16_32(vy0);
        y2 = AE_SEXT32X2D16_32(vy1);
        AE_MULA32X16_L3(ACC0,y1,vx0);
        AE_MULA32X16_L2(ACC0,y2,vx0);
        AE_MULA32X16_L3(ACC1,y1,vx1);
        AE_MULA32X16_L2(ACC1,y2,vx1);
        break;
      case 3:
        AE_L16_IP(vy0,castxcc(ae_int16,py), 2);
        AE_L16_IP(vy1,castxcc(ae_int16,py), 2);
        AE_L16_IP(vy2,castxcc(ae_int16,py), 2);
        AE_LA16X4_IP(vx0, x0_align, px0);
        AE_LA16X4_IP(vx1, x1_align, px1);
        t1 = AE_SEXT32X2D16_32(vy0);
        t2 = AE_SEXT32X2D16_32(vy1);
        y1 = AE_SEL32_HH(t1,t2);
        y2 = AE_SEXT32X2D16_32(vy2);
        AE_MULAAD32X16_H3_L2(ACC0,y1,vx0);
        AE_MULA32X16_L1(ACC0,y2,vx0);
        AE_MULAAD32X16_H3_L2(ACC1,y1,vx1);
        AE_MULA32X16_L1(ACC1,y2,vx1);
        break;
      default: break;
    }

    /* Format values to Q15 and store */
    va0h = AE_TRUNCA32X2F64S(ACC1, ACC0, 33+lsh);
    /* convert to 16-bit data */
    out0 = AE_ROUND16X4F32SASYM(va0h, va0h);
    out1 = AE_SEL16_4321(out0, out0);

    AE_S16_0_IP(out0,pz,2);
    AE_S16_0_IP(out1,pz,2);
  }

  /* Compute last odd value */
  if (M&1)
  {
    px0 = (const ae_int16x4 *)(x+(M-1)*N);
    py  = (const ae_int16x4 *)(y);
    x0_align = AE_LA64_PP(px0);
    y_align  = AE_LA64_PP(py);
    ACC0 = AE_ZERO64();

    for(n=0; n<(N>>2); n++)
    {
      AE_LA16X4_IP(vx0, x0_align, px0);
      AE_LA16X4_IP(vy0, y_align , py );
      y1 = AE_SEXT32X2D16_32(vy0);
      y2 = AE_SEXT32X2D16_10(vy0);
      AE_MULAAD32X16_H3_L2(ACC0,y1,vx0);
      AE_MULAAD32X16_H1_L0(ACC0,y2,vx0);
    }
    
    switch(N&3)
    {
      case 1:
        AE_L16_IP(vx0, castxcc(ae_int16,px0), 2);
        AE_L16_IP(vy0, castxcc(ae_int16,py ), 2);
        y1 = AE_SEXT32X2D16_32(vy0);
        AE_MULA32X16_L0(ACC0,y1,vx0);
        break;
      case 2:
        AE_L16_IP(vy0, castxcc(ae_int16,py), 2);
        AE_L16_IP(vy1, castxcc(ae_int16,py), 2);
        AE_LA16X4_IP(vx0, x0_align, px0); 
        y1 = AE_SEXT32X2D16_32(vy0);
        y2 = AE_SEXT32X2D16_32(vy1);
        AE_MULA32X16_L3(ACC0,y1,vx0);
        AE_MULA32X16_L2(ACC0,y2,vx0);
        break;
      case 3:
        AE_L16_IP(vy0,castxcc(ae_int16,py), 2);
        AE_L16_IP(vy1,castxcc(ae_int16,py), 2);
        AE_L16_IP(vy2,castxcc(ae_int16,py), 2);
        AE_LA16X4_IP(vx0, x0_align, px0);
        t1 = AE_SEXT32X2D16_32(vy0);
        t2 = AE_SEXT32X2D16_32(vy1); 
        y1 = AE_SEL32_HH(t1,t2);
        y2 = AE_SEXT32X2D16_32(vy2); 
        AE_MULAAD32X16_H3_L2(ACC0,y1,vx0);
        AE_MULA32X16_L1(ACC0,y2,vx0);
        break;
      default: break;
    }
    
    /* Format values to Q15 and store */
    va0h = AE_TRUNCA32X2F64S(ACC0, ACC0, 33+lsh);
    /* convert to 16-bit data */
    out0 = AE_ROUND16X4F32SASYM(va0h, va0h);

    AE_S16_0_IP(out0,pz,2);
  }
#else
    int m, n;
    xtbool4 mask;
    const ae_int16x4      *  restrict px0;
    const ae_int16x4      *  restrict px1;
    const ae_int16x4      *  restrict py;
    ae_int16        *  restrict pz = (ae_int16       *)z;
    ae_int32x2 a0;
    ae_int16x4 t, x0, x1, y0;
    ae_int64 B0, B1;
    ae_valign ax0, ax1, ay;
    NASSERT(lsh >= -15 && lsh <= 15);
    if (N <= 0 || M <= 0)    /* exceptional situation */
    {
        for (m = 0; m < M; m++) z[m] = 0;
        return;
    }
    /* mask for last 1...3 values from y */
    mask = (1 << (4 - (N & 3))) - 1;
    /* process by pair of rows */
    for (m = 0; m < (M&~1); m += 2)
    {
        B0 = B1 = AE_ZERO64();
        px0 = (const ae_int16x4 *)x;
        px1 = (const ae_int16x4 *)(x + N);
        py = (const ae_int16x4 *)y;
        ax0 = AE_LA64_PP(px0);
        ax1 = AE_LA64_PP(px1);
        ay = AE_LA64_PP(py);
        for (n = 0; n < (N&~3); n += 4)
        {
            AE_LA16X4_IP(x0, ax0, px0);
            AE_LA16X4_IP(x1, ax1, px1);
            AE_LA16X4_IP(y0, ay, py);
            AE_MULAAAAQ16(B0, x0, y0);
            AE_MULAAAAQ16(B1, x1, y0);
        }
        if (N & 3)
        {
            AE_LA16X4_IP(x0, ax0, px0);
            AE_LA16X4_IP(x1, ax1, px1);
            AE_LA16X4_IP(y0, ay, py);
            AE_MOVT16X4(y0, AE_ZERO16(), mask);
            AE_MULAAAAQ16(B0, x0, y0);
            AE_MULAAAAQ16(B1, x1, y0);
        }
        a0 = AE_TRUNCA32X2F64S(B0, B1, lsh + 33);
        t = AE_ROUND16X4F32SASYM(a0, a0);
        AE_S16_0_IP(AE_SHORTSWAP(t), pz, 2);
        AE_S16_0_IP(t, pz, 2);
        x += N * 2; /* next 2 rows */
    }
    if (M & 1)
    {
        B0 = AE_ZERO64();
        px0 = (const ae_int16x4 *)x;
        py = (const ae_int16x4 *)y;
        ax0 = AE_LA64_PP(px0);
        ay = AE_LA64_PP(py);
        for (n = 0; n < (N&~3); n += 4)
        {
            AE_LA16X4_IP(x0, ax0, px0);
            AE_LA16X4_IP(y0, ay, py);
            AE_MULAAAAQ16(B0, x0, y0);
        }
        if (N & 3)
        {
            AE_LA16X4_IP(x0, ax0, px0);
            AE_LA16X4_IP(y0, ay, py);
            AE_MOVT16X4(y0, AE_ZERO16(), mask);
            AE_MULAAAAQ16(B0, x0, y0);
        }
        a0 = AE_TRUNCA32X2F64S(B0, B0, lsh + 33);
        t = AE_ROUND16X4F32SASYM(a0, a0);
        AE_S16_0_IP(t, pz, 2);
    }
#endif
} /* mtx_vecmpy16x16() */
