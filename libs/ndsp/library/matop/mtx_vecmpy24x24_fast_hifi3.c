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
void mtx_vecmpy24x24_fast (  f24* restrict z,
               const f24* restrict x,
               const f24* restrict y,
               int M, int N, int lsh)

{
  const ae_f24x2 * restrict px0;
  const ae_f24x2 * restrict px1;
  const ae_f24x2 * restrict px2;
  const ae_f24x2 * restrict px3;
  const ae_f24x2 * restrict py;
        ae_f24x2 * restrict pz;
  
  ae_f24x2 vx0,vx1,vx2,vx3,vy0;
  ae_f24x2 vz0,vz1;
  ae_f64   ACC0,ACC1,ACC2,ACC3;
  ae_int32x2 zz;
  ae_valign z_align;
  int m, n;

  zz = AE_MOVI(0);
  NASSERT(x);
  NASSERT(y);
  NASSERT(z);
  NASSERT_ALIGN8(x);
  NASSERT_ALIGN8(y);
  ASSERT((M&1)==0);
  ASSERT((N&1)==0);
  NASSERT(lsh >= -31 && lsh <= 31);
  px3 = (const ae_f24x2 *)x;
  pz  = (      ae_f24x2 *)z;
  z_align = AE_ZALIGN64();
  #ifdef COMPILER_XTENSA
    #pragma concurrent
  #endif
  for(m=M; m>0; m-=4)
  {
    px0 = px3;
    px1 = (const ae_f24x2 *)((f24 *)px0+N);
    px2 = (const ae_f24x2 *)((f24 *)px1+N);
    px3 = (const ae_f24x2 *)((f24 *)px2+N);
    py  = (const ae_f24x2 *)y;

    ACC0 = ACC1 = ACC2 = ACC3 = AE_MOVF64_FROMINT32X2(zz);

    for(n=(N>>1); n>0; n--)
    {
      AE_L32X2F24_IP(vx0, px0, 8);
      AE_L32X2F24_IP(vx1, px1, 8);
      AE_L32X2F24_IP(vx2, px2, 8);
      AE_L32X2F24_IP(vx3, px3, 8);
      AE_L32X2F24_IP(vy0, py , 8);
      AE_MULAAFD24_HH_LL(ACC0,vx0,vy0);
      AE_MULAAFD24_HH_LL(ACC1,vx1,vy0);
      AE_MULAAFD24_HH_LL(ACC2,vx2,vy0);
      AE_MULAAFD24_HH_LL(ACC3,vx3,vy0);
    }
    ACC0 = AE_SLAA64(ACC0,lsh);
    ACC1 = AE_SLAA64(ACC1,lsh);
    ACC2 = AE_SLAA64(ACC2,lsh);
    ACC3 = AE_SLAA64(ACC3,lsh);
    vz0 = AE_ROUND24X2F48SASYM(ACC0, ACC1);
    vz1 = AE_ROUND24X2F48SASYM(ACC2, ACC3);
    AE_SA32X2F24_IP(vz0, z_align, pz);
    AE_SA32X2F24_IP(vz1, z_align, pz);
  }
  AE_SA64POS_FP(z_align, pz);
} /* mtx_vecmpy24x24_fast() */
