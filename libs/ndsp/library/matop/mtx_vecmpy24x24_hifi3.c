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
void mtx_vecmpy24x24 (  f24* restrict z,
                  const f24* restrict x,
                  const f24* restrict y,
                  int M, int N, int lsh)

{
  const ae_f24x2 * restrict px0 = (const ae_f24x2 *)x;
  const ae_f24x2 * restrict px1 = (const ae_f24x2 *)x;
  const ae_f24x2 * restrict py  = (const ae_f24x2 *)y;
        ae_f24   * restrict pz  = (      ae_f24   *)z;
  ae_f24x2  vx0,vx1,vy0;
  ae_f24x2  vz0,vz1;
  ae_f32x2  vxf0,vxf1,vyf0;
  ae_f64    ACC0,ACC1;
  ae_valign al_px0,al_px1,al_py;
  int m,n;
  NASSERT(lsh >= -31 && lsh <= 31);
  if (N<=0 || M<=0)/* exceptional situation */
  {
      for (m=0; m<M; m++) z[m]=0;
      return;
  }
  
  px1 = (const ae_f24x2 *)(x);
  pz  = (      ae_f24   *)z;

  for(m=0; m<(M>>1); m++)
  {
    px0 = px1;
    px1 = (const ae_f24x2 *)((f24 *)px0+N);
    py  = (const ae_f24x2 *)(y);
    al_px0 = AE_LA64_PP(px0);
    al_px1 = AE_LA64_PP(px1);
    al_py  = AE_LA64_PP(py);
    
    ACC0 = ACC1 = AE_ZERO();
    
    for(n=0; n<(N>>1); n++)
    {
      AE_LA32X2F24_IP(vx0, al_px0, px0);
      AE_LA32X2F24_IP(vx1, al_px1, px1);
      AE_LA32X2F24_IP(vy0, al_py , py);
      AE_MULAAFD24_HH_LL(ACC0,vy0,vx0);
      AE_MULAAFD24_HH_LL(ACC1,vy0,vx1);
    }
    if (N&1)
    {
      AE_L32F24_IP(vx0, castxcc(ae_f24,px0), 4);
      AE_L32F24_IP(vx1, castxcc(ae_f24,px1), 4);
      AE_L32F24_IP(vy0, castxcc(ae_f24,py ), 4);
      vxf0 = (vx0);
      vxf1 = (vx1);
      vyf0 = (vy0);
      AE_MULAF32S_LL(ACC0,vxf0,vyf0);
      AE_MULAF32S_LL(ACC1,vxf1,vyf0);
    }

    ACC0 = AE_SLAA64(ACC0,lsh);
    ACC1 = AE_SLAA64(ACC1,lsh);
    vz0 = AE_ROUNDSP24Q48ASYM(ACC0);
    vz1 = AE_ROUNDSP24Q48ASYM(ACC1);
    AE_S32F24_L_IP(vz0,pz,4);
    AE_S32F24_L_IP(vz1,pz,4);
  }
  
  if (M&1)
  {
    px0 = px1;
    py = (const ae_f24x2 *)y;
    al_px0 = AE_LA64_PP(px0);
    al_py  = AE_LA64_PP(py);
    
    ACC0 = AE_ZERO();
    
    for(n=0; n<(N>>1); n++)
    {
      AE_LA32X2F24_IP(vx0, al_px0, px0);
      AE_LA32X2F24_IP(vy0, al_py , py);
      AE_MULAAFD24_HH_LL(ACC0,vy0,vx0);
    }
    if (N&1)
    {
      AE_L32F24_IP(vx0, castxcc(ae_f24,px0), 4);
      AE_L32F24_IP(vy0, castxcc(ae_f24,py), 4);
      vxf0 = (vx0);
      vyf0 = (vy0);
      AE_MULAF32S_LL(ACC0,vxf0,vyf0);
    }

    ACC0 = AE_SLAA64(ACC0,lsh);
    vz0 = AE_ROUNDSP24Q48ASYM(ACC0);
    AE_S32F24_L_IP(vz0,pz,4);
  }
} /* mtx_vecmpy24x24() */
