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

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_matop.h"
/* Common helper macros. */
#include "common.h"
#include "common_fpu.h"

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void ,mtx_vecmpyf,( float32_t * z, const float32_t * x,  const float32_t * y, int M, int N ))
#elif (HAVE_VFPU)
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
void mtx_vecmpyf( float32_t * z, const float32_t * x,  const float32_t * y, int M, int N )
{
  const xtfloatx2 *restrict px0;
  const xtfloatx2 *restrict px1;
  const xtfloatx2 *restrict py;
        xtfloatx2 *restrict pz;
  xtfloatx2 y0, y1, z0;
  xtfloatx2 x10, x11, x00, x01;
  xtfloatx2 acc00, acc01, acc10, acc11;
  xtfloat x0_, y0_, z0_;
  ae_valign vx0, vx1, vy, vz;
  int m, n, N_;

  NASSERT(x);
  NASSERT(y);
  NASSERT(z);
  NASSERT((z != x) && (z != y));

  if (M <= 0)    return;
  /* If N<=0 then clear output vector and return */
  if (N <= 0)
  {
    pz = (xtfloatx2 *)z;
    vz = AE_ZALIGN64();
    z0 = (xtfloatx2)0.0f;
    for (n = 0; n < (M>>1); n++)
    {
      XT_SASX2IP(z0, vz, pz);
    }
    XT_SASX2POSFP(vz, pz);
    if (M&1) XT_SSI(z0, (xtfloat *)pz, 0);
    return;
  }

  N_ = N & (~3);
  pz = (xtfloatx2 *)z;
  vz = AE_ZALIGN64();
  
  /* Compute by 2 values */
  /* Perform all multiplications for each value except last N%4 */
  for (m = 0; m < (M>>1); m++)
  {
    int m2 = m<<1;
    px0 = (const xtfloatx2 *)(x+m2*N);
    px1 = (const xtfloatx2 *)((float32_t *)px0+N);
    py  = (const xtfloatx2 *)(y);

    vx0 = XT_LASX2PP(px0);
    vx1 = XT_LASX2PP(px1);
    vy  = XT_LASX2PP(py );
    acc00 = acc01 = acc10 = acc11 = (xtfloatx2)0.0f;

    for (n = 0; n < (N>>2); n++)
    {
      XT_LASX2IP(x00, vx0, px0);
      XT_LASX2IP(x01, vx0, px0);
      XT_LASX2IP(x10, vx1, px1);
      XT_LASX2IP(x11, vx1, px1);
      XT_LASX2IP( y0,  vy,  py);
      XT_LASX2IP( y1,  vy,  py);

      XT_MADD_SX2(acc00, x00, y0);
      XT_MADD_SX2(acc01, x01, y1);
      XT_MADD_SX2(acc10, x10, y0);
      XT_MADD_SX2(acc11, x11, y1);
    }
    acc00 = acc00 + acc01;
    acc10 = acc10 + acc11;
    y0 = XT_SEL32_HL_SX2(acc00, acc10);
    y1 = XT_SEL32_LH_SX2(acc00, acc10);
    z0 = y0 + y1;

    XT_SASX2IP(z0, vz, pz);
  }
  XT_SASX2POSFP(vz, pz);
  /* Compute last (if odd) output element */
  if (M&1)
  {
    px0 = (const xtfloatx2 *)(x+M*N-N);
    py  = (const xtfloatx2 *)(y);

    vx0 = XT_LASX2PP(px0);
    vy  = XT_LASX2PP(py );
    acc00 = acc01 = (xtfloatx2)0.0f;
    
    for (n = 0; n < (N>>2); n++)
    {
      XT_LASX2IP(x00, vx0, px0);
      XT_LASX2IP(x01, vx0, px0);
      XT_LASX2IP( y0,  vy,  py);
      XT_LASX2IP( y1,  vy,  py);

      XT_MADD_SX2(acc00, x00, y0);
      XT_MADD_SX2(acc01, x01, y1);
    }
    acc00 = acc00 + acc01;

    z0_ = XT_RADD_SX2(acc00);

    XT_SSIP(z0_, castxcc(xtfloat,pz), sizeof(float32_t));
  }
  /* Update output vector with remaining N%4 multiplications */
  py  = (const xtfloatx2 *)(y+N_);
  for (n = N_; n < N; n++)
  {
    px0 = (const xtfloatx2 *)(x+n);
    pz  = (      xtfloatx2 *)(z);

    XT_LSIP(y0_, castxcc(const xtfloat,py), sizeof(float32_t));
    
    for (m = 0; m < M; m++)
    {
      XT_LSXP(x0_, castxcc(const xtfloat,px0), N*sizeof(float32_t));
      z0_ = XT_LSI((xtfloat *)pz, 0);
      XT_MADD_S(z0_, x0_, y0_);
      XT_SSIP(z0_, castxcc(xtfloat,pz), sizeof(float32_t));
    }
  }
} /* mtx_vecmpyf() */
#elif (HAVE_FPU)
void mtx_vecmpyf( float32_t * z, const float32_t * x,  const float32_t * y, int M, int N )
{
    int m, n;
          xtfloat * restrict pZ;
    const xtfloat * restrict pX;
    const xtfloat * restrict pY;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    NASSERT((z != x) && (z != y));

    if (M <= 0)    return;
    /* If N<=0 then clear output vector and return */
    pZ=(xtfloat*)z;
    if (N<0)
    {
        for ( m=0; m<M; m++) XT_SSIP(XT_CONST_S(0),pZ,sizeof(xtfloat));
        return;
    }
    for ( m=0; m<(M&~3); m+=4)
    {
        xtfloat a0,a1,a2,a3,x0,x1,x2,x3,y0;
        a0=a1=a2=a3=XT_CONST_S(0);
        pX=(const xtfloat*)(x+m*N);
        pY=(const xtfloat*)(y);
        for ( n=0; n<N; n++ )
        {
            XT_LSIP(y0,pY,sizeof(xtfloat));
            x1=XT_LSX(pX,1*N*sizeof(xtfloat));
            x2=XT_LSX(pX,2*N*sizeof(xtfloat));
            x3=XT_LSX(pX,3*N*sizeof(xtfloat));
            XT_LSIP(x0,pX,sizeof(xtfloat));
            XT_MADD_S(a0,x0,y0);
            XT_MADD_S(a1,x1,y0);
            XT_MADD_S(a2,x2,y0);
            XT_MADD_S(a3,x3,y0);
        }
        XT_SSIP(a0,pZ,sizeof(xtfloat));
        XT_SSIP(a1,pZ,sizeof(xtfloat));
        XT_SSIP(a2,pZ,sizeof(xtfloat));
        XT_SSIP(a3,pZ,sizeof(xtfloat));
    }
    if (M&2)
    {
        xtfloat a0,a1,x0,x1,y0;
        a0=a1=XT_CONST_S(0);
        pX=(const xtfloat*)(x+m*N);
        pY=(const xtfloat*)(y);
        for ( n=0; n<N; n++ )
        {
            XT_LSIP(y0,pY,sizeof(xtfloat));
            x1=XT_LSX(pX,1*N*sizeof(xtfloat));
            XT_LSIP(x0,pX,sizeof(xtfloat));
            XT_MADD_S(a0,x0,y0);
            XT_MADD_S(a1,x1,y0);
        }
        XT_SSIP(a0,pZ,sizeof(xtfloat));
        XT_SSIP(a1,pZ,sizeof(xtfloat));
        m+=2;
    }
    if (M&1)
    {
        xtfloat a0,x0,y0;
        a0=XT_CONST_S(0);
        pX=(const xtfloat*)(x+m*N);
        pY=(const xtfloat*)(y);
        for ( n=0; n<N; n++ )
        {
            XT_LSIP(y0,pY,sizeof(xtfloat));
            XT_LSIP(x0,pX,sizeof(xtfloat));
            XT_MADD_S(a0,x0,y0);
        }
        XT_SSIP(a0,pZ,sizeof(xtfloat));
    }
}

#endif
