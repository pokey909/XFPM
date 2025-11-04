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
    Real data circular cross-correlation, floating point
    C code optimized for HiFi3
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Circular Correlation
  Estimates the circular cross-correlation between vectors x (of length N) 
  and y (of length M)  resulting in vector r of length N. It is a similar 
  to correlation but x is read in opposite direction.

  Precision: 
  16x16     16x16-bit data, 16-bit outputs
  24x24     24x24-bit data, 24-bit outputs
  32x16     32x16-bit data, 32-bit outputs
  32x32     32x32-bit data, 32-bit outputs
  f         floating point 


  Input:
  x[N]      input data Q15, Q31 or floating point
  y[M]      input data Q15, Q31 or floating point
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
#include "common_fpu.h"

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void,fir_xcorrf,(     float32_t * restrict r,
               const float32_t * restrict x,
               const float32_t * restrict y,
               int N, int M ))
#elif (HAVE_VFPU)

void fir_xcorrf(     float32_t * restrict r,
               const float32_t * restrict x,
               const float32_t * restrict y,
               int N, int M )
{
  //
  // Circular cross-correlation algorithm:
  //
  //   r[n] = sum( x[mod(n+m,N)]*y[m] )
  //        m=0..M-1
  //
  //   where n = 0..N-1
  //
    xtfloatx2 A01,A23,A45,A67,X01,X12,X23,X34,Y0,Y1;
    const xtfloatx2 * restrict pX;
    const xtfloatx2 * restrict pX1;
    const ae_int32* restrict pY;
    xtfloatx2 * restrict pR;
    ae_valign ax1;

    int n, m;

    NASSERT(r);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(r,8);
    NASSERT_ALIGN(x,8);
    NASSERT_ALIGN(y,8);
    NASSERT((M>0) && M%4==0);
    NASSERT((N>0) && N%4==0);

    pR=(      xtfloatx2*)r;
    /* set circular buffer boundaries */
    WUR_AE_CBEGIN0( (uintptr_t)( x + 0 ) );
    WUR_AE_CEND0  ( (uintptr_t)( x + N ) );
    for ( n=0; n<N; n+=4,x+=4 )
    {
        A01=A23=A45=A67=(xtfloat)0.f;
        pX=(const xtfloatx2*)(x);
        pX1=(const xtfloatx2*)(x+1);
        pY=(const ae_int32*)y;
        AE_LA32X2POS_PC(ax1,castxcc(const ae_int32x2,pX1));
        XT_LSX2XC(X01,pX,sizeof(X01));
        XT_LASX2IC(X12,ax1,pX1);
        __Pragma("loop_count min=1");
        for ( m=0; m<(M>>2); m++ )
        {
            ae_int32x2 tmp;
            XT_LSX2XC(X23,pX,sizeof(X23));
            XT_LASX2IC(X34,ax1,pX1);
            AE_L32_IP(tmp,pY,sizeof(float32_t)); Y0=XT_AE_MOVXTFLOATX2_FROMINT32X2(tmp);
            AE_L32_IP(tmp,pY,sizeof(float32_t)); Y1=XT_AE_MOVXTFLOATX2_FROMINT32X2(tmp);
            XT_MADD_SX2(A01,X01,Y0);
            XT_MADD_SX2(A23,X23,Y0);
            XT_MADD_SX2(A45,X12,Y1);
            XT_MADD_SX2(A67,X34,Y1);
            X01=X23;
            X12=X34;

            XT_LSX2XC(X23,pX,sizeof(X23));
            XT_LASX2IC(X34,ax1,pX1);
            AE_L32_IP(tmp,pY,sizeof(float32_t)); Y0=XT_AE_MOVXTFLOATX2_FROMINT32X2(tmp);
            AE_L32_IP(tmp,pY,sizeof(float32_t)); Y1=XT_AE_MOVXTFLOATX2_FROMINT32X2(tmp);
            XT_MADD_SX2(A01,X01,Y0);
            XT_MADD_SX2(A23,X23,Y0);
            XT_MADD_SX2(A45,X12,Y1);
            XT_MADD_SX2(A67,X34,Y1);
            X01=X23;
            X12=X34;
        }
        A01=XT_ADD_SX2(A01,A45);
        A23=XT_ADD_SX2(A23,A67);
        XT_SSX2IP(A01,pR,sizeof(A01));
        XT_SSX2IP(A23,pR,sizeof(A23));
    }
} /* fir_xcorrf() */
#else 
// for scalar FPU
void fir_xcorrf(     float32_t * restrict r,
               const float32_t * restrict x,
               const float32_t * restrict y,
               int N, int M )
{
  /*
  * Circular cross-correlation algorithm:
  *
  *   r[n] = sum( x[mod(n+m,N)]*y[m] )
  *        m=0..M-1
  *
  *   where n = 0..N-1
  */
  xtfloat A0,A1,A2,A3,X0,X1,X2,X3,X4,Y0,Y1;
  const xtfloat * restrict pX;
  const xtfloat* restrict pY;
  xtfloat * restrict pR;
  int n, m;

  NASSERT(r);
  NASSERT(x);
  NASSERT(y);
  NASSERT_ALIGN(r,8);
  NASSERT_ALIGN(x,8);
  NASSERT_ALIGN(y,8);
  NASSERT(M>0 && M%4==0);
  NASSERT(N>0 && N%4==0);
  if (N <= 0 || M <= 0) return;
  pR=(      xtfloat*)r;
  /* set circular buffer boundaries */
  WUR_AE_CBEGIN0( (uintptr_t)( x + 0 ) );
  WUR_AE_CEND0  ( (uintptr_t)( x + N ) );
  for ( n=0; n<N; n+=4,x+=4 )
  {
    A0=A1=A2=A3=XT_CONST_S(0);
    pX = (const xtfloat*)(x);
    pY = (const xtfloat*)y;
    XT_LSXC(X0, pX, 4);
    XT_LSXC(X1, pX, 4);
    XT_LSXC(X2, pX, 4);
    XT_LSXC(X3, pX, 4);
    XT_LSXC(X4, pX, 4);
    __Pragma("loop_count min=1");
    for ( m=0; m<(M>>1); m++ )
    {
      XT_LSIP(Y0, pY, 4);
      XT_LSIP(Y1, pY, 4);
      XT_MADD_S(A0, X0, Y0);
      XT_MADD_S(A1, X1, Y0);
      XT_MADD_S(A2, X2, Y0);
      XT_MADD_S(A3, X3, Y0);
      XT_MADD_S(A0, X1, Y1);
      XT_MADD_S(A1, X2, Y1);
      XT_MADD_S(A2, X3, Y1);
      XT_MADD_S(A3, X4, Y1);
      X0 = X2;
      X1 = X3;
      X2 = X4;
      XT_LSXC(X3, pX, 4);
      XT_LSXC(X4, pX, 4);
    }
    XT_SSIP(A0,pR,4);
    XT_SSIP(A1,pR,4);
    XT_SSIP(A2,pR,4);
    XT_SSIP(A3,pR,4);
  }
} /* fir_xcorrf() */
#endif
