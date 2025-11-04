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
 * Vector Operations vec_cplx2real_multvf:
 * vector complex to real multiplicaiton (scaling) operation for single precision float data
 * Each complex number is scaled by a unique number
 * Input x : single precision float vector complex (real and imag parts are interleaved)
 * Input y : single precision float vector real
 * Output z: single precision float complex (real and imag parts are interleaved)
 *
 */#include "common.h"
#include "common_fpu.h"

#define  __HiFi3_3z_NDSP__ 1

#if !HAVE_VFPU && !HAVE_FPU
DISCARD_FUN(void, vec_cplx2real_multvf, (complex_float* restrict z, complex_float* restrict x, float32_t* restrict y, int N))
DISCARD_FUN(void, vec_cplx2real_multsf, (complex_float* restrict z, complex_float* restrict x, float32_t y, int N))
#elif HAVE_VFPU

#if __HiFi3_3z_NDSP__
void vec_cplx2real_multvf (complex_float* restrict z, complex_float* restrict x, float32_t* restrict y, int N)
{
const xtfloatx2 *pX;
const xtfloat *pY;
      xtfloatx2 *pZ;
   int i;
   xtfloatx2 x00, y00;
   xtfloat ytmp;
   xtfloatx2 z00;

   NASSERT(x);
   NASSERT(y);
   NASSERT(z);
   NASSERT_ALIGN(z, 8);
   NASSERT_ALIGN(x, 8);

   NASSERT(N > 0);

   if(N<=0) return;

   pX = (xtfloatx2*)x;
   pY =  (xtfloat*)y;
   pZ = (xtfloatx2*)z;
   for(i=0; i<N; i++)
   {
       AE_LSX2IP(x00, pX, sizeof(xtfloatx2));   // Load x
       AE_LSIP(ytmp, pY, sizeof(xtfloat));      // Load y
       y00 = AE_MOVXTFLOATX2_FROMXTFLOAT(ytmp);
       z00 = MUL_SX2(y00, x00);
       AE_SSX2IP(z00,pZ, sizeof(xtfloatx2));
   }
}
 /*
  * vec_cplx2real_multsf:
  * vector complex to real multiplicaiton (scaling) operation for single precision float data
  * Entire complex vector is scaled by a fixed number
  * Input x : single precision float vector complex (real and imag parts are interleaved)
  * Input y : single precision float scalar real
  * Output z: single precision float vector complex (real and imag parts are interleaved)
  *
  */
void vec_cplx2real_multsf (complex_float* restrict z, complex_float* restrict x, float32_t y, int N)
{
	const xtfloatx2 *pX;
	      xtfloatx2 *pZ;
	   int i;
	   xtfloatx2 x00, y00;
	   xtfloat ytmp=y;
	   xtfloatx2 z00;

	   NASSERT(x);
	   NASSERT(y);
	   NASSERT(z);
	   NASSERT_ALIGN(z, 8);
	   NASSERT_ALIGN(x, 8);

	   NASSERT(N > 0);

	   if(N<=0) return;

	   pX = (xtfloatx2*)x;
	   pZ = (xtfloatx2*)z;
	   for(i=0; i<N; i++)
	   {
	       AE_LSX2IP(x00, pX, sizeof(xtfloatx2));   // Load x
	       y00 = AE_MOVXTFLOATX2_FROMXTFLOAT(ytmp);
	       z00 = MUL_SX2(y00, x00);
	       AE_SSX2IP(z00,pZ, sizeof(xtfloatx2));
	   }
}

#elif __REFC_OOB__
void vec_cplx2real_multvf (complex_float* restrict z, complex_float* restrict x, float32_t* restrict y, int N)
{
    int i,j;
    float32_t *xptr, *yptr, *zptr;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);

    if(N<=0) return;

    xptr = (float32_t*) x;
    yptr = (float32_t*) y;
    zptr = (float32_t*) z;

    #pragma aligned (x, 16)
    #pragma aligned (y, 16)
    #pragma aligned (z, 16)

    for(i=0,j=0; i<(N*2); i+=2,j++)
    {
        zptr[i] = xptr[i] * yptr[j];
        zptr[i+1] = xptr[i+1] * yptr[j];
    }
}

void vec_cplx2real_multsf (complex_float* restrict z, complex_float* restrict x, float32_t y, int N)
{
	int i;
	float32_t *xptr, *zptr;

    NASSERT(x);
    NASSERT(z);

    if(N<=0) return;

	xptr = (float32_t*) x;
	zptr = (float32_t*) z;

	#pragma aligned (x, 16)
	#pragma aligned (z, 16)

    for(i=0; i<N; i+=2)
    {
    	zptr[i] = xptr[i] * y;
    	zptr[i+1] = xptr[i+1] * y;
    }
}
#endif 

#endif

