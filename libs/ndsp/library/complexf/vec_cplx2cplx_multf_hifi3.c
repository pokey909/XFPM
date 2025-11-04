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
/* complex to complex multiplication of single precision floating point data
 * input  x,y : complex float
 * output z   : complex float
 */
#include "common_fpu.h"

#define __HiFi3_3z_NDSP__ 1

#if !HAVE_VFPU && !HAVE_FPU
DISCARD_FUN(void, vec_cplx2cplx_multf, (complex_float* restrict z, complex_float* restrict x, complex_float* restrict y, int N))
#elif HAVE_VFPU
#if __HiFi3_3z_NDSP__
void vec_cplx2cplx_multf (complex_float* restrict z, complex_float* restrict x, complex_float* restrict y, int N)
{
    const xtfloatx2 * restrict pX;
    const xtfloatx2 * restrict pY;
          xtfloatx2 * restrict pZ;
          xtfloatx2 z0;
          xtfloatx2 x0,y0;
    int i;

    NASSERT(z);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(z, 8);
    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(y, 8);
    if (N<=0) return;

    pX = (const xtfloatx2 *)x;
    pY = (const xtfloatx2 *)y;
    pZ = (xtfloatx2 *)z;

    for(i=0;i<N;i++)
    {
        AE_LSX2IP(x0, pX, sizeof(xtfloatx2));
        AE_LSX2IP(y0, pY, sizeof(xtfloatx2));
        z0 = MULC_S(x0, y0);
        AE_SSX2IP(z0, pZ, sizeof(xtfloatx2));
    }
}
#elif __REFC_OOB__   
void vec_cplx2cplx_multf (complex_float* restrict z, complex_float* restrict x, complex_float* restrict y, int N)
{
	int i;
	float32_t *xptr, *yptr, *zptr;

	xptr = (float32_t*) x;
	yptr = (float32_t*) y;
	zptr = (float32_t*) z;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);

    if(N<=0) return;

	#pragma aligned (xptr, 16)
	#pragma aligned (yptr, 16)
	#pragma aligned (zptr, 16)

    for ( i=0; i<(N*2); i+=2 )
    {
    	zptr[i] 	= (xptr[i] * yptr[i])   - (xptr[i+1] * yptr[i+1]);
    	zptr[i+1] 	= (xptr[i] * yptr[i+1]) + (xptr[i+1] * yptr[i]);
    }
}
#endif
#endif
