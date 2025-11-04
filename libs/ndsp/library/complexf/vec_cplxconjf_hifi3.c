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
 * Vector Operations
 * vector complex conjugate operation for floating point data
 * Input  x: single precision float complex vector  (real and imag parts are interleaved)
 * Output z: single precision float complex vector  (real and imag parts are interleaved)
 *
 */
#include "common_fpu.h"

#if !HAVE_VFPU && !HAVE_FPU
DISCARD_FUN(void, vec_cplx_Conjf, (complex_float * restrict z, const complex_float * restrict x, int N))
#elif HAVE_VFPU

#define __HiFi3_3z_NDSP__ 1

#if __HiFi3_3z_NDSP__
void vec_cplx_Conjf(complex_float * restrict r,
					const complex_float * restrict x,
					int N)
{
	const xtfloatx2 *pX;
   		  xtfloatx2 *pR;
   	int i;
   	xtfloatx2 xt, rt;

	NASSERT(r);
	NASSERT(x);
	NASSERT_ALIGN(r, 16);
	NASSERT_ALIGN(x, 16);
	if(N<0) return;

	pX = (xtfloatx2*)x;
	pR = (xtfloatx2*)r;

	for(i=0;i<(N); i++)
	{
		AE_LSX2IP(xt, pX, sizeof(xtfloatx2));
		rt = XT_CONJC_S(xt);
		AE_SSX2IP(rt, pR, sizeof(xtfloatx2));
	}
}

#else
void vec_cplx_Conjf(complex_float * restrict z, const complex_float * restrict x, int N)
{
	int i;
	float32_t *xptr, *zptr;

	NASSERT(x);
	NASSERT_ALIGN(x, 16);

	if(N<0) return;

	xptr = (float32_t*) x;
	zptr = (float32_t*) z;

	#pragma aligned (x, 16)

    for ( i=0; i<(N*2); i+=2 )
    {
    	zptr[i] = xptr[i];
    	zptr[i+1] = -1.0f * xptr[i+1];
    }
	return;
}

#endif

#endif
