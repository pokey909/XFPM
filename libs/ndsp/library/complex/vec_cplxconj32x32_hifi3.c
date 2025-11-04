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
 * vector complex conjugate operation for 32 bit integer data
 * Input x : int32 complex vector  (real and imag parts are interleaved)
 * Output z: int32 complex vector  (real and imag parts are interleaved)
 *
 */
#include "common.h"

#define __HiFi3_3z_NDSP__ 1
#define __REFC_OOB__ 0 //OOB cannot handle edge cases
#if __HiFi3_3z_NDSP__
void vec_cplxconj32x32(complex32_t* restrict z, const complex32_t* restrict x, int N)
{
	const ae_int32x2* restrict px;
		  ae_int32x2* restrict pz;
	ae_int32x2 xt, yt=0,  zt;
	int i;

	NASSERT(x);
	NASSERT_ALIGN(x, 16);
	NASSERT(z);
	NASSERT_ALIGN(z, 16);
	if(N<=0) return;

	px = (const ae_int32x2*)x;
	pz = (ae_int32x2*) z;

	for ( i=0; i<(N); i++ )
	{
		AE_L32X2_IP(xt,px,sizeof(ae_int32x2));
		zt = AE_ADDSUB32S(yt, xt);
		AE_S32X2_IP(zt,pz,sizeof(ae_int32x2));
	}
	return;
}

#else
void vec_cplxconj32x32(complex32_t* restrict z, const complex32_t* restrict x, int N)
{
	int i;
	int32_t *xptr, *zptr;

	NASSERT(x);
	NASSERT_ALIGN(x, 16);

	if(N<0) return;

	xptr = (int32_t*) x;
	zptr = (int32_t*) z;

	#pragma aligned (x, 16)

    for ( i=0; i<(N*2); i+=2 )
    {
    	zptr[i] = xptr[i];
    	zptr[i+1] = -1 * xptr[i+1];
    }
	return;
}
#endif
