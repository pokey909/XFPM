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
 * vector complex conjugate operation for 16 bit integer data
 * Input x : int16 complex vector  (real and imag parts are interleaved)
 * Output z: int16 complex vector  (real and imag parts are interleaved)
 *
 */
#include "common.h"


#define __HiFi3_3z_NDSP__ 1
#define __REFC_OOB__ 0 //OOB cannot handle edge cases

#if (__HiFi3_3z_NDSP__)
void vec_cplxconj16x16(complex16_t* restrict z, const complex16_t* restrict x, int N)
{
   const ae_int16x4* restrict pX;
	      ae_int16x4* restrict pZ;
	      ae_int16x4 x0r, x0i, t0;
	    int i;

	    NASSERT(x);
	    NASSERT(z);
	    NASSERT_ALIGN(x, 8);
	    NASSERT_ALIGN(z, 8);

	    if(N<=0) return;

	    pX = (const ae_int16x4*)x;
	    pZ = (ae_int16x4*) z;

	    for ( i=0; i<(N); i++ )
	    {
	    	AE_L16_IP(x0r,castxcc(ae_int16,pX),sizeof(int16_t));
	    	AE_L16_IP(x0i,castxcc(ae_int16,pX),sizeof(int16_t));
	        t0 = AE_NEG16S(x0i);
	        AE_S16_0_IP(x0r,castxcc(ae_int16,pZ),sizeof(int16_t));
	        AE_S16_0_IP(t0,castxcc(ae_int16,pZ),sizeof(int16_t));
	    }
}

#else
void vec_cplxconj16x16(complex16_t* restrict z, const complex16_t* restrict x, int N)
{
	int i;
	int16_t *xptr, *zptr;

	NASSERT(x);
	NASSERT_ALIGN(x, 16);
	if(N<=0) return;

	xptr = (int16_t*) x;
	zptr = (int16_t*) z;

	#pragma aligned (x, 16)

    for ( i=0; i<(N*2); i+=2 )
    {
    	zptr[i] = xptr[i];
    	zptr[i+1] = -1 * xptr[i+1];
    }
	return;
}
#endif