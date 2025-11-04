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
 * Copyright (c) 2024 by Cadence Design Systems, Inc.  ALL RIGHTS RESERVED.
 * These coded instructions, statements, and computer programs are the
 * copyrighted works and confidential proprietary information of
 * Cadence Design Systems Inc.  They may be adapted and modified by bona fide
 * purchasers for internal use, but neither the original nor any adapted
 * or modified version may be disclosed or distributed to third parties
 * in any manner, medium, or form, in whole or in part, without the prior
 * written consent of Cadence Design Systems Inc.  This software and its
 * derivatives are to be executed solely on products incorporating a Cadence
 * Design Systems processor.
 */
/*  RMS of a vector
	Input 	: x - int32 vector
	Output 	: z - int32
*/
/*
 *  Conditions :  N is less than 1024
 *  			  x[n] -> No restriction
 *  			  72 bit accumulation
 *  			  if x[n] values are very small and N value is also very small
 */

#include "common.h"
#include "NatureDSP_Signal_math.h"


#define __HiFi3_3z_NDSP__ 1


#if __HiFi3_3z_NDSP__
int32_t vec_rms32x32( const int32_t* restrict x, int N )
{
#if XCHAL_HAVE_HIFI3Z
	const ae_int32x2* restrict px;
	ae_int32x2 xt;
	ae_int64 result=0;
	ae_int64 y0=0;
	int i, rt;

	NASSERT(x);
	if(N<=0) return 0;
	px=(const ae_int32x2 *)x;

	for ( i=0; i<(N>>1); i++ )
	{
		AE_L32X2_IP(xt,(ae_int32x2*)px,sizeof(ae_int32x2));
		AE_MULAAD32_HH_LL(y0,xt,xt);
	}
	if(N&1)
	{
		AE_L32_IP(xt,(ae_int32*)px,sizeof(ae_int32));
		AE_MULA32_HH(y0,xt,xt);
	}
	result = ((long long)y0/ N);
	result = result << 1; 	// To convert from q62 to q63
	rt = scl_sqrt64x32(result);
	return rt;
#else //HiFI3
	const ae_int32x2* restrict px;
	ae_int32x2 x0;
	ae_int64 sum64=0;
	ae_int64x2 sum64x2 = ae_int32x2_rtor_ae_int64x2 ((AE_MOV32(0)));
	int i, rt;

	NASSERT(x);
	if(N<=0) return 0;
	px=(const ae_int32x2 *)x;

	for ( i=0; i<(N>>1); i++ )
	{
		AE_L32X2_IP(x0,(ae_int32x2*)px,sizeof(ae_int32x2));
		AE_MULA32X2_vector (sum64x2,x0,x0);
	}
	if(N&1)
	{
		AE_L32_IP(x0,(ae_int32*)px,sizeof(ae_int32));
		x0 = AE_SEL32_HH(x0, AE_MOV32(0));
		AE_MULA32X2_vector(sum64x2,x0,x0);
	}
	sum64 = AE_INT64X2_RADD(sum64x2);
	sum64 = ((long long)sum64/ N);
	sum64 = sum64 << 1; 	// To convert from q62 to q63
	rt = scl_sqrt64x32(sum64);
	return rt;
#endif
}

#elif __REFC_OOB__
#include "math.h"
int32_t vec_rms32x32( const int32_t* restrict x, int N )
{
	int i;
	//int32_t sumt=0,rms;
	long long tmp, sumt=0;
	int32_t rms;

	NASSERT(x);
	NASSERT_ALIGN(x, 16);

	if(N<=0) return 0.0f;
	#pragma aligned (x, 16)

	for (i=0; i<N; i++)
	{
		tmp = (long long) x[i];
		sumt += ( tmp * tmp );
		//sumt += x[i] * x[i];
	}
	tmp = sumt / N;
	tmp = (tmp << 1);
	rms = (int32_t) scl_sqrt64x32(tmp);
	// rms = rms << 31;

	return rms;
}
#endif
