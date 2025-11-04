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
/*  Variance of a vector
	Input 	: x - int32 vector
	Output 	: z - int32 scalar
*/
#include "common.h"
#include "NatureDSP_Signal_math.h"
#define __HiFi3_3z_NDSP__ 1

#if __HiFi3_3z_NDSP__
int32_t vec_var32x32( const int32_t* restrict x, int N )
{
	const ae_int32x2* restrict px;
		ae_int32x2 xt, meanx2, zt=0;
		ae_int32  res;
		int i, rt;
		ae_int64x2 t =ae_int32x2_rtor_ae_int64x2 ((AE_MOV32(0)));
		ae_int64x2 sum64x2 = ae_int32x2_rtor_ae_int64x2 ((AE_MOV32(0)));
		long long sum64;

		NASSERT(x);
		if(N<2) return 0;

		px=(const ae_int32x2 *)x;

		// Calculate Mean
		for ( i=0; i<(N>>1); i++ )
		{
			AE_L32X2_IP(xt,(ae_int32x2*)px,sizeof(ae_int32x2));
			zt = AE_ADD32S(xt,zt);
		}
		ae_int32x2 ztLH;
		ztLH = AE_SEL32_LH(zt,zt); // swap
		zt = AE_ADD32S(zt, ztLH);
		if(N&1)
		{
			AE_L32_IP(xt,(ae_int32*)px,sizeof(ae_int32));
			zt = AE_ADD32S(xt,zt);
		}
		rt = AE_MOVAD32_L(zt);
		rt = rt / N;

		// Calculate Variance
		meanx2 = AE_MOVINT32X2_FROMF32(rt);
		px=(const ae_int32x2 *)x;
		for (i=0; i<(N>>1); i++ )
		{
			AE_L32X2_IP(xt,(ae_int32x2*)px,sizeof(ae_int32x2));
			zt = AE_SUB32S(xt,meanx2);
			AE_MULA32X2_vector (sum64x2,zt,zt);
		}
		if(N&1)
		{
			AE_L32_IP(xt,(ae_int32*)px,sizeof(ae_int32));
			zt = AE_SUB32S(xt,meanx2);
			zt = AE_SEL32_HH(zt, AE_MOV32(0));
			AE_MULA32X2_vector(t,zt,zt);
		}
		sum64x2 = AE_ADD64X2_vector(sum64x2, t);
		sum64 = AE_INT64X2_RADD(sum64x2);
		res = (int32_t)(sum64/ (N-1) );
	return res;
}


#elif __REFC_OOB__
#include "math.h"
int32_t vec_var32x32 (const int32_t* restrict x, int N)
{
	int i;
	int32_t tmp, sumt=0, mean, var;

	NASSERT(x);
	NASSERT_ALIGN(x, 16);

	if(N<=0) return 0.0f;
	#pragma aligned (x, 16)

	/*  Compute the Mean of all elements */
	for (i=0; i<N; i++)
	{
		sumt += x[i];
	}
	mean = (int32_t) sumt / N;

	/*  Compute  variance */
	sumt = 0;
	for (i=0; i<N; i++)
	{
		tmp = x[i] - mean;
		sumt += tmp * tmp;
	}
	var = (int32_t) sumt / (N-1);

	return var;
}

#endif
