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
/* Calculates sum of Input vector
 * inputs: x  int32 (vector)
 * output: resF int32 scalar value
 * restrictions: input should be 8 byte aligned
 */
#include "common.h"
#define __HiFi3_3z_NDSP__ 1

#if __HiFi3_3z_NDSP__
int32_t vec_sum32x32 (const int32_t* restrict x, int N)
{
	const ae_int32x2* restrict px;
	ae_int32x2 x0,zt=0;
	ae_int32 resF;
	int i;
	NASSERT(x);
	if(N<=0) return 0;
	px=(const ae_int32x2 *)x;

	for ( i=0; i<(N>>1); i++ )
	{
		AE_L32X2_IP(x0,(ae_int32x2*)px,sizeof(ae_int32x2));
		zt = AE_ADD32S(x0,zt);
	}
#if XCHAL_HAVE_HIFI3Z
	zt = AE_SUBADD32S_HL_LH(zt, zt);
#else
	ae_int32x2 zh;
	zh = AE_SEL32_HH(zt,zt);
	zt = AE_SUBADD32S(zt, zh); //only lower part contains valid result
#endif
	if(N&1)
	{
		AE_L32_IP(x0,(ae_int32*)px,sizeof(ae_int32));
		zt = AE_ADD32S(x0,zt);
	}
	resF = AE_MOVAD32_L(zt);//Only lower element contains correct result

	return resF;
}

#elif __REFC_OOB__
#include "math.h"
int32_t vec_sum32x32 (const int32_t* restrict x, int N)
{
	int i;
	int32_t sumt=0;

	NASSERT(x);
	NASSERT_ALIGN(x, 16);

	if(N<0) return 0.0f;
	#pragma aligned (x, 16)

	for (i=0; i<N; i++)
	{
		sumt += (int32_t)x[i];
	}
	return (int32_t)sumt;

}

#endif
