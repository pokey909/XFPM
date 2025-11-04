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
 * inputs: x  int16 (vector)
 * output: zr int16 scalar value
 * kernel can handle saturation cases
 * restrictions: input should be 8 byte aligned
 */
#include "common.h"
#define __HiFi3_3z_NDSP__ 1

#if __HiFi3_3z_NDSP__
int16_t vec_sum16x16 (const int16_t* restrict x, int N)
{
	const ae_int16x4* restrict px;
    	ae_int16x4 xt, zt=0;
    	ae_int16   zf;

	int i;
	NASSERT(x);
	NASSERT_ALIGN(x, 16);
	NASSERT((N&3)==0);
	if(N<=0) return 0;

	px=(const ae_int16x4 *)x;
	for ( i=0; i<(N>>2); i++ )
	{
		AE_L16X4_IP(xt,(ae_int16x4*)px,sizeof(ae_int16x4));
		zt = AE_ADD16S(zt,xt);
	}
	zf = AE_INT16X4_RADD(zt);

	for ( i=0; i<(N&3); i++ )
	{
		ae_int16x4 tmp;
		AE_L16_IP(tmp,castxcc(ae_int16,px),sizeof(int16_t));
		zf = AE_ADD16S((ae_int16x4)zf,(ae_int16x4)tmp);
	}

	return zf;
}

#elif __REFC_OOB__
#include "math.h"
int16_t vec_sum16x16 (const int16_t* restrict x, int N)
{
	int i;
	int16_t sumt=0;

	NASSERT(x);
	NASSERT_ALIGN(x, 16);

	if(N<0) return 0.0f;
	#pragma aligned (x, 16)

	for (i=0; i<N; i++)
	{
		sumt += x[i];
	}
	return sumt;
}
#endif
