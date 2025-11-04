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
/*  Vector RMS: RMS of entire Input vector
	Input 	: x - float vector
	Output 	: z - scalar float
	Restrictions: Inputs should be 8 byte aligned
*/
#include "common_fpu.h"

#define __HiFi3_3z_NDSP__ 1

#if !HAVE_VFPU && !HAVE_FPU
DISCARD_FUN(float32_t, vec_rmsf, ( float32_t * restrict x, int N))
#elif HAVE_VFPU

#if __HiFi3_3z_NDSP__
float32_t vec_rmsf (  const float32_t * restrict x, int N )
{
	xtfloatx2 x0, x1;
	xtfloatx2 * pX;
	xtfloatx2 zt = 0.0f;
	xtfloat n, z;
	int i;
	pX = (xtfloatx2*)x;
    NASSERT(x);
    NASSERT_ALIGN(x, 16);
	if(N<=0) return XT_CONST_S(0);

	for(i=0;i<N>>2;i++)
	{
		XT_LSX2IP(x0,pX,sizeof(xtfloatx2));
		XT_LSX2IP(x1,pX,sizeof(xtfloatx2));
		MADD_SX2(zt,x0,x0);
		MADD_SX2(zt,x1,x1);
	}
	if(N&2)
	{
		AE_LSX2IP(x0,castxcc(xtfloatx2, pX),sizeof(xtfloatx2));
		MADD_SX2(zt,x0,x0);
	}
	z = RADD_SX2(zt);
	if(N&1)
	{
		xtfloat y;
		AE_LSIP(y,castxcc(xtfloat,pX),sizeof(xtfloat));
		MADD_S(z,y,y);
	}
	n = (xtfloat)(N);
	z = DIV_S( z, n);
	z = SQRT_S(z);

	return z;
}

#elif __REFC_OOB__
#include "math.h"
float32_t vec_rmsf (const float32_t * restrict x, int N)
{
	int i;
	float32_t sumt=0,rms;

	NASSERT(x);
	NASSERT_ALIGN(x, 16);

	if(N<=0) return 0.0f;
	#pragma aligned (x, 16)

	for (i=0; i<N; i++)
	{
		sumt += x[i] * x[i];
	}
	rms = sumt / N;
	rms = sqrt(rms);

	return rms;
}

#endif
#endif
