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
#define __PERF_OPT__ 1
/*  Vector Variance: Variance of entire Input vector
	Input 	: x - float vector
	Output 	: z - scalar float
	Restrictions: Inputs should be 8 byte aligned
*/
#include "common_fpu.h"
#define __HiFi3_3z_NDSP__ 1


#if !HAVE_VFPU && !HAVE_FPU
DISCARD_FUN(float32_t, vec_varf, ( const float32_t * restrict x, int N ))
#elif HAVE_VFPU

#if __HiFi3_3z_NDSP__

#if __PERF_OPT__
float32_t vec_varf ( const float32_t * restrict x, int N )
{
	xtfloatx2 xtt;
	const xtfloatx2 * px;
	xtfloatx2 z0, zt=0.0f, st=0.0f;
	xtfloat xt, y=0, n, z;
	int i;

    NASSERT(x);
    NASSERT_ALIGN(x, 16);
    if(N<=1) return XT_CONST_S(0);

	px = (xtfloatx2*)x;

    // Calculate Mean
	for(i=0;i<(N>>1);i++)
	{
		XT_LSX2IP(xtt, px, sizeof(xtfloatx2));
		zt = ADD_SX2(xtt,zt);
	}
	y = RADD_SX2(zt);

	if(N&1)
	{
		XT_LSIP(xt, (xtfloat*)px, sizeof(xtfloat));
		y = XT_ADD_S(y,xt);
	}

	n = (xtfloat)(N);
	z = DIV_S( y, n);
	y = z;

	// Calculate Variance
	zt = AE_MOVXTFLOATX2_FROMXTFLOAT(z);
	px = (xtfloatx2*)x;

    for (i=0;i<(N>>1);i++)
    {
    	XT_LSX2IP(xtt, px, sizeof(xtfloatx2));
    	z0 = XT_SUB_SX2(xtt,zt);
    	MADD_SX2(st,z0,z0);
    }
    z = RADD_SX2(st);

	if(N&1)
	{
	   xtfloat xf,zf;
	   AE_LSIP(xf,castxcc(xtfloat,px),sizeof(xtfloat));
	   zf = SUB_S(xf,y);
	   MADD_S(z,zf,zf);
	}
    z = DIV_S( z, (float)(N-1)); 	// Normalize with N-1
    return z;

}

#else	// HIGH ACCURACY
float32_t vec_varf ( const float32_t * restrict x, int N )
{
	xtfloatx2 xtt;
	const xtfloatx2 * px;
	xtfloatx2 z0, zt=0.0f, st=0.0f;
	xtfloat xt, y=0, z;
	int i;

    NASSERT(x);
    NASSERT_ALIGN(x, 16);
	if(N<=1) return XT_CONST_S(0);

	px = (xtfloatx2*)x;

	// Calculate Mean
	for ( i=0; i<(N); i++ )
	{
		XT_LSIP(xt, (xtfloat*)px, sizeof(xtfloat));
		y = XT_ADD_S(y, xt);
	}
	y = DIV_S( y, N);

	// Calculate Variance
	zt = AE_MOVXTFLOATX2_FROMXTFLOAT(y);
	px = (xtfloatx2*)x;

    for (i=0;i<(N>>1);i++)
    {
    	XT_LSX2IP(xtt, px, sizeof(xtfloatx2));
    	z0 = XT_SUB_SX2(xtt,zt);
    	MADD_SX2(st,z0,z0);
    }
    z = RADD_SX2(st);

	if(N&1)
	{
	   xtfloat xf,zf;
	   AE_LSIP(xf,castxcc(xtfloat,px),sizeof(xtfloat));
	   zf = SUB_S(xf,y);
	   MADD_S(z,zf,zf);
	}
    z = DIV_S( z, (float)(N-1)); 	// Normalize with N-1
    return z;

}
#endif // __PERF_OPT__

#elif __REFC_OOB__
#include "math.h"
float32_t vec_varf (const float32_t * restrict x, int N)
{
	int i;
	float32_t tmp, sumt=0, mean, var;

	NASSERT(x);
	NASSERT_ALIGN(x, 16);

	if(N<=0) return 0.0f;
	#pragma aligned (x, 16)

	/*  Compute the Mean of all elements */
	if (N==1) return 0.0f;
	for (i=0; i<N; i++)
	{
		sumt += x[i];
	}
	mean = sumt / (float)N;

	/*  Compute  variance */
	sumt = 0;
	for (i=0; i<N; i++)
	{
		tmp = x[i] - mean;
		sumt += tmp * tmp;
	}
	var = sumt / (float)(N-1);

	return var;
}

#endif

#endif
