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
 * Copyright (c) 2021 by Cadence Design Systems, Inc.  ALL RIGHTS RESERVED.
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
/*  Element-wise multiplication of two vectors
	Input 	: x,y - float vector
	Output 	: z - float vector
	Restrictions: x,y and z should not overlap
 * 		        : Inputs should be 8 byte aligned
*/
#include "common_fpu.h"

#if !HAVE_VFPU && !HAVE_FPU
DISCARD_FUN(void, vec_elemultf, (float32_t* restrict z, float32_t* restrict x, float32_t* restrict y, int N))
#elif HAVE_VFPU

#define __HiFi3_3z_NDSP__ 1

#if __HiFi3_3z_NDSP__
void vec_elemultf (float32_t* restrict z, float32_t* restrict x, float32_t* restrict y, int N)
{
	const xtfloatx2* restrict px;
	const xtfloatx2* restrict py;
	xtfloatx2* restrict pz;
	xtfloatx2 xtt, ytt, ztt;
	xtfloat xt,yt,zt;
	int i;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    if(N<=0) return;

    px = (xtfloatx2*)x;
    py = (xtfloatx2*)y;
    pz = (xtfloatx2*)z;

    for ( i=0; i<(N>>1); i++ )
    {
    	XT_LSX2IP(xtt, px, sizeof(xtfloatx2));
    	XT_LSX2IP(ytt, py, sizeof(xtfloatx2));
    	ztt = XT_MUL_SX2(xtt,ytt);
    	XT_SSX2IP(ztt, pz, sizeof(xtfloatx2));
    }
   if(N&1)
    {
    	AE_LSIP(xt, (xtfloat*)px, sizeof(xtfloat));
    	AE_LSIP(yt, (xtfloat*)py, sizeof(xtfloat));
    	zt=XT_MUL_S(xt,yt);
    	AE_SSIP(zt, (xtfloat*)pz, sizeof(xtfloat));
    }
}

#elif __REFC_OOB__
#include "math.h"
void vec_elemultf (float32_t* restrict z, float32_t* restrict x, float32_t* restrict y, int N)
{
	int i;

    NASSERT(x);
    NASSERT(z);
    if(N<=0) return;

	#pragma aligned (x, 16)
	#pragma aligned (z, 16)

    for ( i=0; i<N; i++ )
    {
    	z[i] = x[i] * y[i];
    }
}

#endif	// __HiFi3_3z_NDSP__

#endif	//	HAVE_VFPU
