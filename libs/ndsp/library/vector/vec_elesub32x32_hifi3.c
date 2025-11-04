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
/*  Element-wise subtraction of two vectors
	Input 	: x - int32 vector
	Input 	: y - int32 vector
	Output 	: z - int32 vector
	restrictions: inputs should be 8 byte aligned
				: x and y should not overlap
*/
#include "common.h"
#define __HiFi3_3z_NDSP__ 1
#if __HiFi3_3z_NDSP__
void vec_elesub32x32 (int32_t* restrict z, const int32_t* restrict x, int32_t* restrict y, int N)
{
	const ae_int32x2* restrict px;
	const ae_int32x2* restrict py;
	ae_int32x2* restrict pz;
	ae_int32x2 xt, yt, zt;
	int i;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    if(N<=0) return;

    px=(const ae_int32x2 *)x;
    py=(const ae_int32x2 *)y;
    pz=(ae_int32x2 *)z;

	for (i=0; i<(N>>1); i++ )
	{
		AE_L32X2_IP(xt,(ae_int32x2*)px,sizeof(ae_int32x2));
		AE_L32X2_IP(yt,(ae_int32x2*)py,sizeof(ae_int32x2));
		zt = AE_SUB32S(xt,yt);
		AE_S32X2_IP(zt,(ae_int32x2*)pz,sizeof(ae_int32x2));
	}
	if(N&1)
	{
		AE_L32_IP(xt,(ae_int32*)px,sizeof(ae_int32));
		AE_L32_IP(yt,(ae_int32*)py,sizeof(ae_int32));
		zt = AE_SUB32S(xt,yt);
		AE_S32_L_IP(zt,(ae_int32*)pz,sizeof(ae_int32));
	}
}

#elif __REFC_OOB__
#include "math.h"
void vec_elesub32x32 (int32_t* restrict z, int32_t* restrict x, int32_t* restrict y, int N)
{
	int i;

    NASSERT(x);
    NASSERT(z);
    if(N<=0) return;

	#pragma aligned (x, 16)
	#pragma aligned (z, 16)

    for ( i=0; i<N; i++ )
    {
    	z[i] = x[i] - y[i];
    }
}
#endif
