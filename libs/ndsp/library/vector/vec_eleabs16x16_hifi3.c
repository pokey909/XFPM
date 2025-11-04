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
/* Calculates element-wise absolute of Input vector
 * inputs: x  int16 (vector)
 * output: z  int16 (vector)
 * kernel  can handle saturation cases
 * restrictions: input should be 8 byte aligned
 */
#include "common.h"
#define __HiFi3_3z_NDSP__ 1

#if __HiFi3_3z_NDSP__
void vec_eleabs16x16 (const int16_t* restrict x, int16_t* restrict z,  int N)
{
	const ae_int16x4* restrict px;
		  ae_int16x4* restrict pz;
    ae_int16x4 xt, zt;

	int i;
	NASSERT(x);
	NASSERT(z);
	if(N<=0) return;

	px=(const ae_int16x4 *)x;
	pz=(ae_int16x4 *)z;	

	for ( i=0; i<(N>>2); i++ )
	{
		AE_L16X4_IP(xt,(ae_int16x4*)px,sizeof(ae_int16x4));
		zt = AE_ABS16S(xt);
		AE_S16X4_IP(zt,(ae_int16x4*)pz,sizeof(ae_int16x4));
	}

	for ( i=0; i<(N&3); i++ )
	{
		AE_L16_IP(xt,castxcc(ae_int16,px),sizeof(int16_t));
		zt = AE_ABS16S(xt);
		AE_S16_0_IP(zt,castxcc(ae_int16,pz),sizeof(int16_t));
	}
}

#elif __REFC_OOB__
#include "math.h"
#include <stdlib.h>
void vec_eleabs16x16 (const int16_t* restrict x, int16_t* restrict z,  int N)
{
	int i;

    NASSERT(x);
    NASSERT(z);
    if(N<=0) return;

	#pragma aligned (x, 16)
	#pragma aligned (z, 16)

    for ( i=0; i<N; i++ )
    {
    	z[i] = (int16_t)abs(x[i]);
    }
}

#endif
