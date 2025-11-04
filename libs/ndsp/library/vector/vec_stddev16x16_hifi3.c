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
/*  Standard Deviation of a vector
	Input 	: x - int16 vector
	Output 	: res - int16 scalar
*/
#include "common.h"
#include "NatureDSP_Signal_math.h"
#define __HiFi3_3z_NDSP__ 1


#if __HiFi3_3z_NDSP__
int16_t vec_stddev16x16 (const int16_t* restrict x, int N)
{
	const ae_int16x4* restrict px;
	ae_int16x4 xt, zt=0, meanx4;
	ae_int16   zf=0;
	int16_t res;
	ae_int16x4 sumx4=0, sum=0;
	long long sum64Bit;
	int32_t sqrt32;

	int i;
	NASSERT(x);
	NASSERT_ALIGN(x, 16);
	if(N<2) return 0;

	// Calculate Mean
	px=(const ae_int16x4 *)x;
	for ( i=0; i<(N>>2); i++ )
	{
		AE_L16X4_IP(xt,(ae_int16x4*)px,sizeof(ae_int16x4));
		zt = AE_ADD16S(zt,xt);
	}
	zf = AE_INT16X4_RADD(zt);
	for ( i=0; i<(N&3); i++ )
	{
		ae_int16x4 x0;
		AE_L16_IP(x0,castxcc(ae_int16,px),sizeof(int16_t));
		zf = AE_ADD16S((ae_int16x4)zf,(ae_int16x4)x0);
	}
	res = zf;
	res = (int16_t) (res/N);
	meanx4 = AE_MOVINT16X4_FROMF16(res);

	/*  Compute  variance */
	px=(const ae_int16x4 *)x;
	for ( i=0; i<(N>>2); i++ )
	{
		AE_L16X4_IP(xt,(ae_int16x4*)px,sizeof(ae_int16x4));
		zt = AE_SUB16S_vector(xt, meanx4);
		AE_MULAAR16P16X4S_vector(sumx4,zt,zt);
	}

	for ( i=0; i<(N&3); i++ )
	{
		ae_int16x4 x0;
		AE_L16_IP(x0,castxcc(ae_int16,px),sizeof(int16_t));
		x0 = AE_SUB16S_vector((ae_int16x4)x0, (ae_int16x4)meanx4);
		AE_MULAAR16P16X4S_vector(sum,x0,x0);
	}

	sumx4 = AE_INT16X4_RADD(sumx4);
	sumx4 = AE_ADD16S_vector(sumx4,sum);
	sum64Bit = AE_MOVINT64_FROMINT16X4(sumx4);
	sum64Bit = AE_SRAI64(sum64Bit,48); //Discard lower 48 bits
	sum64Bit = (sum64Bit/ (N-1));
	// Calculate standard deviation
	sum64Bit = sum64Bit << 1; 				  	// To convert from q62 to q63
	sqrt32 = scl_sqrt32x32(sum64Bit);
	res = (int16_t) (sqrt32 >> 16);
	return res;
}


#elif __REFC_OOB__
#include "math.h"
int16_t vec_stddev16x16 (const int16_t* restrict x, int N)
{
	int i;
	int16_t sum, stddev, mean, var;
	long long sumt = 0;

	NASSERT(x);
	NASSERT_ALIGN(x, 16);

	if(N<=0) return 0;
	#pragma aligned (x, 16)

	/*  Compute the Mean of all elements */
	for (i=0; i<N; i++)
	{
		sumt += x[i];
	}
	mean = (int16_t) (sumt / N);

	/*  Compute  variance */
	sumt = 0;
	for (i=0; i<N; i++)
	{
		sum = x[i] - mean;
		sumt += (long long) (sum*sum);
	}
	var = (int16_t) sumt / (N-1);
	stddev = (int16_t) sqrt(var);

	return stddev;
}
#endif
