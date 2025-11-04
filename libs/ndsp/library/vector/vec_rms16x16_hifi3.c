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
/* Calculates RMS of Input vector
 * inputs: x int16 (vector)
 * output: rmsF  int16 scalar value
 * restrictions: input should be 8 byte aligned
 */
#include "common.h"
#include "NatureDSP_Signal_math.h"

#define __HiFi3_3z_NDSP__ 1

#if __HiFi3_3z_NDSP__
int16_t vec_rms16x16( const int16_t* restrict x, int N )
{
# if XCHAL_HAVE_HIFI3Z
	const ae_int16x4* restrict px;
	ae_int16x4 xt;
	ae_int64 tmp=0;
	int i, rt;
	long long tmpL;
	int16_t rmsF;

	NASSERT(x);
	if(N<=0) return 0;
	px=(const ae_int16x4 *)x;

	for ( i=0; i<(N>>2); i++ )
	{
		AE_L16X4_IP(xt,(ae_int16x4*)px,sizeof(ae_int16x4));
		AE_MULAAAAQ16(tmp,xt,xt);
	}
	for ( i=0; i<(N&3); i++ )
	{
		ae_int16x4 d;
		AE_L16_IP(d,castxcc(ae_int16,px),sizeof(int16_t));
		AE_MULAAAAQ16(tmp,d,d);
	}

	tmpL = tmp;
	tmpL = (tmpL/ N);
	tmpL = tmpL << 1; 				  	// To convert from q62 to q63
	rt = scl_sqrt32x32(tmpL);
	rmsF = (int16_t) (rt >> 16);

#else //HiFi3
	const ae_int16x4* restrict px;
		ae_int16x4 xt;
		ae_int16x4 sumx4=0;
		ae_int16x4 sum=0;
		int i,sqrt32 ;
		long long sum64Bit;
		int16_t rmsF;

		NASSERT(x);
		if(N<=0) return 0;
		px=(const ae_int16x4 *)x;

		for ( i=0; i<(N>>2); i++ )
		{
			AE_L16X4_IP(xt,(ae_int16x4*)px,sizeof(ae_int16x4));
			AE_MULAAR16P16X4S_vector(sumx4,xt,xt);
		}
		for ( i=0; i<(N&3); i++ )
		{
			ae_int16x4 d;
			AE_L16_IP(d,castxcc(ae_int16,px),sizeof(int16_t));
			AE_MULAAR16P16X4S_vector(sum,d,d);
		}
		sumx4 = AE_INT16X4_RADD(sumx4);
		sumx4 = AE_ADD16S_vector(sumx4,sum);
		sum64Bit = AE_MOVINT64_FROMINT16X4(sumx4);
		sum64Bit = AE_SRAI64(sum64Bit,48); //Discard lower 48 bits
		sum64Bit = (sum64Bit/ N);
		sum64Bit = sum64Bit << 1; 				  	// To convert from q62 to q63
		sqrt32 = scl_sqrt32x32(sum64Bit);
		rmsF = (int16_t) (sqrt32 >> 16);
#endif

	return  rmsF;
}

#elif __REFC_OOB__
#include "math.h"
int16_t vec_rms16x16( const int16_t* restrict x, int N )
{
	int i;
	long long sumt=0;
	int32_t rms;
	int16_t rmsF;

	NASSERT(x);
	NASSERT_ALIGN(x, 16);

	if(N<=0) return 0;
	#pragma aligned (x, 16)

	for (i=0; i<N; i++)
	{
		sumt += (long long) (x[i] * x[i]);			// q15 * q15 -> q30
	}
	rms = (int32_t) (sumt / N);						//
	rms = rms << 1; 								// from q30 to q31
	rms = scl_sqrt32x32(rms);						//
	rmsF = (int16_t) (rms >> 16);

	return rmsF;
}

#elif __HiFi3_3z_NDSP__

#elif __HiFi5NDSP__

#endif
