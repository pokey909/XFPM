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
 * Vector Operations vec_cplx2real_multv16x16:
 * vector complex to real multiplicaiton (scaling) operation for 16 bit integer data
 * Each complex number is scaled by a unique number
 * Input x : int16 vector complex (real and imag parts are interleaved)
 * Input y : int16 vector real
 * Output z: int16 vector complex (real and imag parts are interleaved)
 *
 */

#include "common.h"
#define __HiFi3_3z_NDSP__ 1
#define __REFC_OOB__ 0

#if __HiFi3_3z_NDSP__
void vec_cplx2real_multv16x16 (complex16_t* restrict z, complex16_t* restrict x, int16_t* restrict y, int N)
{
    const ae_int16x4* restrict px;
    const ae_int16* restrict py;
    ae_int16* restrict pz;
    ae_int16x4 xt, yt, y0, y1, zt;
    ae_int32x2 d0, d1;
    int i;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT_ALIGN(z, 16);
    if(N<=0) return;

    px = (const ae_int16x4*)x;
    py = (const ae_int16*)y;
    pz = (ae_int16*) z;

    for ( i=0; i<(N>>1); i++ )
    {
		AE_L16X4_IP(xt,(ae_int16x4*)px,sizeof(ae_int16x4));
		AE_L16_IP(y0,castxcc(ae_int16,py),sizeof(int16_t));
		AE_L16_IP(y1,castxcc(ae_int16,py),sizeof(int16_t));
		yt = AE_SEL16_7632(y0, y0);
		yt = AE_SEL16_7632(yt, y1);
		AE_MUL16X4(d0, d1, xt, yt);
		zt = AE_SAT16X4(d0, d1);
		AE_S16X4_IP(zt,(ae_int16x4*)pz,sizeof(ae_int16x4));
    }
    if(N&1)
    {
    	AE_L16X4_IP(xt,(ae_int16x4*)px,sizeof(ae_int16x4));
    	AE_L16_IP(yt,castxcc(ae_int16,py),sizeof(int16_t));
		AE_MUL16X4(d0, d1, xt, yt);
		zt = AE_SAT16X4(d0, d1);
		y0 = AE_SEL16_6543(zt, zt);
		AE_S16_0_IP(y0,castxcc(ae_int16,pz),sizeof(int16_t));
		y1 = AE_SEL16_7632(zt, zt);
		AE_S16_0_IP(y1,castxcc(ae_int16,pz),sizeof(int16_t));
    }
    return;
}

/*
 * vec_cplx2real_mults16x16:
 * vector complex to real multiplicaiton (scaling) operation for 16 bit integer data
 * Entire complex vector is scaled by a fixed number
 * Input x : int16 vector complex (real and imag parts are interleaved)
 * Input y : int16 scalar real
 * Output z: int16 vector complex (real and imag parts are interleaved)
 *
 */
void vec_cplx2real_mults16x16 (complex16_t* restrict z, const complex16_t* restrict x, const int16_t y, int N)
{
    const ae_int16x4* restrict px;
    ae_int16* restrict pz;
    ae_int16x4 xt, yt=y, y0, y1, zt;
    ae_int32x2 d0, d1;
    int i;

    NASSERT(x);
    NASSERT(z);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(z, 16);
    if(N<=0) return;

    px = (const ae_int16x4*)x;
    pz = (ae_int16*) z;

    for ( i=0; i<(N>>1); i++ )
    {
		AE_L16X4_IP(xt,(ae_int16x4*)px,sizeof(ae_int16x4));
		AE_MUL16X4(d0, d1, xt, yt);
		zt = AE_SAT16X4(d0, d1);
		AE_S16X4_IP(zt,(ae_int16x4*)pz,sizeof(ae_int16x4));
    }
    if(N&1)
    {
    	AE_L16X4_IP(xt,(ae_int16x4*)px,sizeof(ae_int16x4));
		AE_MUL16X4(d0, d1, xt, yt);
		zt = AE_SAT16X4(d0, d1);
		y0 = AE_SEL16_6543(zt, zt);
		AE_S16_0_IP(y0,castxcc(ae_int16,pz),sizeof(int16_t));
		y1 = AE_SEL16_7632(zt, zt);
		AE_S16_0_IP(y1,castxcc(ae_int16,pz),sizeof(int16_t));
    }
    return;
}

#else 
void vec_cplx2real_multv16x16 (complex16_t* restrict z, complex16_t* restrict x, int16_t* restrict y, int N)
{
    int i,j;
    int16_t *xptr, *yptr, *zptr;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);

    if(N<=0) return;

    xptr = (int16_t*) x;
    yptr = (int16_t*) y;
    zptr = (int16_t*) z;

    #pragma aligned (x, 16)
    #pragma aligned (y, 16)
    #pragma aligned (z, 16)

    for(i=0,j=0; i<(N*2); i+=2,j++)
    {
        zptr[i] = xptr[i] * yptr[j];
        zptr[i+1] = xptr[i+1] * yptr[j];
    }
}

void vec_cplx2real_mults16x16 (complex16_t* restrict z, complex16_t* restrict x, int16_t y, int N)
{
    int i;
    int16_t *xptr, *zptr;

    NASSERT(x);
    NASSERT(z);

    if(N<=0) return;

    xptr = (int16_t*) x;
    zptr = (int16_t*) z;

    #pragma aligned (x, 16)
    #pragma aligned (z, 16)

    for(i=0; i<(N*2); i+=2)
    {
    	zptr[i] = xptr[i] * y;
    	zptr[i+1] = xptr[i+1] * y;
    }
}
#endif

