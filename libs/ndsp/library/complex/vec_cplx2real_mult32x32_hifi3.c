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
 * Vector Operations vec_cplx2real_multv32x32:
 * vector complex to real multiplicaiton (scaling) operation for 32 bit integer data
 * Each complex number is scaled by a unique number
 * Input x : int32 vector complex (real and imag parts are interleaved)
 * Input y : int32 vector real
 * Output z: int32 vector complex (real and imag parts are interleaved)
 *
 */
#include "common.h"

#define __HiFi3_3z_NDSP__ 1
#define __REFC_OOB__ 0   //OOB is better

#if __HiFi3_3z_NDSP__
void vec_cplx2real_multv32x32 (complex32_t* restrict z, complex32_t* restrict x, int32_t* restrict y, int N)
{
    const ae_int32x2* restrict pX;
    const ae_int32* restrict pY;
    ae_int32x2* restrict pZ;
    ae_int32x2 x0;
    ae_int32x2 y00;
    ae_int32x2 z0,zt;
    ae_int64 z0Real, z0Imag;
    int i;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(y, 8);
    NASSERT_ALIGN(z, 8);
    if(N<=0) return;

    pX = (const ae_int32x2*)x;
    pY = (const ae_int32*)y;
    pZ = (ae_int32x2*) z;

    for ( i=0; i<N; i++ )
    {
        AE_L32X2_IP(x0,pX,sizeof(ae_int32x2));
        AE_L32_IP(y00,pY,sizeof(ae_int32));
        z0Real = AE_MUL32_HH(x0, y00);
        z0Imag = AE_MUL32_LL(x0, y00);
        z0 = AE_MOVINT32X2_FROMINT64(z0Real);
        zt = AE_MOVINT32X2_FROMINT64(z0Imag);
        z0 = AE_SEL32_LL(z0, zt);
        AE_S32X2_IP(z0,pZ,sizeof(ae_int32x2));
    }
    return;
}
/*
 * vec_cplx2real_mults32x32:
 * vector complex to real multiplicaiton (scaling) operation for 32 bit integer data
 * Entire complex vector is scaled by a fixed number
 * Input x : int32 vector complex (real and imag parts are interleaved)
 * Input y : int32 scalar real
 * Output z: int32 vector complex (real and imag parts are interleaved)
 *
 */
void vec_cplx2real_mults32x32 (complex32_t* restrict z, complex32_t* restrict x, int32_t y, int N)
{
    const ae_int32x2* restrict pX;
    ae_int32x2* restrict pZ;
    ae_int32x2 x0;
    ae_int32x2 y00=y;
    ae_int32x2 z0,zt;
    ae_int64 z0Real, z0Imag;
    int i;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(y, 8);
    NASSERT_ALIGN(z, 8);
    if(N<=0) return;

    pX = (const ae_int32x2*)x;
    pZ = (ae_int32x2*) z;

    for ( i=0; i<N; i++ )
    {
        AE_L32X2_IP(x0,pX,sizeof(ae_int32x2));
        z0Real = AE_MUL32_HH(x0, y00);
        z0Imag = AE_MUL32_LL(x0, y00);
        z0 = AE_MOVINT32X2_FROMINT64(z0Real);
        zt = AE_MOVINT32X2_FROMINT64(z0Imag);
        z0 = AE_SEL32_LL(z0, zt);
        AE_S32X2_IP(z0,pZ,sizeof(ae_int32x2));
    }
    return;
}

#elif __REFC_OOB__
void vec_cplx2real_multv32x32 (complex32_t* restrict z, complex32_t* restrict x, int32_t* restrict y, int N)
{
    int i,j;
    int32_t *xptr, *yptr, *zptr;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);

    if(N<=0) return;

    xptr = (int32_t*) x;
    yptr = (int32_t*) y;
    zptr = (int32_t*) z;

    #pragma aligned (x, 16)
    #pragma aligned (y, 16)
    #pragma aligned (z, 16)

    for(i=0,j=0; i<(N*2); i+=2,j++)
    {
        zptr[i] = xptr[i] * yptr[j];
        zptr[i+1] = xptr[i+1] * yptr[j];
    }
}

void vec_cplx2real_mults32x32 (complex32_t* restrict z, complex32_t* restrict x, int32_t y, int N)
{
    int i;
    int32_t *xptr, *zptr;

    NASSERT(x);
    NASSERT(z);

    if(N<=0) return;

    xptr = (int32_t*) x;
    zptr = (int32_t*) z;

    #pragma aligned (x, 16)
    #pragma aligned (z, 16)

    for(i=0; i<(N*2); i+=2)
    {
    	zptr[i] = xptr[i] * y;
    	zptr[i+1] = xptr[i+1] * y;
    }
}
#endif
