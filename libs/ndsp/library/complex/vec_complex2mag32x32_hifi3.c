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
/* complex to complex magnitude calculation of 32 bit data
 * input  x : complex int32
 * output z : real int32
 */
#include "common.h"
#include "NatureDSP_Signal_math.h"
#define __HiFi3_3z_NDSP__ 1

#if __HiFi3_3z_NDSP__
void vec_complex2mag32x32 (int32_t* restrict z, complex32_t* restrict x, int N)
{
    const ae_int32x2* restrict pX;
          ae_int32* restrict pZ;
          ae_int32x2 x0;
          ae_int64  sumt0=0;
    int i, rt0;

    NASSERT(x);
    NASSERT(z);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(z, 16);
    if(N<=0) return;

    pX = (const ae_int32x2*)x;
    pZ = (ae_int32*) z;

    for ( i=0; i<N; i++ )
    {
        AE_L32X2_IP(x0,pX,sizeof(ae_int32x2));

		#if XCHAL_HAVE_HIFI3Z
        	sumt0 = AE_MULZAAFD32S_HH_LL(x0, x0);
		#else
        sumt0 = AE_MULF32S_HH(x0, x0);
		AE_MULAF32S_LL(sumt0, x0, x0);
		#endif
		rt0 = scl_sqrt64x32(sumt0);
        ae_int32 out = AE_MOVDA32(rt0);
        AE_S32_L_IP(out, pZ, sizeof(ae_int32));
    }
    //TBD for OPT:
    /*  Create a scratch buffer and store all the squared values.
     *  Call vector_sqrt instead of scl_sqrt
     */
    return;
}

#elif __REFC_OOB__
#include "math.h"
void vec_complex2mag32x32 (int32_t* restrict z, complex32_t* restrict x, int N)
{
    int i,j=0;
    int32_t *xptr, *zptr, res;
    long long sqsumt;

    NASSERT(x);
    NASSERT(z);

    if(N<=0) return;

    xptr = (int32_t*) x;
    zptr = (int32_t*) z;

    #pragma aligned (x, 16)
    #pragma aligned (z, 16)

    for ( i=0; i<(N*2); i+=2 )
    {
        sqsumt = 0;
        sqsumt = (long long) xptr[i] * xptr[i];
        sqsumt += (long long) xptr[i+1] * xptr[i+1];
        res = (int32_t) sqrt(sqsumt);
        zptr[j++] = res;
    }
}
#endif
