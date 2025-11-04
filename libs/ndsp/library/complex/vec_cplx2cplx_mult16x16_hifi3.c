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
/*  complex to complex multiplication of 16 bit integers
 * Input x,y: int16 vector complex (real and imag parts are interleaved)
 * Output z : int16 vector complex (real and imag parts are interleaved)
 */
#include "common.h"

#if(XCHAL_HAVE_HIFI3Z)
void vec_cplx2cplx_mult16x16 (complex16_t* restrict z, complex16_t* restrict x, complex16_t* restrict y, int N)
{
const ae_int16x4* restrict pX;
const ae_int16x4* restrict pY;
      ae_int16x4* restrict pZ;
      ae_int16x4 x0,y0, y0ir;
      ae_int16x4 z0;
      ae_int16x4 z0Real,z0Img,z0Imgtmp,z0Imgneg, z0_H, z0_L;
    int i;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(y, 8);
    NASSERT_ALIGN(z, 8);
    if(N<=0) return;

    pX = (const ae_int16x4*)x;
    pY = (const ae_int16x4*)y;
    pZ = (ae_int16x4*) z;

    for ( i=0; i<(N>>1); i++ )
    {
        AE_L16X4_IP(x0,castxcc(ae_int16x4,pX),sizeof(ae_int16x4));
        AE_L16X4_IP(y0,castxcc(ae_int16x4,pY),sizeof(ae_int16x4));
        y0ir = AE_SEL16_2301(y0,y0);

        z0Real = AE_MUL16X4_vector(x0,y0); //xr*yr , xi*yi
        z0Img = AE_MUL16X4_vector(x0,y0ir); //xr*yi , xi*yr

        z0_H = AE_SEL16_7362(z0Real,z0Img);
        z0_L = AE_SEL16_5140(z0Real,z0Img);

        z0Real = AE_SEL16_7632(z0_H,z0_L);
        z0Img  = AE_SEL16_5410(z0_H,z0_L);

        z0Imgneg = AE_INT16X4_NEG(z0Img);
        z0Imgtmp  = AE_SEL16_7520(z0Imgneg,z0Img);
        z0Img     = AE_SEL16_7520(z0Imgtmp,z0Imgtmp);

        z0 = AE_ADD16 (z0Real,z0Img);

        AE_S16X4_IP(z0,castxcc(ae_int16x4,pZ),sizeof(ae_int16x4));
    }
    if(N&1)
    {
        AE_L16X4_IP(x0,castxcc(ae_int16x4,pX),sizeof(ae_int16x4));
        AE_L16X4_IP(y0,castxcc(ae_int16x4,pY),sizeof(ae_int16x4));
        y0ir = AE_SEL16_2301(y0,y0);

        z0Real = AE_MUL16X4_vector(x0,y0); //xr*yr , xi*yi
        z0Img = AE_MUL16X4_vector(x0,y0ir); //xr*yi , xi*yr

        z0_H = AE_SEL16_7362(z0Real,z0Img);
        z0_L = AE_SEL16_5140(z0Real,z0Img);

        z0Real = AE_SEL16_7632(z0_H,z0_L);
        z0Img  = AE_SEL16_5410(z0_H,z0_L);

        z0Imgneg = AE_INT16X4_NEG(z0Img);
        z0Imgtmp  = AE_SEL16_7520(z0Imgneg,z0Img);
        z0Img     = AE_SEL16_7520(z0Imgtmp,z0Imgtmp);

        z0 = AE_ADD16 (z0Real,z0Img);

		y0 = AE_SEL16_6543(z0, z0);
		AE_S16_0_IP(y0,castxcc(ae_int16,pZ),sizeof(int16_t));
		y0 = AE_SEL16_7632(z0, z0);
		AE_S16_0_IP(y0,castxcc(ae_int16,pZ),sizeof(int16_t));

    }
    return;
}
#else
void vec_cplx2cplx_mult16x16 (complex_fract16* restrict z, complex_fract16* restrict x, complex_fract16* restrict y, int N)
{
    int i;
    int16_t *xptr, *yptr, *zptr;

    xptr = (int16_t*) x;
    yptr = (int16_t*) y;
    zptr = (int16_t*) z;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);

    if(N<=0) return;

    #pragma aligned (xptr, 16)
    #pragma aligned (yptr, 16)
    #pragma aligned (zptr, 16)

    for ( i=0; i<(N*2); i+=2 )
    {
        zptr[i]     = (xptr[i] * yptr[i])   - (xptr[i+1] * yptr[i+1]);
        zptr[i+1]   = (xptr[i] * yptr[i+1]) + (xptr[i+1] * yptr[i]);
    }
}
#endif

