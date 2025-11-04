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
/* complex to complex multiplication of 32 bit data
 * input  x,y : complex int32
 * output z   : complex int32
 */
#include "common.h"

#define __HiFi3_3z_NDSP__ 1
#define __REFC_OOB__ 0  //ISA oob both OK, but ISA not on hf3

#if __HiFi3_3z_NDSP__
void vec_cplx2cplx_mult32x32 (complex32_t* restrict z, complex32_t* restrict x, complex32_t* restrict y, int N)
{
#if 1
const ae_int32x2* restrict pX;
const ae_int32x2* restrict pY;
      ae_int32x2* restrict pZ;
      ae_int32x2 x0,x1,y0,y1;
      ae_int32x2 z0, z1;
      ae_int32x2 x0x1RR, x0x1II, y0y1RR, y0y1II;
      ae_int32x2 zR, zI;
    int i;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(y, 8);
    NASSERT_ALIGN(z, 8);
    if(N<=0) return;

    pX = (const ae_int32x2*)x;
    pY = (const ae_int32x2*)y;
    pZ = (ae_int32x2*) z;

    for ( i=0; i<N>>1; i++ )
    {
        AE_L32X2_IP(x0,pX,sizeof(ae_int32x2));
        AE_L32X2_IP(x1,pX,sizeof(ae_int32x2));
        AE_L32X2_IP(y0,pY,sizeof(ae_int32x2));
        AE_L32X2_IP(y1,pY,sizeof(ae_int32x2));

        x0x1RR = AE_SEL32_HH(x0, x1);
        x0x1II = AE_SEL32_LL(x0, x1);

        y0y1RR = AE_SEL32_HH(y0, y1);
        y0y1II = AE_SEL32_LL(y0, y1);

        zR = AE_MULP32X2(x0x1RR, y0y1RR);    //R*R
        AE_MULSP32X2(zR, x0x1II, y0y1II); //xRyR-xIyI

        zI = AE_MULP32X2(x0x1RR, y0y1II);    //R*I
        AE_MULAP32X2(zI, x0x1II, y0y1RR); //xRyI-xIyR

        z0 = AE_SEL32_HH(zR, zI);
        z1 = AE_SEL32_LL(zR, zI);
        AE_S32X2_IP(z0,pZ,sizeof(ae_int32x2));
        AE_S32X2_IP(z1,pZ,sizeof(ae_int32x2));
    }
    //Add tail handling
    if(N&1) //last 1 complex sample
    {
    	ae_int32x2 y0IR,zRRII,zRIIR,z0R,z0I;
    	 AE_L32X2_IP(x0,pX,sizeof(ae_int32x2));
    	 AE_L32X2_IP(y0,pY,sizeof(ae_int32x2));
    	 y0IR = AE_SEL32_LH(y0, y0);

    	 zRRII = AE_MULP32X2(x0, y0);    //R*R, I*I
    	 zRIIR = AE_MULP32X2(x0, y0IR);   //R*I, I*R

#if XCHAL_HAVE_HIFI3Z  //check if we use saturating ISA in main and tail loop
    	 z0R = AE_SUBADD32_HL_LH(zRRII, zRRII); //neglect lower
    	 z0I = AE_SUBADD32_HL_LH(zRIIR, zRIIR); //neglect Upper
#else
    	 ae_int32x2 zIIRR,zIRRI;
    	 zIIRR = AE_SEL32_LH(zRRII, zRRII);
    	 zIRRI = AE_SEL32_LH(zRIIR, zRIIR);
		 z0R = AE_SUBADD32(zRRII, zIIRR); //neglect lower
    	 z0I = AE_SUBADD32(zRIIR, zIRRI); //neglect Upper

#endif
    	 z0 = AE_SEL32_HL(z0R, z0I);
    	 AE_S32X2_IP(z0,pZ,sizeof(ae_int32x2));
    }

#else  //non optimized way
    const ae_int32x2* restrict pX;
    const ae_int32x2* restrict pY;
          ae_int32x2* restrict pZ;
    ae_int32x2 x0, y0,y0IR, z0;
    ae_int32x2 z0Real, z0Img, z0_H, z0_L;
    int i;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(y, 8);
    NASSERT_ALIGN(z, 8);
    if(N<=0) return;

    pX = (const ae_int32x2*)x;
    pY = (const ae_int32x2*)y;
    pZ = (ae_int32x2*) z;

    for ( i=0; i<N; i++ )
       {
		   AE_L32X2_IP(x0,pX,sizeof(ae_int32x2));
		   AE_L32X2_IP(y0,pY,sizeof(ae_int32x2));
           y0IR = AE_SEL32_LH(y0,y0);

           z0Real = AE_MULP32X2(x0,y0); //xr*yr , xi*yi
           z0Img = AE_MULP32X2(x0,y0IR); //xr*yi , xi*yr

           z0_H = AE_SEL32_HH(z0Real,z0Img);
           z0_L = AE_SEL32_LL(z0Real,z0Img);
           z0 = AE_SUBADD32(z0_H,z0_L);
           AE_S32X2_IP(z0,pZ,sizeof(ae_int32x2));
       }
#endif

    return;
}

#elif __REFC_OOB__
void vec_cplx2cplx_mult32x32 (complex32_t* restrict z, complex32_t* restrict x, complex32_t* restrict y, int N)
{
    int i;
    int32_t *xptr, *yptr, *zptr;

    xptr = (int32_t*) x;
    yptr = (int32_t*) y;
    zptr = (int32_t*) z;

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

