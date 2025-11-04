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
/* complex to complex magnitude calculation of 16 bit data
 * input  x : complex int16
 * output z : real int16
 */
#include "common.h"
#include "polyrsqrtq23_tbl.h"
int16_t scl_sqrt32x16(int32_t x);
#define __HiFi3_3z_NDSP__ 1

#if __HiFi3_3z_NDSP__
void vec_complex2mag16x16 (int16_t* restrict z, complex16_t* restrict x, int N)
{
    const ae_int16x4* restrict pX;
          ae_int32* restrict pZ;
          ae_int16x4 x0;
          ae_int32  sumt0, sumt1;
    int16_t rt0, rt1;
    int i;

    NASSERT(x);
    NASSERT(z);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(z, 16);
    if(N<=0) return;

    pX = (const ae_int16x4*)x;
    pZ = (ae_int32*) z;

    for ( i=0; i<(N>>2); i++ )
    {
        AE_L16X4_IP(x0,pX,sizeof(ae_int16x4));
        sumt0 = AE_MULZAAFD16SS_33_22(x0, x0);
        rt0 = scl_sqrt32x16(sumt0);
        sumt1 = AE_MULZAAFD16SS_11_00(x0, x0);
        rt1 = scl_sqrt32x16(sumt1);

        ae_int16x4 out = AE_MOVDA16X2(rt1, rt0);
        ae_int32x2 outTmp = AE_MOVINT32X2_FROMINT16X4(out);
        AE_S32_L_IP(outTmp, pZ, sizeof(ae_int32));

        AE_L16X4_IP(x0,pX,sizeof(ae_int16x4));
        sumt0 = AE_MULZAAFD16SS_33_22(x0, x0);
        rt0 = scl_sqrt32x16(sumt0);
        sumt1 = AE_MULZAAFD16SS_11_00(x0, x0);
        rt1 = scl_sqrt32x16(sumt1);

        out = AE_MOVDA16X2(rt1, rt0);
        outTmp = AE_MOVINT32X2_FROMINT16X4(out);
        AE_S32_L_IP(outTmp, pZ, sizeof(ae_int32));

    }
    for ( i=0; i<(N&3); i++ )
    {
        ae_int16x4 x1,zt;
        AE_L16_IP(x0,castxcc(ae_int16,pX),sizeof(int16_t));
        AE_L16_IP(x1,castxcc(ae_int16,pX),sizeof(int16_t));
        x0 = AE_SEL16_7362(x0, x1);
        sumt0 = AE_MULZAAFD16SS_33_22(x0, x0);
        rt0 = scl_sqrt32x16(sumt0);
        zt = AE_MOVINT16X4_FROMF16(rt0);
        AE_S16_0_IP(zt,castxcc(ae_int16,pZ),sizeof(int16_t));
    }
    return;
}

#elif __REFC_OOB__
#include "math.h"
void vec_complex2mag16x16 (int16_t* restrict z, complex16_t* restrict x, int N)
{
    int i,j=0;
    int16_t *xptr, *zptr, res;
    long sqsumt;

    NASSERT(x);
    NASSERT(z);

    if(N<=0) return;

    xptr = (int16_t*) x;
    zptr = (int16_t*) z;

    #pragma aligned (x, 16)
    #pragma aligned (z, 16)

    for ( i=0; i<(N*2); i+=2 )
    {
        sqsumt = 0;
        sqsumt = (long) xptr[i] * xptr[i];
        sqsumt += (long) xptr[i+1] * xptr[i+1];
        res = (int16_t) sqrt(sqsumt);
        zptr[j++] = res;
    }
}
#endif
int16_t scl_sqrt32x16(int32_t x)
{
    ae_int32x2 X,Y,D,R;
    ae_f32x2 f;
    xtbool2 lezero;
    int sh;
    /* load, take exponent */
    X=(x);
    sh=AE_NSAZ32_L(X)&~1;
    X=AE_SLAA32(X,sh-8);
    lezero=AE_LT32(X,0);
    /* compute rsqrt */
    R=polyrsqrtq23[0];
    f=polyrsqrtq23[1]; AE_MULAFP24X2RA(f,AE_MOVF24X2_FROMINT32X2(X),AE_MOVF24X2_FROMINT32X2(R)); R=f;
    f=polyrsqrtq23[2]; AE_MULAFP24X2RA(f,AE_MOVF24X2_FROMINT32X2(X),AE_MOVF24X2_FROMINT32X2(R)); R=f;
    f=polyrsqrtq23[3]; AE_MULAFP24X2RA(f,AE_MOVF24X2_FROMINT32X2(X),AE_MOVF24X2_FROMINT32X2(R)); R=f;
    f=polyrsqrtq23[4]; AE_MULAFP24X2RA(f,AE_MOVF24X2_FROMINT32X2(X),AE_MOVF24X2_FROMINT32X2(R)); R=f;
    /* reiterate rsqrt */
    R=AE_SLAI24S(AE_MOVF24X2_FROMINT32X2(R),3);
    D=AE_MULFP24X2RA(AE_MOVF24X2_FROMINT32X2(R),AE_MOVF24X2_FROMINT32X2(R));
    f=0x80000; AE_MULSFP24X2RA(f,AE_MOVF24X2_FROMINT32X2(D),AE_MOVF24X2_FROMINT32X2(X)); D=f;
    D=AE_MULFP24X2RA(AE_MOVF24X2_FROMINT32X2(D),AE_MOVF24X2_FROMINT32X2(R));
    D=AE_SLAI32(D,3);
    R=AE_ADD24S(R,D);
    /* compute sqrt and rescale back */
    Y=AE_MULFP24X2RA(AE_MOVF24X2_FROMINT32X2(X),AE_MOVF24X2_FROMINT32X2(R));
    X=AE_SLAA32S(Y,14-(sh>>1)-4);
    AE_MOVT32X2(X,0x80000000,lezero);
    return (int16_t)(AE_MOVAD32_L(X)>>16);
}
