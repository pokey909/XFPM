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
  NatureDSP Signal Processing Library. FFT part
    C code optimized for HiFi3
  Integrit, 2006-2017
*/
/*===========================================================================

    fft_unpack24to32_ie:
    Unpack 24 bits packed data to 32 bits unpacked.

    fft_pack32to24_ie:
    Pack 32-bits data into 24-bits packed. 

    N - number of input/output words.

    Called from: 
    fft_real_32x16_ie_24p
    fft_real_24x24_ie_24p
    ifft_real_32x16_ie_24p
    ifft_real_24x24_ie_24p

===========================================================================*/
#include "NatureDSP_Signal_fft.h"
#include "common.h"

void fft_pack32to24_ie(uint8_t *x,  uint8_t *y, int N)
{
    ae_f24x2 d; 
    ae_valign v = AE_ZALIGN64();
    const ae_int32x2 *restrict pin  = (const ae_int32x2 *)x; 
    void *restrict pout = (ae_f24x2 *)y; 
    int i; 

    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT((N&1)==0);

    for(i=0; i<N/2; i++)
    {
        ae_int32x2 t;
        AE_L32X2_IP(t, pin, sizeof(*pin) );
        /* Do rounding before packing */
        t = AE_ADD32(t, 0x80);
        t = AE_SRAI32(t, 8);
        d = AE_MOVF24X2_FROMINT32X2(t);
        AE_SA24X2_IP(d, v, pout); 
    }
    AE_SA64POS_FP(v, pout); 
} /* fft_pack32to24_ie() */

void fft_unpack24to32_ie(uint8_t *x,  uint8_t *y, int N)
{
    ae_int24x2 d; 
    ae_valign v = AE_LA64_PP(x);
    void *pin  = x; 
    ae_f24x2 *pout = (ae_f24x2 *)y; 
    int i; 

    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT((N&1)==0);

    for(i=0; i<N/2; i++)
    {
        AE_LA24X2_IP(d, v, pin);
        AE_S32X2F24_IP(d, pout, sizeof(*pout)); 
    }
} /* fft_unpack24to32_ie() */
