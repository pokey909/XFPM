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
    Discrete Cosine Transform, Type II
    C code optimized for HiFi3
   Integrit, 2006-2017
*/

/* Signal Processing Library API. */
#include "NatureDSP_Signal_fft.h"
#include "common.h"
#include "dct2_twd.h"

/*
    Reference Matlab code for Radix4 DCT-II

    y(1:N/2)=x(1:2:N);
    y(N:-1:N/2+1)=x(2:2:N);
    y=y(1:2:N)+i*y(2:2:N);
    Y=fft(y);
    cosi=exp(2*i*pi*(0:N/4-1)/(N));
    w=exp(-i*pi/2*(0:N-1)/N);
    % DCT split algorithm
    Y0=Y(1);
    T0=real(Y0)+imag(Y0);
    T1=real(Y0)-imag(Y0);
    z(1      )= real(T0)*sqrt(2)/2;
    z(N/2+1  )= real(T1)*sqrt(2)/2;
    for k=2:N/4
        Y0=Y(k);
        Y1=Y(N/2+2-k);
        COSI=cosi(k);
        W1=w(k);
        W2=w(N/2+2-k);
        S=Y0+Y1;
        D=Y0-Y1;
        T0=i*real(D)+imag(S);
        T1=i*imag(D)+real(S);
        T1=T1/2;
        Y0= (imag(T0)*imag(COSI)-real(T0)*real(COSI)) + ...
           i*(real(T0)*imag(COSI)+imag(T0)*real(COSI));
        Y0=Y0/2;
        T0=T1-Y0;
        T1=T1+Y0;
        z(k      )= real(T0)*real(W1)-imag(T0)*imag(W1);
        z(N+2-k  )=-real(T0)*imag(W1)-imag(T0)*real(W1);
        z(N/2+2-k)= real(T1)*real(W2)+imag(T1)*imag(W2);
        z(N/2+k  )=-real(T1)*imag(W2)+imag(T1)*real(W2);
    end
    W1=w(N/4+1);
    T0=Y(N/4+1);
    z(N/4+1  )= real(T0)*real(W1)+imag(T0)*imag(W1);
    z(N+1-N/4)=-real(T0)*imag(W1)+imag(T0)*real(W1);
    z=z*sqrt(2/N);

*/

/*
   scaled fft with reordering
   NOTE: y is input and output, x - temporary
*/
void fft16_24x24(int32_t *y, int32_t *x, const int32_t *ptwd);/* N=16 */
void fft32_24x24(int32_t *y, int32_t *x, const int32_t *ptwd);/* N=32 */
/* pointer to complex fft function with reordering */
typedef void(*cfftProc_func)(int32_t *y, int32_t *x, const int32_t *ptwd);

/*-------------------------------------------------------------------------
  Discrete Cosine Transform.
  These functions apply DCT (Type II, Type IV) to input.
  Scaling:
      +-----------------------+--------------------------------------+
      |      Function         |           Scaling options            |
      +-----------------------+--------------------------------------+
      |       dct_16x16       |  3 - fixed scaling before each stage |
      |       dct_24x24       |  3 - fixed scaling before each stage |
      |       dct_32x16       |  3 - fixed scaling before each stage |
      |       dct_32x32       |  3 - fixed scaling before each stage |
      |       dct4_24x24      |  3 - fixed scaling before each stage |
      |       dct4_32x16      |  3 - fixed scaling before each stage |
      |       dct4_32x32      |  3 - fixed scaling before each stage |
      +-----------------------+--------------------------------------+
  NOTES:
     1. DCT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED after 
     the call.
     2. N - DCT size (depends on selected DCT handle)

  Precision: 
  16x16  16-bit input/outputs, 16-bit twiddles
  24x24  24-bit input/outputs, 24-bit twiddles
  32x16  32-bit input/outputs, 16-bit twiddles
  32x32  32-bit input/outputs, 32-bit twiddles
  f      floating point

  Input:
  x[N]        input signal
  h           DCT handle
  scalingOpt  scaling option (see table above) 
              not applicable to the floating point function
  Output:
  y[N]        transform output
  
  Returned value:
              total number of right shifts occurred during scaling 
              procedure 
  Restriction:
  x,y         should not overlap
  x,y         aligned on 8-bytes boundary
-------------------------------------------------------------------------*/
int dct_24x24(int32_t* y,
              int32_t* x,
              dct_handle_t h,
              int scalingOpt)
{
    int i, N, N4, cfftIx;
    const cfftProc_func cfftTbl[] =
    {
        fft16_24x24,
        fft32_24x24
    };
    const tdct2_twd *descr=(const tdct2_twd *)h;
    const ae_f24x2 * restrict p_y0, * restrict p_y1,
                   * restrict p_cos, * restrict p_w;
          ae_f24x2 * restrict p_z0, * restrict p_z1;
    ae_f24x2    vF0, vF1, vC0, vW0, vW1;
    ae_f32x2    vFR;
    ae_int32x2  vA0, vA1, vB0, vB1, vR0;

    const int step_back = -(int)sizeof(ae_f24x2);

    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT(scalingOpt==3);
    NASSERT(descr->magic==MAGIC_DCT2_32);
    N = descr->N;
    NASSERT(N==32 || N==64);

    /* fft of half-size with reordering */
    cfftIx = 25-NSA(N);
    cfftTbl[cfftIx](x, y, (const int32_t *)descr->fft_twd);

    /* split algorithm */
    N4 = N >> 2;
    p_cos = (const ae_f24x2 *)(descr->rfft_split_twd)+1;
    p_w   = (const ae_f24x2 *)(descr->dct_twd);
    p_y0  = (const ae_f24x2 *)(x);
    p_y1  = (const ae_f24x2 *)(x+N-2);
    p_z0  = (      ae_f24x2 *)(y);
    p_z1  = (      ae_f24x2 *)(y+N-2);

    /* load data and prepare pointers for pre-increment
       first and last samples */
    AE_L32X2F24_IP(vF0, p_y0, 8);
    AE_L32X2F24_IP(vW0, p_w, 8);
    vA0 = (vF0);
    vA1 = AE_SEL32_LH(vA0, vA0);
    vB0 = AE_ADDSUB32S(vA1, vA0);
  
    vF0 = AE_MOVF24X2_FROMINT32X2(vB0);
    vFR = AE_MULFP24X2RA(vF0, vW0);
    vR0 = (vFR);

    AE_L32X2F24_IP(vF0, p_y0, 8);
    AE_L32X2F24_XP(vF1, p_y1, step_back);
    AE_L32X2F24_IP(vC0, p_cos, 8);
    AE_L32X2F24_IP(vW0, p_w, 8);
    AE_L32X2F24_IP(vW1, p_w, 8);

    vB0 = AE_SUBADD32S(vF0, vF1);
    vB1 = AE_ADDSUB32S(vF0, vF1);

    vF0 = AE_MOVF24X2_FROMINT32X2(vB0);
    vFR = AE_MULFC24RA(vF0, vC0);

    vA0 = AE_SRAI32(vFR, 1);
    vA1 = AE_SRAI32(vB1, 1);

    vB0 = AE_SUB32S(vA1, vA0);
    vB1 = AE_ADD32S(vA1, vA0);

    vF0 = AE_MOVF24X2_FROMINT32X2(vB0);
    vF1 = AE_MOVF24X2_FROMINT32X2(vB1);

    /*
      z(k      )= real(T0)*real(W1)-imag(T0)*imag(W1);
      z(N+2-k  )=-real(T0)*imag(W1)-imag(T0)*real(W1);
      z(N/2+2-k)= real(T1)*real(W2)+imag(T1)*imag(W2);
      z(N/2+k  )=-real(T1)*imag(W2)+imag(T1)*real(W2);
    */
    vA0 = AE_MULFC24RA(vF0, vW0);
    vA1 = AE_NEG24S(vA0);
    vB1 = AE_MULFC24RA(vF1, vW1);

    vA0 = AE_SEL32_HH(vR0, vA0);
    vB0 = AE_SEL32_LL(vR0, vB1);
    vR0 = AE_SEL32_HL(vB1, vA1);

    vF0 = AE_MOVF24X2_FROMINT32X2(vA0);
    vF1 = AE_MOVF24X2_FROMINT32X2(vB0);

    AE_S32X2F24_X (vF1, p_z0, N4*sizeof(*p_z0));
    AE_S32X2F24_IP(vF0, p_z0, sizeof(*p_z0));

    for (i = 1; i < (N>>3); i++)
    {
      /*
        Y0=Y(k);
        Y1=Y(N/2+2-k);
        COSI=cosi(k);
        W1=w(k);
        W2=w(N/2+2-k);
        S=Y0+Y1;
        D=Y0-Y1;
        T0=i*real(D)+imag(S);
        T1=i*imag(D)+real(S);
        T1=T1/2;
        Y0= (imag(T0)*imag(COSI)-real(T0)*real(COSI)) + ...
              i*(real(T0)*imag(COSI)+imag(T0)*real(COSI));
        Y0=Y0/2;
        T0=T1-Y0;
        T1=T1+Y0;
      */
      AE_L32X2F24_IP(vF0, p_y0, 8);
      AE_L32X2F24_XP(vF1, p_y1, step_back);
      AE_L32X2F24_IP(vC0, p_cos, 8);
      AE_L32X2F24_IP(vW0, p_w, 8);
      AE_L32X2F24_IP(vW1, p_w, 8);

      vB0 = AE_SUBADD32S(vF0, vF1);
      vB1 = AE_ADDSUB32S(vF0, vF1);

      vF0 = AE_MOVF24X2_FROMINT32X2(vB0);
      vA0 = AE_MULFC24RA(vF0, vC0);

      vA0 = AE_SRAI32(vA0, 1);
      vA1 = AE_SRAI32(vB1, 1);

      vB0 = AE_SUB32S(vA1, vA0);
      vB1 = AE_ADD32S(vA1, vA0);

      vF0 = AE_MOVF24X2_FROMINT32X2(vB0);
      vF1 = AE_MOVF24X2_FROMINT32X2(vB1);

      /*
        z(k      )= real(T0)*real(W1)-imag(T0)*imag(W1);
        z(N+2-k  )=-real(T0)*imag(W1)-imag(T0)*real(W1);
        z(N/2+2-k)= real(T1)*real(W2)+imag(T1)*imag(W2);
        z(N/2+k  )=-real(T1)*imag(W2)+imag(T1)*real(W2);
      */
      vA0 = AE_MULFC24RA(vF0, vW0);
      vA1 = AE_NEG24S(vA0);

      vB1 = AE_MULFC24RA(vF1, vW1);

      vA1 = AE_SEL32_LL(vA1, vR0);
      vB0 = AE_SEL32_HH(vB1, vR0);
      vR0 = AE_SEL32_HL(vA0, vB1);

      vF0 = AE_MOVF24X2_FROMINT32X2(vB0);
      vF1 = AE_MOVF24X2_FROMINT32X2(vA1);
      AE_S32X2F24_X (vF0, p_z1, -N4*(int)sizeof(*p_z1));
      AE_S32X2F24_XP(vF1, p_z1, step_back);

      AE_L32X2F24_IP(vF0, p_y0, 8);
      AE_L32X2F24_XP(vF1, p_y1, step_back);
      AE_L32X2F24_IP(vC0, p_cos, 8);
      AE_L32X2F24_IP(vW0, p_w, 8);
      AE_L32X2F24_IP(vW1, p_w, 8);

      vB0 = AE_SUBADD32S(vF0, vF1);
      vB1 = AE_ADDSUB32S(vF0, vF1);

      vF0 = AE_MOVF24X2_FROMINT32X2(vB0);
      vA0 = AE_MULFC24RA(vF0, vC0);

      vA0 = AE_SRAI32(vA0, 1);
      vA1 = AE_SRAI32(vB1, 1);

      vB0 = AE_SUB32S(vA1, vA0);
      vB1 = AE_ADD32S(vA1, vA0);

      vF0 = AE_MOVF24X2_FROMINT32X2(vB0);
      vF1 = AE_MOVF24X2_FROMINT32X2(vB1);

      /*
        z(k      )= real(T0)*real(W1)-imag(T0)*imag(W1);
        z(N+2-k  )=-real(T0)*imag(W1)-imag(T0)*real(W1);
        z(N/2+2-k)= real(T1)*real(W2)+imag(T1)*imag(W2);
        z(N/2+k  )=-real(T1)*imag(W2)+imag(T1)*real(W2);
      */
      vA0 = AE_MULFC24RA(vF0, vW0);
      vA1 = AE_NEG24S(vA0);

      vB1 = AE_MULFC24RA(vF1, vW1);
 
      vA0 = AE_SEL32_HH(vR0, vA0);
      vB0 = AE_SEL32_LL(vR0, vB1);
      vR0 = AE_SEL32_HL(vB1, vA1);

      vF0 = AE_MOVF24X2_FROMINT32X2(vA0);
      vF1 = AE_MOVF24X2_FROMINT32X2(vB0);

      AE_S32X2F24_X (vF1, p_z0, N4*sizeof(*p_z0));
      AE_S32X2F24_IP(vF0, p_z0, sizeof(*p_z0));
    }

    /*** middle sample ***/
    /*
      W1=w(N/4+1);
      T0=Y(N/4+1);
      z(N/4+1  )= real(T0)*real(W1)+imag(T0)*imag(W1);
      z(N+1-N/4)=-real(T0)*imag(W1)+imag(T0)*real(W1);
    */
    AE_L32X2F24_IP(vF0, p_y0, 8);
    AE_L32X2F24_IP(vW0, p_w, 8);

    vB0 = AE_MULFC24RA(vF0, vW0);

    vA0 = AE_SEL32_HH(vB0, vR0);
    vA1 = AE_SEL32_LL(vB0, vR0);

    vF0 = AE_MOVF24X2_FROMINT32X2(vA0);
    vF1 = AE_MOVF24X2_FROMINT32X2(vA1);

    AE_S32X2F24_I(vF0, p_z0, 0);
    AE_S32X2F24_I(vF1, p_z1, 0);

    return 30-NSA(N);
} /* dct_24x24() */
