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

#include "NatureDSP_Signal_fft.h"
#include "NatureDSP_Signal_vector.h"
#include "fft_twiddles32x32.h"
#include "common.h"


inline_ void _cmult32x32(ae_int32x2 *result, ae_int32x2 *x, ae_int32x2 *y)
{
#if (XCHAL_HAVE_HIFI3Z)
    ae_f32x2 z;

    z = AE_MULFCI32RAS(AE_MOVF32X2_FROMINT32X2(*x), AE_MOVF32X2_FROMINT32X2(*y));
    AE_MULFCR32RAS(z, AE_MOVF32X2_FROMINT32X2(*x), AE_MOVF32X2_FROMINT32X2(*y));
    *result = AE_MOVINT32X2_FROMF32X2(z);
#else

    ae_int32x2 b, c, d, ac, bd;
    c = AE_SEL32_HH(*y, *y);
    d = AE_SEL32_LL(*y, *y);
    b = AE_SEL32_LH(*x, *x);
    ac = AE_MULFP32X2RAS(*x, c);
    bd = AE_MULFP32X2RAS(b, d);
    ac = AE_SUBADD32S(ac, bd);
    *result = ac;
#endif
}


#define _CONJ32(_x) {_x = AE_SEL32_HL(_x, AE_NEG32S(_x) ); }

static int SpectrConv(complex_fract32 *y, int N, const complex_fract32 *twiddle_table, int twiddle_stride, int scalingOpt, int bexp)
{
    ae_int32x2  vA0, vA1, vB0, vB1, tw;
    ae_int32x2 * restrict p_x0,
        *restrict p_x1,
        *restrict ptw = (ae_int32x2*)(twiddle_table + twiddle_stride),
        *restrict p_y0,
        *restrict p_y1;
    /*----------------------------------------------------------------------------
    Apply the in-place real-to-complex spectrum conversion

    MATLAB code:
    twd = exp(-2*pi*1j*(0:N/4-1)/N);
    a0 = y(1:N/4);
    a1 = [y(1),wrev(y(N/4+2:N/2))];
    b0 = 1/2*(a0+conj(a1));
    b1 = 1/2*(a0-conj(a1))*-1j.*twd;
    a0 = b0+b1;
    a1 = b0-b1;
    y(1:N) = [a0,conj(y(N/4+1)),wrev(conj(a1(2:N/4))), ...
    a1,     y(N/4+1) ,wrev(conj(a0(2:N/4)))];
    */
    int shift = (scalingOpt == 3) ?
        0 :
        (bexp < 1) ? 1 - bexp : 0;

    int  n;
    p_x0 = (ae_int32x2 *)(y);
    p_x1 = (ae_int32x2 *)(y + N / 2 - 1);
    p_y0 = (ae_int32x2 *)(y);
    p_y1 = (ae_int32x2 *)(y + N / 2);

    /*
    b0 = y[0];
    b0.s.re >>= shift;
    b0.s.im >>= shift;

    a0.s.re = L_add_ll(b0.s.re, b0.s.im);
    a0.s.im = 0;
    a1.s.re = L_sub_ll(b0.s.re, b0.s.im);
    a1.s.im = 0;

    a1 = conj_fr32c(a1);
    y[0] = a0;
    y[N / 2] = a1;
    */
    AE_L32X2_IP(vB0, p_x0, 8);
    vB0 = AE_SRAA32(vB0, shift);
    vB1 = AE_SEL32_HH(vB0, vB0);
    vB0 = AE_SEL32_LL(vB0, AE_NEG32S(vB0));

    vB0 = AE_ADD32S(vB0, vB1);

    vA0 = AE_SEL32_HH(vB0, AE_MOVI(0));
    vA1 = AE_SEL32_LL(vB0, AE_MOVI(0));

    AE_S32X2_IP(vA0, p_y0, sizeof(complex_fract32));
    AE_S32X2_XP(vA1, p_y1, -(int)sizeof(complex_fract32));

    /*15 cycles per pipeline stage in steady state with unroll=2*/
    __Pragma("loop_count min=3");
    for (n = 1; n < N / 4; n++)
    {
        /*  a0 = y[++k0];
        a1 = y[++k1];

        a0.s.re >>= shift;
        a0.s.im >>= shift;
        a1.s.re >>= shift;
        a1.s.im >>= shift;

        // b0 <- 1/2*(a0+conj(a1));
        b0.s.re = (fract32)(((int64_t)a0.s.re + a1.s.re) >> 1);
        b0.s.im = (fract32)(((int64_t)a0.s.im - a1.s.im) >> 1);
        // b1 <- 1/2*(a0-conj(a1))*-1j
        b1.s.re = (fract32)(((int64_t)a0.s.im + a1.s.im) >> 1);
        b1.s.im = (fract32)(((int64_t)a1.s.re - a0.s.re) >> 1);

        b0.s.re = (fract32)vB0._[0];
        b0.s.im = (fract32)vB0._[1];
        b1.s.re = (fract32)vB1._[0];
        b1.s.im = (fract32)vB1._[1];

        // b1 <- b1*twd
        b1 = mpy_CQ31CQ31(b1, twiddle_table[n * twiddle_stride]);        */

        AE_L32X2_IP(vB0, p_x0, 8);
        AE_L32X2_XP(vB1, p_x1, -8);
        AE_L32X2_XP(tw, ptw, twiddle_stride*sizeof(complex_fract32));

        vB0 = AE_SRAA32(vB0, shift);
        vB1 = AE_SRAA32(vB1, shift);

        // ADD/SUBB
        vA0 = AE_ADD32S(vB0, vB1);
        vA1 = AE_SUB32S(vB0, vB1);

        vB1 = AE_SEL32_LH(vA0, AE_NEG32S(vA1));
        vB0 = AE_SEL32_HL(vA0, vA1);

        vB1 = AE_SRAI32(vB1, 1);
        vB0 = AE_SRAI32(vB0, 1);


        _cmult32x32(&vB1, &vB1, &tw);
        /*  a0 <- b0+b1
        a0.s.re = L_add_ll(b0.s.re, b1.s.re);
        a0.s.im = L_add_ll(b0.s.im, b1.s.im);
        a1 <- b0-b1
        a1.s.re = L_sub_ll(b0.s.re, b1.s.re);
        a1.s.im = L_sub_ll(b0.s.im, b1.s.im);
        a1 = conj_fr32c(a1);
        y[k0] = a0;
        y[k1] = a1;  */

        vA0 = AE_ADD32S(vB0, vB1);
        vA1 = AE_SUB32S(vB0, vB1);
        _CONJ32(vA1);  // a1 = conj(a1)

        AE_S32X2_IP(vA0, p_y0, sizeof(complex_fract32));
        AE_S32X2_XP(vA1, p_y1, -(int)sizeof(complex_fract32));
    }

    /*
    a0 = y[++k0];
    a0.s.re >>= shift;
    a0.s.im >>= shift;

    y[k0] = conj_fr32c(a0);
    */
    AE_L32X2_IP(vA0, p_x0, 8);
    vA0 = AE_SRAA32(vA0, shift);
    _CONJ32(vA0);  /* a0 = conj(a0) */
    AE_S32X2_IP(vA0, p_y0, sizeof(complex_fract32));

    return shift;
} /* SpectrConv */

/*-------------------------------------------------------------------------
  FFT on Real Data with Optimized Memory Usage
  These functions make FFT on real data forming half of spectrum with
  optimized memory usage.
  Scaling: 
      +-----------------------+--------------------------------------+
      |      Function         |           Scaling options            |
      +-----------------------+--------------------------------------+
      |  fft_real16x16_ie     |  2 - 16-bit dynamic scaling          |
      |  fft_real32x16_ie     |  3 - fixed scaling before each stage |
      |  fft_real24x24_ie     |  3 - fixed scaling before each stage |
      |  fft_real24x24_ie_24p |  3 - fixed scaling before each stage |
      |  fft_real32x16_ie_24p |  1 - 24-bit scaling                  |
      |  fft_real32x32_ie     |  2 - 32-bit dynamic scaling          |    
      |                       |  3 - fixed scaling before each stage |   
      +-----------------------+--------------------------------------+
    
  NOTES:
  1. Bit-reversing reordering is done here.
  2. INPUT DATA MAY APPEAR DAMAGED after the call.
  3. FFT functions may use input and output buffers for temporal storage
     of intermediate 32-bit data, so FFT functions with 24-bit packed
     I/O (Nx3-byte data) require that the buffers are large enough to 
     keep Nx4-byte data.
  4. FFT of size N may be supplied with constant data (twiddle factors) 
     of a larger-sized FFT = N*twdstep.

  Precision:
  16x16_ie      16-bit input/outputs, 16-bit data, 16-bit twiddles
  24x24_ie      24-bit input/outputs, 24-bit data, 24-bit twiddles
  32x16_ie      32-bit input/outputs, 32-bit data, 16-bit twiddles
  32x32_ie      32-bit input/outputs, 32-bit data, 32-bit twiddles
  24x24_ie_24p  24-bit packed input/outputs, 24-bit data, 24-bit twiddles
  32x16_ie_24p  24-bit packed input/outputs, 32-bit data, 16-bit twiddles
  f_ie          floating point

  Input:
  x                     input signal: 
  --------------+----------+-----------------+----------+
  Function      |   Size   |  Allocated Size |  type    |
  --------------+----------+-----------------+-----------
  16x16_ie      |     N    |      N          |  int16_t |
  24x24_ie      |     N    |      N          |   f24    |
  32x16_ie      |     N    |      N          |  int32_t |
  32x32_ie      |     N    |      N          |  int32_t |
  24x24_ie_24p  |     3*N  |      4*N+8      |  uint8_t |
  32x16_ie_24p  |     3*N  |      4*N+8      |  uint8_t |
  --------------+----------+-----------------+----------+

  twd[N*twdstep*3/4]    twiddle factor table of a complex-valued 
                        FFT of size N*twdstep
  N                     FFT size
  twdstep               twiddle step
  scalingOpt            scaling option (see table above), not applicable 
                        to the floating point function

  Output:
  y                     output spectrum. Real and imaginary data 
                        are interleaved and real data goes first:
  --------------+----------+-----------------+---------------+
  Function      |   Size   |  Allocated Size |  type         |
  --------------+----------+-----------------+----------------
  16x16_ie      |   N/2+1  |      N/2+1      |complex_fract16|
  24x24_ie      |   N/2+1  |      N/2+1      |complex_fract32|
  32x16_ie      |   N/2+1  |      N/2+1      |complex_fract32|
  32x32_ie      |   N/2+1  |      N/2+1      |complex_fract32|
  24x24_ie_24p  |  3*(N+2) |      4*N+8      |  uint8_t      |
  32x16_ie_24p  |  3*(N+2) |      4*N+8      |  uint8_t      |
  f_ie          |   N/2+1  |      N/2+1      |complex_float  |
  --------------+----------+-----------------+---------------+

  Returned value: total number of right shifts occurred during scaling
  procedure

  Restrictions:
  x,y     should not overlap
  x,y     aligned on 8-bytes boundary
-------------------------------------------------------------------------*/

int fft_real32x32_ie(complex_fract32* y, int32_t  * x, const complex_fract32* twd, int twdstep, int N, int scalingOpt)
{
    int shift ; 
    int bexp; 

    NASSERT(scalingOpt==2 || scalingOpt==3); 
    NASSERT((void*)x != (void*)y);
    NASSERT_ALIGN8(x); 
    NASSERT_ALIGN8(y);
    NASSERT_ALIGN8(twd);

    shift = fft_cplx32x32_ie(y, (complex_fract32*)x, twd, twdstep*2, N/2, scalingOpt);

    if (scalingOpt == 2)
    {
        bexp = vec_bexp32((int32_t*)y, N);
    }
    else
    {
        bexp = 0;
    }

    shift += SpectrConv(y, N, twd, 3*twdstep, scalingOpt, bexp);

    return shift;  
}
