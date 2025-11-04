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
#include "NatureDSP_Signal_fft.h"
#include "NatureDSP_Signal_vector.h"
#include "common.h"
#include "fft_twiddles32x32.h"


#define _CONJ32(_x) {_x = AE_SEL32_HL(_x, AE_NEG32S(_x) ); }

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
 /* 
    The real - to - complex spectrum conversion 
    MATLAB code:
    % N - Size of transform
    % X - input complex vector 1 x (N/2+1 )
    % x = real(N*ifft([X, conj(wrev(X(2:N/2))] ) )- output real vector 1xN
    
    twd = exp(-2*pi*1j*(0:N4-1)/N);
    a0 = X(1:N4);
    a1 = X(N/2+1:-1:N/2-N4+2); 
    b0 = a0+conj(a1);
    b1 = (a0-conj(a1))*1j.*conj(twd);
    a0 = b0+b1;
    a1 = conj(b0-b1);
    if(mod(N,4))
        x = [a0,  wrev(a1(2:N4))]; % N/2 complex samples
    else
        x = [a0,2*conj(X(N4+1)),wrev(a1(2:N4))]; % N/2 complex samples
    end

    tmp = N/2*ifft(x); 
    x = zeros(1, N);
    x(1:2:end) = real(tmp);
    x(2:2:end) = imag(tmp);
 */
static int iSpectrConv(complex_fract32 *x, int N, const complex_fract32 *twiddle_table, int twiddle_stride, int scalingOpt, int bexp)
{
    const int N4 = (N / 2 + 1) >> 1; /* Works for all even N */
    ae_int32x2  vA0, vA1, vB0, vB1, tw;
    ae_int32x2 * restrict p_x0,
        *restrict p_x1,
        *restrict ptw = (ae_int32x2*)(twiddle_table + twiddle_stride),
        *restrict p_y0,
        *restrict p_y1;
    int n;
    int shift = (scalingOpt == 3) ? 2 : 2 - bexp;
    ASSERT(shift>-32 && shift<32);

    p_x0 = (ae_int32x2 *)(x);
    p_x1 = (ae_int32x2 *)(x + N / 2 );
    p_y0 = (ae_int32x2 *)(x);
    p_y1 = (ae_int32x2 *)(x + N / 2 - 1);

    AE_L32X2_IP(vB0, p_x0, 8);
    AE_L32X2_XP(vB1, p_x1, -8);

    vB0 = AE_SRAA32(vB0, shift);
    vB1 = AE_SRAA32(vB1, shift);
    vA0 = AE_ADD32S(vB0, vB1);
    vA1 = AE_SUB32S(vB0, vB1);

    vA0 = AE_SEL32_HH(vA0, vA1);
    AE_S32X2_IP(vA0, p_y0, sizeof(complex_fract32)); 

    __Pragma("loop_count min=3");
    for (n = 1; n < N4; /*k0++, k1--, */n++)
    {
        /* 15 cycles per pipeline stage in steady state with unroll=2 */
        AE_L32X2_IP(vB0, p_x0, 8);
        AE_L32X2_XP(vB1, p_x1, -8);
        AE_L32X2_XP(tw, ptw, twiddle_stride*sizeof(complex_fract32));

        vB0 = AE_SRAA32(vB0, shift);
        vB1 = AE_SRAA32(vB1, shift);

        // ADD/SUBB
        vA0 = AE_ADD32S(vB0, vB1);
        vA1 = AE_SUB32S(vB0, vB1);

        vB0 = AE_SEL32_HL(vA0, vA1);
        vB1 = AE_SEL32_LH(AE_NEG32S(vA0), vA1);

        _CONJ32(tw); 
        _cmult32x32(&vB1, &vB1, &tw);
        
        vA0 = AE_ADD32S(vB0, vB1); 
        vA1 = AE_SUB32S(vB0, vB1); 
        _CONJ32(vA1); 

        AE_S32X2_IP(vA0, p_y0, sizeof(complex_fract32)); 
        AE_S32X2_XP(vA1, p_y1, -(int)sizeof(complex_fract32));
    }

    if (N & 3)
    {
        /* When N is not multiple of 4 */
        return shift;
    }

    /* 2*conj(x(N/4+1)) */
    AE_L32X2_IP(vB0, p_x0, 8);
    vB0 = AE_SRAA32(vB0, shift-1);
    _CONJ32(vB0); 
    AE_S32X2_IP(vB0, p_y0, sizeof(complex_fract32));

    return shift;
} /* iSpectrConv */

/*-------------------------------------------------------------------------
  Inverse FFT on Real Data
  These functions make inverse FFT on half spectral data forming real
  data samples.
      Scaling: 
      +-------------------+----------------------------------------+
      |      Function     |           Scaling options              |
      +-------------------+----------------------------------------+
      |  ifft_real16x16   |  2 - 16-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  ifft_real32x32   |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  ifft_real32x16   |  3 - fixed scaling before each stage   | 
      |  ifft_real24x24   |  0 - no scaling                        | 
      |                   |  1 - 24-bit scaling                    |
      |                   |  2 - 32-bit scaling on the first stage |
      |                   |  and 24-bit scaling later              |
      |                   |  3 - fixed scaling before each stage   |
      +-------------------+----------------------------------------+

  NOTES:
  1. Bit-reversing reordering is done here. 
  2. IFFT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED after
     the call.
  3. Inverse FFT function for real signal transforms the input spectrum  
     and then calls ifft_cplx() with FFT size set to N/2.
  4. 32x32 FFTs support mixed radix transforms
  5. N - FFT size

  Precision:
  32x32  32-bit input/outputs, 32-bit twiddles
  24x24  24-bit input/outputs, 24-bit twiddles
  32x16  32-bit input/outputs, 16-bit twiddles
  16x16  16-bit input/outputs, 16-bit twiddles

  Input:
  x[(N/2+1)*2]	input spectrum. Real and imaginary data are interleaved  
                and real data goes first
  scalingOpt	scaling option (see table above)

  Output:			
  y[N]	        real output signal

  Returned value: total number of right shifts occurred during scaling 
                  procedure

  Restrictions:
  x,y           should not overlap
  x,y           aligned on a 8-bytes boundary
-------------------------------------------------------------------------*/
int ifft_real32x32(int32_t* y, int32_t* x, fft_handle_t h, int scalingOpt)
{
    int shift;
    int N, bexp;
    fft_real32x32_descr_t *hr = (fft_real32x32_descr_t *)h;
    fft_cplx32x32_descr_t *hc = (fft_cplx32x32_descr_t *)hr->cfft_hdl;
    N = 2 * hc->N;

    NASSERT(scalingOpt == 2 || scalingOpt == 3);
    NASSERT(x != y);
    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);

    if (scalingOpt == 2)
    {
        bexp = vec_bexp32(x, N + 2);
    }
    else
    {
        bexp = 0;
    }

    shift = iSpectrConv((complex_fract32*)x, N, (complex_fract32*)hr->twd, 1, scalingOpt, bexp);
    shift += ifft_cplx32x32(y, x, hr->cfft_hdl, scalingOpt);

    return shift;
} /* ifft_real32x32 */

