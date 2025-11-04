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
    Inverse FFT on Real Data with Optimized Memory Usage
    C code optimized for HiFi3
	Integrit, 2006-2017
*/

#include "NatureDSP_Signal_fft.h"
#include "common.h"

extern void fft_pack32to24_ie(uint8_t *x,  uint8_t *y, int N); 
extern     int fft_unpack24to32_s1_ie(uint8_t *x, uint8_t *y, int N); 
extern int fft_cplx_24x24_s1_ie(f24* y, f24* x, const     f24* twd, int twdstep, int N, int exp); 
extern void fft_revorder_ie(int32_t *x, int N); 

/*
in-place inverse split part of FFT:
x[N+2]  input (N+2 samples)/output(N samples)
N       size of FFT
*/
static int isplitPart_x2_24x24_s1(f24 *X, int N, f24 *tw, int step, int shift)
{
    int i;

    const int step_back = -8;
    ae_int32x2  maxV = 0, minV = 0;
    ae_f24x2 * restrict p_x0, *restrict p_x1;
    const ae_f24x2 * restrict p_twd;

    ae_int32x2  vA0, vA1, vB0, vB1, vC0, vC1;
    ae_f24x2    vT, t0, t1;
    ae_f32x2    vF0, vF1;

    NASSERT_ALIGN8(X);  

    step *= 2 * sizeof(f24);

    p_twd = (const ae_f24x2 *)tw;
    p_x0 = (ae_f24x2 *)X;
    p_x1 = (ae_f24x2 *)(X + N);

    WUR_AE_SAR(shift);

    // first point
    vA0 = AE_L32X2F24_I(p_x0, 0);
    vA1 = AE_L32X2F24_I(p_x1, 0);

    vB0 = AE_ADD32S(vA0, vA1);
    vB1 = AE_SUB32S(vA0, vA1);
    // Additional scaling need for prevent overflow
    vB0 = AE_SRAS32(vB0);
    vB1 = AE_SRAS32(vB1);

    vB0 = AE_SEL32_HH(vB0, vB1);
    vB1 = 0;
    vA1 = 1;
    vB0 = AE_ADD32S(vB0, vA1);
    vB0 = AE_SRAI32(vB0, 1);


    t0 = AE_MOVF24X2_FROMINT32X2(vB0);    t1 = AE_MOVF24X2_FROMINT32X2(vB1);
    minV = AE_MIN32(minV, vB0); maxV = AE_MAX32(maxV, vB0);
    minV = AE_MIN32(minV, vB1); maxV = AE_MAX32(maxV, vB1);

    AE_S32X2F24_IP(t0, p_x0, 8);
    AE_S32X2F24_XP(t1, p_x1, step_back);

    vF0 = AE_L32X2F24_I(p_x0, 0);
    vF1 = AE_L32X2F24_I(p_x1, 0);

    AE_L32X2F24_XP(vT, p_twd, step);

    //  10 cycles per pipeline stage in steady state with unroll=1
    for (i = 1; i < N >> 2; i++)
    {
        ae_int32x2 tmp;
        // load twiddle
        AE_L32X2F24_XP(vT, p_twd, step);

        // conj (tw)
        tmp = AE_MOVINT32X2_FROMF24X2(vT);
        tmp = AE_ADDSUB32S(0, tmp);
        tmp = AE_SEL32_LH(tmp, tmp); 
        vT = AE_MOVF24X2_FROMINT32X2(tmp);

        vA0 = (vF0);
        vA1 = (vF1);

        // ADD/SUBB
        vB0 = AE_ADD32S(vA0, vA1);
        vB1 = AE_SUB32S(vA0, vA1);

        // Additional scaling need for prevent overflow
        vB0 = AE_SRAS32(vB0);
        vB1 = AE_SRAS32(vB1);

        vA0 = AE_SEL32_LH(vB0, vB1);
        vB1 = AE_SEL32_HL(vB0, vB1);

        // do rotation        
        vB0 = AE_MULFC24RA( AE_MOVF24X2_FROMINT32X2(vA0), vT);
        vB0 = AE_NEG32S(vB0);
        vB0 = AE_SEL32_LH(vB0, vB0); 

        // load next data
        vF0 = AE_L32X2F24_I(p_x0, 8);
        vF1 = AE_L32X2F24_X(p_x1, step_back);

        // SUM/DIFF
        vC0 = AE_ADD32S(vB1, vB0);
        vC1 = AE_ADDSUB32S(0, AE_SUB32S(vB1, vB0));

        vC0 = AE_SRAA32RS(vC0, 1);
        vC1 = AE_SRAA32RS(vC1, 1);

        minV = AE_MIN32(minV, vC0); maxV = AE_MAX32(maxV, vC0);
        minV = AE_MIN32(minV, vC1); maxV = AE_MAX32(maxV, vC1);

        t0 = AE_MOVF24X2_FROMINT32X2(vC0);
        AE_S32X2F24_IP(t0, p_x0, 8);
        t1 = AE_MOVF24X2_FROMINT32X2(vC1);
        AE_S32X2F24_XP(t1, p_x1, step_back);
    }

    // Additional scaling need for prevent overflow
    vA0 = AE_SRAS32(vF0);
    vA1 = AE_SRAS32(vF1);

    // middle sample
    vB0 = AE_NEG32S(vA0);
    vC0 = AE_SEL32_HL(vA0, vB0);
    minV = AE_MIN32(minV, vC0); maxV = AE_MAX32(maxV, vC0);
    t0 = AE_MOVF24X2_FROMINT32X2(vC0);
    AE_S32X2F24_I(t0, p_x0, 0);

    maxV = AE_MAX32( maxV, AE_NEG32S(minV) );
    maxV = AE_MAX32( maxV, AE_SEL32_LH(maxV, maxV));

    return AE_NSAZ32_L(maxV);
} /* isplitPart_x2_24x24_s1() */


/*-------------------------------------------------------------------------
  Inverse FFT on Real Data with Optimized Memory Usage
  These functions make inverse FFT on real data from half of spectrum with
  optimized memory usage.
  Scaling: 
      +-----------------------+--------------------------------------+
      |      Function         |           Scaling options            |
      +-----------------------+--------------------------------------+
      | ifft_real16x16_ie     |  2 - 16-bit dynamic scaling          |
      | ifft_real32x16_ie     |  3 - fixed scaling before each stage |
      | ifft_real24x24_ie     |  3 - fixed scaling before each stage |
      | ifft_real24x24_ie_24p |  3 - fixed scaling before each stage |
      | ifft_real32x16_ie_24p |  1 - 24-bit scaling                  |
      | ifft_real32x32_ie     |  2 - 32-bit dynamic scaling          |    
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
  x                     input spectrum (positive side). Real and imaginary
                        data are interleaved and real data goes first:
  --------------+----------+-----------------+----------------
  Function      |   Size   |  Allocated Size |       type    |
  --------------+----------+-----------------+----------------
  16x16_ie      |   N/2+1  |      N/2+1      |complex_fract16|
  24x24_ie      |   N/2+1  |      N/2+1      |complex_fract32|
  32x16_ie      |   N/2+1  |      N/2+1      |complex_fract32|
  32x32_ie      |   N/2+1  |      N/2+1      |complex_fract32|
  24x24_ie_24p  |   3*(N+2)|      4*N+8      |       uint8_t |
  32x16_ie_24p  |   3*(N+2)|      4*N+8      |       uint8_t |
  f_ie          |   N/2+1  |      N/2+1      | complex_float |
  --------------+----------+-----------------+----------------

  twd[2*N*twdstep*3/4]  twiddle factor table of a complex-valued FFT
                        of size N*twdstep
  N                     FFT size
  twdstep               twiddle step
  scalingOpt            scaling option (see table above), not applicable 
                        to the floating point function
  Output:
  y                     output spectrum. Real and imaginary data are 
                        interleaved and real data goes first:
  --------------+----------+-----------------+-----------
  Function      |   Size   |  Allocated Size |  type    |
  --------------+----------+-----------------+-----------
  16x16_ie      |     N    |      N          |  int16_t |
  24x24_ie      |     N    |      N          |   f24    |
  32x16_ie      |     N    |      N          |  int32_t |
  32x32_ie      |     N    |      N          |  int32_t |
  24x24_ie_24p  |    3*N   |      4*N+8      |  uint8_t |
  32x16_ie_24p  |    3*N   |      4*N+8      |  uint8_t |
  f_ie          |      N   |      N          | float32_t|
  --------------+----------+-----------------+-----------

  Returned value: total number of right shifts occurred during scaling
  procedure

  Restrictions:
  x,y   should not overlap
  x,y   aligned on 8-bytes boundary
-------------------------------------------------------------------------*/

int ifft_real24x24_ie_24p(  uint8_t* y,       uint8_t * x, const complex_fract32* twd, int twdstep, int N, int scalingOpt)
{
    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT(scalingOpt==1);
    NASSERT(N==256||N==512||N==1024);
#if 1
    int exp; 
    int shift; 
    exp = fft_unpack24to32_s1_ie(x, y, N + 2);
    shift = 2 + 8 - exp; 

    // Scaling of isplitPart_x2_24x24 is shift+1
    exp = isplitPart_x2_24x24_s1((f24*)y, N, (f24*)twd, twdstep, shift);
    fft_revorder_ie((int32_t*)y, N/2);
    shift += fft_cplx_24x24_s1_ie((f24*)x, (f24*)y, (const int32_t*)twd, twdstep * 2, N / 2, exp);
    fft_pack32to24_ie(x, y, N); 
    return shift + 1; 
#else
    int s; 
    fft_unpack24to32_ie(x, y, N+2); 
    s = ifft_real24x24_ie((f24*)x, (f24*)y, twd, twdstep, N, scalingOpt);     
    fft_pack32to24_ie(x, y, N); 
    return s;
#endif

} /* ifft_real24x24_ie_24p() */
