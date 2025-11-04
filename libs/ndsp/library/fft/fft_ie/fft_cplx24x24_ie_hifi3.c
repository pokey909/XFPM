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
    FFT on Complex Data with Optimized Memory Usage
    C code optimized for HiFi3
	Integrit, 2006-2017
*/

#include "NatureDSP_Signal_fft.h"
#include "common.h"

extern int fft_stage_last_ie( int32_t *x, 
                               int32_t *y, 
                               int N); 


/*-------------------------------------------------------------------------
  FFT on Complex Data with Optimized Memory Usage
  These functions make FFT on complex data with optimized memory usage.
  Scaling: 
      +-------------------+----------------------------------------+
      |      Function     |           Scaling options              |
      +-------------------+----------------------------------------+
      |  fft_cplx16x16_ie |  2 - 16-bit dynamic scaling            | 
      |  fft_cplx24x24_ie |  3 - fixed scaling before each stage   | 
      |  fft_cplx32x16_ie |  3 - fixed scaling before each stage   | 
      |  fft_cplx32x32_ie |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      +-------------------+----------------------------------------+
  NOTES:
  1. Bit-reversing reordering is done here.
  2. FFT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED after 
     the call.
  3. FFT of size N may be supplied with constant data
     (twiddle factors) of a larger-sized FFT = N*twdstep.

  Precision: 
  16x16_ie      16-bit input/outputs, 16-bit twiddles
  24x24_ie      24-bit input/outputs, 24-bit twiddles
  32x16_ie      32-bit input/outputs, 16-bit twiddles
  32x32_ie      32-bit input/outputs, 32-bit twiddles
  f_ie          floating point
 
  Input:
  x[N]                  complex input signal. Real and imaginary data 
                        are interleaved and real data goes first
  twd[N*twdstep*3/4]    twiddle factor table of a complex-valued FFT of 
                        size N*twdstep
  N                     FFT size
  twdstep               twiddle step 
  scalingOpt            scaling option (see table above), not applicable
                        to the floating point function 
  Output:
  y[N]                  output spectrum. Real and imaginary data are 
                        interleaved and real data goes first

  Returned value: total number of right shifts occurred during scaling 
                  procedure. Floating point function always return 0.

  Restrictions:
  x,y   should not overlap
  x,y   aligned on 8-bytes boundary
-------------------------------------------------------------------------*/
int fft_cplx24x24_ie(complex_fract32* y,complex_fract32* x, const complex_fract32* twd, int twdstep, int N, int scalingOpt)
#if 0
{
    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT(scalingOpt==3);
    NASSERT(N==128||N==256||N==512||N==1024);

    int shift = 0; 
    int s;    
    int stride = N/4;     
    int M = N*twdstep; 

    ae_f24x2   * restrict py24; 
    ae_f24x2   * restrict px24; 

    ae_f24x2   * restrict ptw1 = (ae_f24x2*)  twd; 
    ae_f24x2   * restrict ptw2 = (ae_f24x2*) (twd+1*M/4);
    ae_f24x2   * restrict ptw3 = (ae_f24x2*) (twd+2*M/4);
 
    ae_int32x2  vA0, vA1, vA2, vA3, vB0, vB1, vB2, vB3, vC0, vC1, vC2, vC3;

    ae_f24x2    vT1, vT2, vT3;
    int i; 
    int log2n = 0; 

    WUR_AE_CBEGIN0((unsigned)x); 
    WUR_AE_CEND0((unsigned)(uintptr_t)(x+(N-1)));     

    s = 3; //Scaling for first stage
    while( stride > 1 )
    {
        
        WUR_AE_SAR(s); 
        unsigned acc = 0; 
        const unsigned acc_inc  = (log2n==0)? 0: (0x80000000 >> (log2n-1)); 

        shift += s;       
        log2n += 2;   
        px24 = py24 = (ae_f24x2*)x; 

        ptw1 = (ae_f24x2*)  twd; 
        ptw2 = (ae_f24x2*) (twd+1*M/4);
        ptw3 = (ae_f24x2*) (twd+2*M/4);

        i = N/4;
        //12 cycles per pipeline stage
        do 
        {
            ae_f24x2 t1, t2, t3, t0;
            int offset_inc = 0; 
            acc+=acc_inc;
            XT_MOVEQZ(offset_inc, twdstep*8, acc); 
            

            vA3 = AE_L32X2F24_X ( px24, 3*8*stride); 
            vA2 = AE_L32X2F24_X ( px24, 2*8*stride); 
            vA1 = AE_L32X2F24_X ( px24, 1*8*stride); 
            AE_L32X2F24_XC(t0, px24,   4*8*stride); 
            vA0 = t0; 

            vB0 = AE_ADD32S(vA0, vA2);
            vB2 = AE_SUB32S(vA0, vA2);
            vB1 = AE_ADD32S(vA1, vA3);
            vB3 = AE_SUB32S(vA1, vA3);
            
            vB3 = AE_SEL32_LH(vB3, vB3);
            
            vC0 = AE_ADD32S(vB0, vB1);
            vC2 = AE_SUB32S(vB0, vB1);
            vC3 = AE_SUBADD32S(vB2, vB3); 
            vC1 = AE_ADDSUB32S(vB2, vB3); 
    #if 1        
            vC0 = AE_SRAS32(vC0); 
            vC2 = AE_SRAS32(vC2); 
            vC3 = AE_SRAS32(vC3); 
            vC1 = AE_SRAS32(vC1); 
    #else
            // Better precision, 
            // but 16 cycles per pipeline stage
            vC0 = AE_SRAA32RS(vC0, s); 
            vC2 = AE_SRAA32RS(vC2, s); 
            vC3 = AE_SRAA32RS(vC3, s); 
            vC1 = AE_SRAA32RS(vC1, s); 
    #endif
            AE_L32X2F24_XP(vT1, ptw1, offset_inc);
            AE_L32X2F24_XP(vT2, ptw2, offset_inc);
            AE_L32X2F24_XP(vT3, ptw3, offset_inc);

            t1 =  AE_MOVF24X2_FROMINT32X2(vC1);
            t2 =  AE_MOVF24X2_FROMINT32X2(vC2);
            t3 =  AE_MOVF24X2_FROMINT32X2(vC3);

            vC1 = AE_MULFC24RA(t1, vT1);
            vC2 = AE_MULFC24RA(t2, vT2); 
            vC3 = AE_MULFC24RA(t3, vT3); 

            t0 = AE_MOVF24X2_FROMINT32X2(vC0);
            t1 = AE_MOVF24X2_FROMINT32X2(vC1);
            t2 = AE_MOVF24X2_FROMINT32X2(vC2);
            t3 = AE_MOVF24X2_FROMINT32X2(vC3);
        
            AE_S32X2F24_X(t3,  py24, 3*8*stride);
            AE_S32X2F24_X(t2,  py24, 1*8*stride);
            AE_S32X2F24_X(t1,  py24, 2*8*stride);
            AE_S32X2F24_XC(t0, py24, 4*8*stride);     

        } while(--i);

        s = 2;  //Scaling for other stages
        stride>>=2;  
        twdstep<<=2;
    }    //while( stride > 1 )

    shift += fft_stage_last_ie((int32_t*)x, (int32_t*)y, N); 

    return shift;
} /* fft_cplx24x24_ie() */
#else
{
    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT(scalingOpt==3);
    NASSERT(N==128||N==256||N==512||N==1024);

    int shift = 0; 
    int stride = N/4;     
    int M = N*twdstep; 

          ae_f24x2   * restrict py; 
    const ae_f24x2   * restrict px; 
          ae_f24x2   * restrict py1; 
    const ae_f24x2   * restrict px1; 
    ae_f24x2   * restrict ptw1 = (ae_f24x2*)  twd; 
    ae_f24x2   * restrict ptw2 = (ae_f24x2*) (twd+1*M/4);
    ae_f24x2   * restrict ptw3 = (ae_f24x2*) (twd+2*M/4);
    int32_t acc = 0;
    uint32_t acc_inc; 
    ae_int32x2  vA0, vA1, vA2, vA3, vB0, vB1, vB2, vB3, vC0, vC1, vC2, vC3;
    ae_f24x2    vT1, vT2, vT3;
    int i; 

    WUR_AE_CBEGIN0((unsigned)x); 
    WUR_AE_CEND0((unsigned)(uintptr_t)(x+(N-1)));     

    // on the first stage we are reading all the twiddles, on all subsequent 
    // stages we will read just each of 2
    {
        WUR_AE_SAR(3); 
        shift = 3;       
        px = py = (ae_f24x2*)x; 

        ptw1 = (ae_f24x2*)  twd; 
        ptw2 = (ae_f24x2*) (twd+1*M/4);
        ptw3 = (ae_f24x2*) (twd+2*M/4);

        for (i=0; i<(N>>2); i++)
        {
            ae_f24x2 t1, t2, t3, t0;
            int offset_inc; 
            offset_inc=twdstep*8; 

            vA3 = AE_L32X2F24_X ( px, 3*8*stride); 
            vA2 = AE_L32X2F24_X ( px, 2*8*stride); 
            vA1 = AE_L32X2F24_X ( px, 1*8*stride); 
            AE_L32X2F24_XC(t0, px,   4*8*stride); 
            vA0 = t0; 

            vB0 = AE_ADD32S(vA0, vA2);
            vB2 = AE_SUB32S(vA0, vA2);
            vB1 = AE_ADD32S(vA1, vA3);
            vB3 = AE_SUB32S(vA1, vA3);
            
            vB3 = AE_SEL32_LH(vB3, vB3);
            
            vC0 = AE_ADD32S(vB0, vB1);
            vC2 = AE_SUB32S(vB0, vB1);
            vC3 = AE_SUBADD32S(vB2, vB3); 
            vC1 = AE_ADDSUB32S(vB2, vB3); 
            vC0 = AE_SRAS32(vC0); 
            vC2 = AE_SRAS32(vC2); 
            vC3 = AE_SRAS32(vC3); 
            vC1 = AE_SRAS32(vC1); 
            AE_L32X2F24_XP(vT1, ptw1, offset_inc);
            AE_L32X2F24_XP(vT2, ptw2, offset_inc);
            AE_L32X2F24_XP(vT3, ptw3, offset_inc);

            t1 =  AE_MOVF24X2_FROMINT32X2(vC1);
            t2 =  AE_MOVF24X2_FROMINT32X2(vC2);
            t3 =  AE_MOVF24X2_FROMINT32X2(vC3);

            vC1 = AE_MULFC24RA(t1, vT1);
            vC2 = AE_MULFC24RA(t2, vT2); 
            vC3 = AE_MULFC24RA(t3, vT3); 

            t0 = AE_MOVF24X2_FROMINT32X2(vC0);
            t1 = AE_MOVF24X2_FROMINT32X2(vC1);
            t2 = AE_MOVF24X2_FROMINT32X2(vC2);
            t3 = AE_MOVF24X2_FROMINT32X2(vC3);
        
            AE_S32X2F24_X(t3,  py, 3*8*stride);
            AE_S32X2F24_X(t2,  py, 1*8*stride);
            AE_S32X2F24_X(t1,  py, 2*8*stride);
            AE_S32X2F24_XC(t0, py, 4*8*stride);     
        }
        stride>>=2;  
        twdstep<<=2;
    }
    // intermediate stages
    acc_inc  = 0x80000000; 
    while( stride > 1 )
    {
        __Pragma("no_reorder")
        acc = 0; 
        WUR_AE_SAR(2); 
        shift += 2;       
        py = (ae_f24x2*)x; 
        py1=py+4*stride;
        ptw1 = (ae_f24x2*)  twd; 
        ptw2 = (ae_f24x2*) (twd+1*M/4);
        ptw3 = (ae_f24x2*) (twd+2*M/4);

        for (i=0; i<(N>>3); i++)
        {
            ae_f24x2 t1, t2, t3, t0;
            int offset_inc = 0; 
            acc+=acc_inc;
            XT_MOVEQZ(offset_inc, twdstep*8, acc); 
            px=(const ae_f24x2*)py;
            vA3 = AE_L32X2F24_X ( px, 3*8*stride); 
            vA2 = AE_L32X2F24_X ( px, 2*8*stride); 
            vA1 = AE_L32X2F24_X ( px, 1*8*stride); 
            t0  = AE_L32X2F24_I ( px, 0); 
            vA0 = t0; 

            vB0 = AE_ADD32S(vA0, vA2);
            vB2 = AE_SUB32S(vA0, vA2);
            vB1 = AE_ADD32S(vA1, vA3);
            vB3 = AE_SUB32S(vA1, vA3);
            
            vB3 = AE_SEL32_LH(vB3, vB3);
            
            vC0 = AE_ADD32S(vB0, vB1);
            vC2 = AE_SUB32S(vB0, vB1);
            vC3 = AE_SUBADD32S(vB2, vB3); 
            vC1 = AE_ADDSUB32S(vB2, vB3); 
            vC0 = AE_SRAS32(vC0); 
            vC2 = AE_SRAS32(vC2); 
            vC3 = AE_SRAS32(vC3); 
            vC1 = AE_SRAS32(vC1); 
            vT2=AE_L32X2F24_X  (ptw1, 1*M*2);
            vT3=AE_L32X2F24_X  (ptw1, 2*M*2);
            AE_L32X2F24_XP(vT1, ptw1, offset_inc);

            t1 =  AE_MOVF24X2_FROMINT32X2(vC1);
            t2 =  AE_MOVF24X2_FROMINT32X2(vC2);
            t3 =  AE_MOVF24X2_FROMINT32X2(vC3);

            vC1 = AE_MULFC24RA(t1, vT1);
            vC2 = AE_MULFC24RA(t2, vT2); 
            vC3 = AE_MULFC24RA(t3, vT3); 

            t0 = AE_MOVF24X2_FROMINT32X2(vC0);
            t1 = AE_MOVF24X2_FROMINT32X2(vC1);
            t2 = AE_MOVF24X2_FROMINT32X2(vC2);
            t3 = AE_MOVF24X2_FROMINT32X2(vC3);
        
            AE_S32X2F24_X(t3,  py,  3*8*stride);
            AE_S32X2F24_X(t2,  py,  1*8*stride);
            AE_S32X2F24_X(t1,  py,  2*8*stride);
            AE_S32X2F24_XC(t0, py,2*4*8*stride);     

            px1=(const ae_f24x2*)py1;
            vA3 = AE_L32X2F24_X ( px1, 3*8*stride); 
            vA2 = AE_L32X2F24_X ( px1, 2*8*stride); 
            vA1 = AE_L32X2F24_X ( px1, 1*8*stride); 
            t0  = AE_L32X2F24_I ( px1, 0); 
            vA0 = t0; 

            vB0 = AE_ADD32S(vA0, vA2);
            vB2 = AE_SUB32S(vA0, vA2);
            vB1 = AE_ADD32S(vA1, vA3);
            vB3 = AE_SUB32S(vA1, vA3);
            
            vB3 = AE_SEL32_LH(vB3, vB3);
            
            vC0 = AE_ADD32S(vB0, vB1);
            vC2 = AE_SUB32S(vB0, vB1);
            vC3 = AE_SUBADD32S(vB2, vB3); 
            vC1 = AE_ADDSUB32S(vB2, vB3); 
            vC0 = AE_SRAS32(vC0); 
            vC2 = AE_SRAS32(vC2); 
            vC3 = AE_SRAS32(vC3); 
            vC1 = AE_SRAS32(vC1); 

            t1 =  AE_MOVF24X2_FROMINT32X2(vC1);
            t2 =  AE_MOVF24X2_FROMINT32X2(vC2);
            t3 =  AE_MOVF24X2_FROMINT32X2(vC3);

            vC1 = AE_MULFC24RA(t1, vT1);
            vC2 = AE_MULFC24RA(t2, vT2); 
            vC3 = AE_MULFC24RA(t3, vT3); 

            t0 = AE_MOVF24X2_FROMINT32X2(vC0);
            t1 = AE_MOVF24X2_FROMINT32X2(vC1);
            t2 = AE_MOVF24X2_FROMINT32X2(vC2);
            t3 = AE_MOVF24X2_FROMINT32X2(vC3);
        
            AE_S32X2F24_X(t3,  py1,  3*8*stride);
            AE_S32X2F24_X(t2,  py1,  1*8*stride);
            AE_S32X2F24_X(t1,  py1,  2*8*stride);
            AE_S32X2F24_XC(t0, py1,2*4*8*stride);     
        }
        stride>>=2;  
        twdstep<<=2;
        acc_inc>>=2; 
    }

    shift += fft_stage_last_ie((int32_t*)x, (int32_t*)y, N); 

    return shift;
} /* fft_cplx24x24_ie() */
#endif
