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
    Modified Discrete Cosine Transform 
    C code optimized for HiFi3
   Integrit, 2006-2017
*/
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fft.h"
#include "dct4_twd.h"
#include "common.h"

/*-------------------------------------------------------------------------
  Modified Discrete Cosine Transform.
  These functions apply Modified DCT to input (convert 2N real data to N 
  spectral components) and make inverse conversion forming 2N numbers from 
  N inputs. Normally, combination of MDCT and DCT is invertible if applied 
  to subsequent data blocks with overlapping.
  Scaling:
      +-----------------------+--------------------------------------+
      |      Function         |           Scaling options            |
      +-----------------------+--------------------------------------+
      |       mdct_24x24      |  3 - fixed scaling before each stage |
      |       mdct_32x16      |  3 - fixed scaling before each stage |
      |       mdct_32x32      |  3 - fixed scaling before each stage |
      |      imdct_24x24      |  3 - fixed scaling before each stage |
      |      imdct_32x16      |  3 - fixed scaling before each stage |
      |      imdct_32x32      |  3 - fixed scaling before each stage |
      +-----------------------+--------------------------------------+
  NOTES:
     1. MDCT/IMDCT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED 
     after the call.
     2. N - MDCT size (depends on selected MDCT handle)

  Precision: 
  24x24  24-bit input/outputs, 24-bit twiddles
  32x16  32-bit input/outputs, 16-bit twiddles
  32x32  32-bit input/outputs, 32-bit twiddles

  -------------------------------------------------------------------------
  For MDCT:
  Input:
  x[2*N]      input signal
  h           MDCT handle
  scalingOpt  scaling option (see table above)
  Output:
  y[N]        output of transform 
  -------------------------------------------------------------------------
  For IMDCT:
  Input:
  x[N]        input signal
  h           IMDCT handle
  scalingOpt  scaling option (see table above)
  Output:
  y[2*N]      output of transform
  -------------------------------------------------------------------------
  Returned value:
              total number of right shifts occurred during scaling 
              procedure 
  Restriction:
  x,y         should not overlap
  x,y         aligned on 8-bytes boundary
-------------------------------------------------------------------------*/
int mdct_24x24(f24* y, f24* x, dct_handle_t h, int scalingOpt)
{  
    const tdct4_twd_fr32 *ptwd=(const tdct4_twd_fr32 *)h;
    const ae_int32x2 * pxrd_pos0;
    const ae_int64   * pxrd_neg0;
    const ae_int32x2 * pxrd_pos1;
    const ae_int64   * pxrd_neg1;
          ae_f24x2   * restrict pxwr_pos;
          ae_f24x2   * restrict pxwr_neg;
    ae_int32x2 t0, t1, t2, t3, t4, t5, t6, t7;
    ae_f24x2   v0, v1, v2, v3, v4, v5, v6, v7;
    ae_int64   t64;
    int N,n,scl;
    NASSERT_ALIGN8(x);
    NASSERT_ALIGN8(y);
    NASSERT(scalingOpt == 3);
    N=ptwd->N;

    pxrd_pos0 = (const ae_int32x2 *)(x        );
    pxrd_neg0 = (const ae_int64   *)(x+N/2-2  );
    pxrd_pos1 = (const ae_int32x2 *)(x+N      );
    pxrd_neg1 = (const ae_int64   *)(x+N+N/2-2);
    pxwr_pos  = (      ae_f24x2   *)(x        );
    pxwr_neg  = (      ae_f24x2   *)(x+N/2-2  );

    /* Transform 2*N samples into N for use in DCT-IV */
    __Pragma("loop_count min=2");
    for(n=0;n<(N>>3);n++)
    {
        t7 = AE_L32X2_X(pxrd_pos0, N/2*sizeof(int32_t));
        t1 = AE_L32X2_X(pxrd_pos1, N/2*sizeof(int32_t));
        AE_L32X2_IP(t4, pxrd_pos0, 2*sizeof(int32_t));
        AE_L32X2_IP(t2, pxrd_pos1, 2*sizeof(int32_t));

        t64 = AE_L64_X(pxrd_neg0, N/2*sizeof(int32_t));     t5 = AE_MOVINT32X2_FROMINT64(t64);
        t64 = AE_L64_X(pxrd_neg1, N/2*sizeof(int32_t));     t3 = AE_MOVINT32X2_FROMINT64(t64);
        AE_L64_IP(t64, pxrd_neg0, -2*(int)sizeof(int32_t)); t6 = AE_MOVINT32X2_FROMINT64(t64);
        AE_L64_IP(t64, pxrd_neg1, -2*(int)sizeof(int32_t)); t0 = AE_MOVINT32X2_FROMINT64(t64);

        v0 = AE_MOVF24X2_FROMINT32X2(AE_SRAI32R(t0, 9));
        v1 = AE_MOVF24X2_FROMINT32X2(AE_SRAI32R(t1, 9));
        v2 = AE_MOVF24X2_FROMINT32X2(AE_SRAI32R(t2, 9));
        v3 = AE_MOVF24X2_FROMINT32X2(AE_SRAI32R(t3, 9));
        v4 = AE_MOVF24X2_FROMINT32X2(AE_SRAI32R(t4, 9));
        v5 = AE_MOVF24X2_FROMINT32X2(AE_SRAI32R(t5, 9));
        v6 = AE_MOVF24X2_FROMINT32X2(AE_SRAI32R(t6, 9));
        v7 = AE_MOVF24X2_FROMINT32X2(AE_SRAI32R(t7, 9));

        v4 = AE_SUB24S(v4, v5);
        v6 = AE_SUB24S(v6, v7);
        v0 = AE_ADD24S(v0, v1);
        v2 = AE_ADD24S(v2, v3);
        v0 = AE_NEG24S(v0);
        v2 = AE_NEG24S(v2);
        v6 = AE_SELP24_LH(v6, v6);
        v2 = AE_SELP24_LH(v2, v2);
        /* make N/2...N-1 samples */
        AE_S32X2F24_X(v4, pxwr_pos, N/2*sizeof(int32_t));
        AE_S32X2F24_X(v6, pxwr_neg, N/2*sizeof(int32_t));
        /* make 0...N/2-1 samples */
        AE_S32X2F24_IP(v0, pxwr_pos, 2*sizeof(int32_t));
        AE_S32X2F24_XP(v2, pxwr_neg, -2*(int)sizeof(int32_t));
    }

    scl=dct4_24x24(y, x, h, scalingOpt);

    return scl+1;

}/* mdct_24x24() */
