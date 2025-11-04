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
* Test procedures for complex FFT functions.
*/
#include "../../cfft/common/test_cfft.h"

static const void* GetFFT_handle(int N, int FFT_Id)
{
    int i, j;
    /* 32-bit handlers for complex FFT */
    const  N_handle_pair_t  cfft32x32_hndls[] =
    {
        { 12  , cnfft32_12  }, { 24  , cnfft32_24  }, { 36  , cnfft32_36  },
        { 48 , cnfft32_48  }, { 60  , cnfft32_60  }, { 72  , cnfft32_72  }, { 96  , cnfft32_96  },
        { 108, cnfft32_108 }, { 120 , cnfft32_120 }, { 144 , cnfft32_144 }, { 180 , cnfft32_180 },
        { 192, cnfft32_192 }, { 216 , cnfft32_216 }, { 240 , cnfft32_240 }, { 288 , cnfft32_288 },
        { 300, cnfft32_300 }, { 324 , cnfft32_324 }, { 360 , cnfft32_360 }, { 432 , cnfft32_432 },
        { 480, cnfft32_480 }, { 540 , cnfft32_540 }, { 576 , cnfft32_576 }, { 768 , cnfft32_768 },
        { 960, cnfft32_960 }, { 600 , cnfft32_600 }, { 400 , cnfft32_400 }, { 384 , cnfft32_384 },
        { 200, cnfft32_200 }, { 160 , cnfft32_160 }, { 100 , cnfft32_100 }, { 80  , cnfft32_80  }
    };
    const  N_handle_pair_t  cifft32x32_hndls[] =
    {
        { 12  , cinfft32_12  }, { 24  , cinfft32_24  }, { 36  , cinfft32_36  },
        { 48 , cinfft32_48  }, { 60  , cinfft32_60  }, { 72  , cinfft32_72  }, { 96  , cinfft32_96  },
        { 108, cinfft32_108 }, { 120 , cinfft32_120 }, { 144 , cinfft32_144 }, { 180 , cinfft32_180 },
        { 192, cinfft32_192 }, { 216 , cinfft32_216 }, { 240 , cinfft32_240 }, { 288 , cinfft32_288 },
        { 300, cinfft32_300 }, { 324 , cinfft32_324 }, { 360 , cinfft32_360 }, { 432 , cinfft32_432 },
        { 480, cinfft32_480 }, { 540 , cinfft32_540 }, { 576 , cinfft32_576 }, { 768 , cinfft32_768 },
        { 960, cinfft32_960 }, { 600 , cinfft32_600 }, { 400 , cinfft32_400 }, { 384 , cinfft32_384 },
        { 200, cinfft32_200 }, { 160 , cinfft32_160 }, { 100 , cinfft32_100 }, { 80  , cinfft32_80  }
    };
#define SZ(arr) (int)(sizeof(arr)/sizeof(arr[0]))
    const FFT_handle_tab_t h_tab[] =
    {
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_CPLX, FMT_FRACT32), (int)(sizeof(cfft32x32_hndls)/sizeof(cfft32x32_hndls[0])), cfft32x32_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_CPLX, FMT_FRACT32),  (int)(sizeof(cifft32x32_hndls)/sizeof(cifft32x32_hndls[0])), cifft32x32_hndls
        }
    };

    for (i=0; i<SZ(h_tab); i++)
    {
        if (h_tab[i].FFT_Id == FFT_Id)
        {
            for (j = 0; j<h_tab[i].numN; j++)
            {
                if (h_tab[i].pNh[j].N == N)
                {
                    return h_tab[i].pNh[j].h; 
                }
            }
        }
    }
    
    return NULL;
}/* GetFFT_handle() */

TestDef_t testTbl_cnfft[] =
{
    // Mixed radix forward 32x32
    { 1, "cnfft1/fft12_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft24_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft36_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft48_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft60_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft72_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft80_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft96_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft100_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft108_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft120_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft144_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft160_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft180_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft192_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft200_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft216_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft240_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft288_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft300_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft324_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft360_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft384_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft400_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft432_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft480_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft540_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft576_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft600_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft768_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft960_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft12_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft24_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft36_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft48_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft60_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft72_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft80_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft96_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft100_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft108_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft120_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft144_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft160_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft180_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft192_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft200_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft216_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft240_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft288_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft300_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft324_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft360_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft384_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft400_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft432_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft480_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft540_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft576_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft600_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft768_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cnfft1/fft960_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    // Mixed radix inverse 32x32
    { 1, "cnfft1/ifft12_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft24_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft36_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft48_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft60_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft72_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft80_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft96_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft100_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft108_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft120_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft144_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft160_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft180_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft192_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft200_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft216_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft240_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft288_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft300_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft324_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft360_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft384_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft400_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft432_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft480_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft540_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft576_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft600_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft768_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft960_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft12_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft24_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft36_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft48_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft60_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft72_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft80_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft96_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft100_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft108_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft120_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft144_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft160_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft180_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft192_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft200_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft216_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft240_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft288_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft300_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft324_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft360_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft384_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft400_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft432_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft480_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft540_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft576_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft600_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft768_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cnfft1/ifft960_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },

    { 0 } /* End of table */
}; // testTbl_cnfft

/* Perform all tests for mixed radix fft_cplx32x32, ifft_cplx32x32 API functions. */
int main_cnfft(int phaseNum, int isFull, int isVerbose, int breakOnError)
{
    int ix, res;
    for (ix = 0, res = 1; testTbl_cnfft[ix].seqFilename && (res || !breakOnError); ix++)
    {
        if (phaseNum == 0 || phaseNum == testTbl_cnfft[ix].phaseNum)
        {
            tTestEngTarget target = testTbl_cnfft[ix].target;

            if (!IS_PRESENT(testTbl_cnfft[ix].desc.frwTransFxn) &&
                !IS_PRESENT(testTbl_cnfft[ix].desc.invTransFxn))
            {
                target = (tTestEngTarget)testTbl_cnfft[ix].desc.frwTransFxn;
            }

            res &= (0 != TestEngRun(target,
                &testTbl_cnfft[ix].desc.desc,
                testTbl_cnfft[ix].seqFilename,
                isFull, isVerbose, breakOnError,0));
        }
    }
    return (res);

} /* main_cnfft() */

int func_cnfft1(int isFull, int isVerbose, int breakOnError)
{
    return main_cnfft(1, isFull, isVerbose, breakOnError);
}
