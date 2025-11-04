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
#include "../../rfft/common/test_rfft.h"

static const void* GetFFT_handle(int N, int FFT_Id)
{
    int i, j;
    /* 32-bit handlers for real FFT */
    const  N_handle_pair_t  rfft32x32_hndls[] =
    {
        { 12  , rnfft32_12   }, { 24  , rnfft32_24   }, { 36  , rnfft32_36   },
        { 48  , rnfft32_48   }, { 60  , rnfft32_60   }, { 72  , rnfft32_72   }, { 96  , rnfft32_96   },
        { 108 , rnfft32_108  }, { 120 , rnfft32_120  }, { 144 , rnfft32_144  }, { 180 , rnfft32_180  },
        { 192 , rnfft32_192  }, { 216 , rnfft32_216  }, { 240 , rnfft32_240  }, { 288 , rnfft32_288  },
        { 300 , rnfft32_300  }, { 324 , rnfft32_324  }, { 360 , rnfft32_360  }, { 432 , rnfft32_432  },
        { 480 , rnfft32_480  }, { 540 , rnfft32_540  }, { 576 , rnfft32_576  }, { 768 , rnfft32_768  },
        { 960 , rnfft32_960  }, { 30  , rnfft32_30   }, { 90  , rnfft32_90   }, { 384 , rnfft32_384  },
        { 720 , rnfft32_720  }, { 1152, rnfft32_1152 }, { 1440, rnfft32_1440 }, { 1536, rnfft32_1536 }, 
        { 1920, rnfft32_1920 }
    };
    const  N_handle_pair_t  rifft32x32_hndls[] =
    {
        { 12  , rinfft32_12   }, { 24  , rinfft32_24   }, { 36  , rinfft32_36   },
        { 48  , rinfft32_48   }, { 60  , rinfft32_60   }, { 72  , rinfft32_72   }, { 96  , rinfft32_96   },
        { 108 , rinfft32_108  }, { 120 , rinfft32_120  }, { 144 , rinfft32_144  }, { 180 , rinfft32_180  },
        { 192 , rinfft32_192  }, { 216 , rinfft32_216  }, { 240 , rinfft32_240  }, { 288 , rinfft32_288  },
        { 300 , rinfft32_300  }, { 324 , rinfft32_324  }, { 360 , rinfft32_360  }, { 432 , rinfft32_432  },
        { 480 , rinfft32_480  }, { 540 , rinfft32_540  }, { 576 , rinfft32_576  }, { 768 , rinfft32_768  },
        { 960 , rinfft32_960  }, { 30  , rinfft32_30   }, { 90  , rinfft32_90   }, { 384 , rinfft32_384  },
        { 720 , rinfft32_720  }, { 1152, rinfft32_1152 }, { 1440, rinfft32_1440 }, { 1536, rinfft32_1536 },
        { 1920, rinfft32_1920 }
    };
#define SZ(arr) (int)(sizeof(arr)/sizeof(arr[0]))
    const FFT_handle_tab_t h_tab[] =
    {
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_REAL, FMT_FRACT32), SZ(rfft32x32_hndls), rfft32x32_hndls
        }, 
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_REAL, FMT_FRACT32), SZ(rifft32x32_hndls), rifft32x32_hndls
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
static TestDef_t testTbl_rnfft[] =
{
    // Mixed radix rfft
    { 1, "rnfft1/fft12_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft24_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft30_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft36_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft48_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft60_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft72_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft90_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft96_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft108_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft120_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft144_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft180_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft192_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft216_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft240_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft288_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft300_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft324_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft360_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft384_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft432_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft480_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft540_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft576_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft720_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft768_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft960_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft1152_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft1440_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft1536_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft1920_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },

    { 1, "rnfft1/fft12_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft24_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft30_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft36_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft48_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft60_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft72_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft90_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft96_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft108_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft120_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft144_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft180_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft192_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft216_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft240_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft288_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft300_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft324_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft360_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft384_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft432_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft480_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft540_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft576_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft720_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft768_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft960_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft1152_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft1440_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft1536_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rnfft1/fft1920_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },

    { 1, "rnfft1/ifft12_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft24_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft30_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft36_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft48_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft60_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft72_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft90_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft96_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft108_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft120_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft144_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft180_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft192_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft216_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft240_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft288_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft300_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft324_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft360_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft384_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft432_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft480_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft540_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft576_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft720_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft768_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft960_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft1152_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft1440_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft1536_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft1920_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },

    { 1, "rnfft1/ifft12_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft24_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft30_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft36_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft48_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft60_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft72_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft90_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft96_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft108_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft120_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft144_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft180_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft192_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft216_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft240_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft288_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft300_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft324_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft360_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft384_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft432_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft480_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft540_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft576_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft720_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft768_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft960_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft1152_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft1440_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft1536_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rnfft1/ifft1920_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 0 } /* End of table */
}; //testTbl_rnfft

/* Perform all tests for fft_realMxN, ifft_realMxN API functions. */
int main_rnfft( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
  int ix, res;
  for ( ix=0,res=1; testTbl_rnfft[ix].seqFilename && ( res || !breakOnError ); ix++ )
  {
    if ( phaseNum == 0 || phaseNum == testTbl_rnfft[ix].phaseNum )
    {
        tTestEngTarget target = testTbl_rnfft[ix].target;
        /* Make sure that all functions is present */
        if (!IS_PRESENT(testTbl_rnfft[ix].desc.frwTransFxn) &&
            !IS_PRESENT(testTbl_rnfft[ix].desc.invTransFxn))
        {
            target = (tTestEngTarget)testTbl_rnfft[ix].desc.frwTransFxn;
        }

        res &= ( 0 != TestEngRun(  target,
                                   &testTbl_rnfft[ix].desc.desc,
                                   testTbl_rnfft[ix].seqFilename,
                                   isFull, isVerbose, breakOnError,0 ) );
    }
  }
  return (res);
} /* main_rnfft() */

int func_rnfft1(int isFull, int isVerbose, int breakOnError)
{
    return main_rnfft(1, isFull, isVerbose, breakOnError);
}
