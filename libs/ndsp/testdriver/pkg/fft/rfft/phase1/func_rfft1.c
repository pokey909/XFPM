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
#include "../common/test_rfft.h"

static const void* GetFFT_handle(int N, int FFT_Id)
{
    int i, j;
    /* 32-bit handlers for real FFT */
    const  N_handle_pair_t  rfft32x32_hndls[] =
    {
        { 32  , rfft32_32    }, { 64  , rfft32_64    }, { 128 , rfft32_128   }, { 256 , rfft32_256   },
        { 512 , rfft32_512   }, { 1024, rfft32_1024  }, { 2048, rfft32_2048  }, { 4096, rfft32_4096  },
        { 8192, rfft32_8192  }
    };
    const  N_handle_pair_t  rifft32x32_hndls[] =
    {
        { 32  , rifft32_32    }, { 64  , rifft32_64    }, { 128 , rifft32_128   }, { 256 , rifft32_256   },
        { 512 , rifft32_512   }, { 1024, rifft32_1024  }, { 2048, rifft32_2048  }, { 4096, rifft32_4096  },
        { 8192, rifft32_8192  }
    };

    /* 16-bit handlers for real FFT */
    const  N_handle_pair_t  rfft16x16_hndls[] =
    {
        { 32  , rfft16_32    }, { 64  , rfft16_64    }, { 128 , rfft16_128   }, { 256 , rfft16_256   },
        { 512 , rfft16_512   }, { 1024, rfft16_1024  }, { 2048, rfft16_2048  }, { 4096, rfft16_4096  },
        { 8192, rfft16_8192  }
    };
    const  N_handle_pair_t  rifft16x16_hndls[] =
    {
        { 32  , rifft16_32    }, { 64  , rifft16_64    }, { 128 , rifft16_128   }, { 256 , rifft16_256   },
        { 512 , rifft16_512   }, { 1024, rifft16_1024  }, { 2048, rifft16_2048  }, { 4096, rifft16_4096  },
        { 8192, rifft16_8192  }
    };
    /* 32x16-bit handlers for real FFT */
    const  N_handle_pair_t  rfft32x16_hndls[] =
    {
        { 32  , rfft16_32    }, { 64  , rfft16_64    }, { 128 , rfft16_128   }, { 256 , rfft16_256   },
        { 512 , rfft16_512   }, { 1024, rfft16_1024  }, { 2048, rfft16_2048  }, { 4096, rfft16_4096  },
        { 8192, rfft16_8192  }
    };
    const  N_handle_pair_t  rifft32x16_hndls[] =
    {
        { 32  , rifft16_32    }, { 64  , rifft16_64    }, { 128 , rifft16_128   }, { 256 , rifft16_256   },
        { 512 , rifft16_512   }, { 1024, rifft16_1024  }, { 2048, rifft16_2048  }, { 4096, rifft16_4096  },
        { 8192, rifft16_8192  }
    };
    /* 24-bit handlers for real FFT */
    const  N_handle_pair_t  rifft24x24_hndls[] =
    {
        { 32  , rifft24_32   },
        { 64  , rifft24_64   },
        { 128 , rifft24_128  },
        { 256 , rifft24_256  },
        { 512 , rifft24_512  },
        { 1024, rifft24_1024 },
        { 2048, rifft24_2048 },
        { 4096, rifft24_4096 },
        { 8192, rifft24_8192 }
    };
    const  N_handle_pair_t  rfft24x24_hndls[] =
    {
        { 32  , rfft24_32   },
        { 64  , rfft24_64   },
        { 128 , rfft24_128  },
        { 256 , rfft24_256  },
        { 512 , rfft24_512  },
        { 1024, rfft24_1024 },
        { 2048, rfft24_2048 },
        { 4096, rfft24_4096 },
        { 8192, rfft24_8192 }
    };
#define SZ(arr) (int)(sizeof(arr)/sizeof(arr[0]))
    const FFT_handle_tab_t h_tab[] =
    {
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_REAL, FMT_FRACT32), SZ(rfft32x32_hndls), rfft32x32_hndls
        }, 
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_REAL, FMT_FRACT32), SZ(rifft32x32_hndls), rifft32x32_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_REAL, FMT_FRACT16), SZ(rfft16x16_hndls), rfft16x16_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_REAL, FMT_FRACT16), SZ(rifft16x16_hndls), rifft16x16_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_REAL | TE_FFT_32X16, FMT_FRACT32), SZ(rfft32x16_hndls), rfft32x16_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_REAL | TE_FFT_32X16, FMT_FRACT32), SZ(rifft32x16_hndls), rifft32x16_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_REAL | TE_FFT_UNPACKED24, FMT_FRACT32), SZ(rfft24x24_hndls), rfft24x24_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_REAL | TE_FFT_UNPACKED24, FMT_FRACT32), SZ(rifft24x24_hndls), rifft24x24_hndls             
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

static TestDef_t testTbl_rfft[] =
{
#if 1
    { 1, "rfft1/fft32_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft64_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft128_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft256_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft512_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft1024_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft2048_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft4096_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft8192_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },

    { 1, "rfft1/fft32_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft64_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft128_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft256_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft512_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft1024_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft2048_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft4096_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft8192_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL,GetFFT_handle) },

    { 1, "rfft1/ifft32_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft64_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft128_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft256_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft512_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft1024_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft2048_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft4096_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft8192_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },

    { 1, "rfft1/ifft32_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft64_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft128_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft256_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft512_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft1024_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft2048_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft4096_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
    { 1, "rfft1/ifft8192_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16,GetFFT_handle) },
#endif

    { 1, "rfft1/fft32_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft64_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft128_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft256_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft512_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft1024_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft2048_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft4096_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft8192_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft32_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft64_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft128_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft256_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft512_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft1024_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft2048_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft4096_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft8192_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft32_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft64_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft128_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft256_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft512_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft1024_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft2048_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft4096_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft8192_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft32_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft64_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft128_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft256_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft512_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft1024_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft2048_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft4096_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/fft8192_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft32_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft64_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft128_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft256_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft512_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft1024_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft2048_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft4096_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft8192_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft32_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft64_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft128_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft256_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft512_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft1024_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft2048_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft4096_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft8192_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft32_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft64_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft128_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft256_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft512_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft1024_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft2048_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft4096_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft8192_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft32_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft64_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft128_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft256_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft512_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft1024_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft2048_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft4096_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },
    { 1, "rfft1/ifft8192_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24,GetFFT_handle) },

    { 1, "rfft1/fft32_real32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft64_real32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft128_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft256_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft512_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft1024_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft2048_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft4096_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft8192_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
#if 0 // not available in the HiFi3/3z APIs
    { 1, "rfft1/fft32_real32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft64_real32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft128_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft256_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft512_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft1024_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft2048_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft4096_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
    { 1, "rfft1/fft8192_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL,GetFFT_handle) },
#endif
    { 1, "rfft1/ifft32_real32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft64_real32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft128_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft256_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft512_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft1024_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft2048_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft4096_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft8192_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
#if 0 // not available in the HiFi3/3z APIs
    { 1, "rfft1/ifft32_real32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft64_real32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft128_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft256_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft512_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft1024_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft2048_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft4096_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
    { 1, "rfft1/ifft8192_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16,GetFFT_handle) },
#endif

    { 1, "rfft1/fft32_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft64_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft128_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft256_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft512_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft1024_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft2048_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft4096_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft8192_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },

    { 1, "rfft1/fft32_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft64_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft128_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft256_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft512_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft1024_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft2048_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft4096_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },
    { 1, "rfft1/fft8192_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL,GetFFT_handle) },

    { 1, "rfft1/ifft32_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rfft1/ifft64_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rfft1/ifft128_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rfft1/ifft256_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rfft1/ifft512_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rfft1/ifft1024_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rfft1/ifft2048_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rfft1/ifft4096_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },

    { 1, "rfft1/ifft32_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rfft1/ifft64_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rfft1/ifft128_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rfft1/ifft256_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rfft1/ifft512_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rfft1/ifft1024_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rfft1/ifft2048_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rfft1/ifft4096_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 1, "rfft1/ifft8192_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32,GetFFT_handle) },
    { 0 } /* End of table */
}; //testTbl_rfft

/* Perform all tests for fft_realMxN, ifft_realMxN API functions. */
int main_rfft( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
  int ix, res;
  for ( ix=0,res=1; testTbl_rfft[ix].seqFilename && ( res || !breakOnError ); ix++ )
  {
    if ( phaseNum == 0 || phaseNum == testTbl_rfft[ix].phaseNum )
    {
        tTestEngTarget target = testTbl_rfft[ix].target;
        /* Make sure that all functions is present */
        if (!IS_PRESENT(testTbl_rfft[ix].desc.frwTransFxn) &&
            !IS_PRESENT(testTbl_rfft[ix].desc.invTransFxn))
        {
            target = (tTestEngTarget)testTbl_rfft[ix].desc.frwTransFxn;
        }

        res &= ( 0 != TestEngRun(  target,
                                   &testTbl_rfft[ix].desc.desc,
                                   testTbl_rfft[ix].seqFilename,
                                   isFull, isVerbose, breakOnError,0 ) );
    }
  }
  return (res);
} /* main_rfft() */

int func_rfft1(int isFull, int isVerbose, int breakOnError)
{
    return main_rfft(1, isFull, isVerbose, breakOnError);
}
