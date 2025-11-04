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
#include "../common/test_cfft.h"

static const void* GetFFT_handle(int N, int FFT_Id)
{
    int i, j;
    /* 32-bit handlers for complex FFT */
    const  N_handle_pair_t  cfft32x32_hndls[] =
    {
        { 32 , cfft32_32   }, { 64  , cfft32_64   }, { 128 , cfft32_128  }, { 256 , cfft32_256  },
        { 512, cfft32_512  }, { 1024, cfft32_1024 }, { 2048, cfft32_2048 }, { 4096, cfft32_4096 },
        { 16 , cfft32_16   }};
    const  N_handle_pair_t  cifft32x32_hndls[] =
    {
        { 32 , cifft32_32   }, { 64  , cifft32_64   }, { 128 , cifft32_128  }, { 256 , cifft32_256  },
        { 512, cifft32_512  }, { 1024, cifft32_1024 }, { 2048, cifft32_2048 }, { 4096, cifft32_4096 },
        { 16 , cifft32_16   }};
    /* 16-bit handlers for complex FFT */
    const  N_handle_pair_t  cfft16x16_hndls[] =
    {
        { 16  , cfft16_16    }, { 32  , cfft16_32    }, { 64  , cfft16_64    }, { 128 , cfft16_128   },
        { 256 , cfft16_256   }, { 512 , cfft16_512   }, { 1024, cfft16_1024  }, { 2048, cfft16_2048  },
        { 4096, cfft16_4096  },
    };
    const  N_handle_pair_t  cifft16x16_hndls[] =
    {
        { 16  , cifft16_16    }, { 32  , cifft16_32    }, { 64  , cifft16_64    }, { 128 , cifft16_128   },
        { 256 , cifft16_256   }, { 512 , cifft16_512   }, { 1024, cifft16_1024  }, { 2048, cifft16_2048  },
        { 4096, cifft16_4096  }
    };
    /* 32x16-bit handlers for complex FFT */
    const  N_handle_pair_t  cfft32x16_hndls[] =
    {
        { 16  , cfft16_16    }, { 32  , cfft16_32    }, { 64  , cfft16_64    }, { 128 , cfft16_128   },
        { 256 , cfft16_256   }, { 512 , cfft16_512   }, { 1024, cfft16_1024  }, { 2048, cfft16_2048  },
        { 4096, cfft16_4096  }
    };
    const  N_handle_pair_t  cifft32x16_hndls[] =
    {
        { 16  , cifft16_16    }, { 32  , cifft16_32    }, { 64  , cifft16_64    }, { 128 , cifft16_128   },
        { 256 , cifft16_256   }, { 512 , cifft16_512   }, { 1024, cifft16_1024  }, { 2048, cifft16_2048  },
        { 4096, cifft16_4096  }
    };
    /* 24-bit handlers for complex FFT */
    const  N_handle_pair_t  cfft24x24_hndls[] =
    {
        { 16  , cfft24_16   },
        { 32  , cfft24_32   },
        { 64  , cfft24_64   },
        { 128 , cfft24_128  },
        { 256 , cfft24_256  },
        { 512 , cfft24_512  },
        { 1024, cfft24_1024 },
        { 2048, cfft24_2048 },
        { 4096, cfft24_4096 }
    };
    const  N_handle_pair_t  cifft24x24_hndls[] =
    {
        { 16  , cifft24_16   },
        { 32  , cifft24_32   },
        { 64  , cifft24_64   },
        { 128 , cifft24_128  },
        { 256 , cifft24_256  },
        { 512 , cifft24_512  },
        { 1024, cifft24_1024 },
        { 2048, cifft24_2048 },
        { 4096, cifft24_4096 }
    };
#define SZ(arr) (int)(sizeof(arr)/sizeof(arr[0]))
    const FFT_handle_tab_t h_tab[] =
    {
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_CPLX, FMT_FRACT32), SZ(cfft32x32_hndls), cfft32x32_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_CPLX, FMT_FRACT32), SZ(cifft32x32_hndls), cifft32x32_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_CPLX, FMT_FRACT16), SZ(cfft16x16_hndls), cfft16x16_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_CPLX, FMT_FRACT16), SZ(cifft16x16_hndls), cifft16x16_hndls
        },        
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_CPLX | TE_FFT_32X16, FMT_FRACT32), SZ(cfft32x16_hndls), cfft32x16_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_CPLX | TE_FFT_32X16, FMT_FRACT32), SZ(cifft32x16_hndls), cifft32x16_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, (TE_FFT_CPLX | TE_FFT_UNPACKED24), FMT_FRACT32), SZ(cfft24x24_hndls), cfft24x24_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, (TE_FFT_CPLX | TE_FFT_UNPACKED24), FMT_FRACT32), SZ(cifft24x24_hndls), cifft24x24_hndls
        },
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

TestDef_t testTbl_cfft[] =
{
#if 1
    { 1, "cfft1/fft16_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft32_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft64_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft128_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft256_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft512_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft1024_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft2048_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft4096_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft16_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft32_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft64_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft128_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft256_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft512_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft1024_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft2048_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/fft4096_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft16_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft32_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft64_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft128_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft256_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft512_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft1024_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft2048_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft4096_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },

    { 1, "cfft1/ifft16_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft32_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft64_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft128_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft256_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft512_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft1024_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft2048_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
    { 1, "cfft1/ifft4096_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16,GetFFT_handle) },
#endif
#if 1 //HiFi3/3z API
    { 1, "cfft1/fft16_cplx24x24s0.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft32_cplx24x24s0.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft64_cplx24x24s0.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft128_cplx24x24s0.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft256_cplx24x24s0.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft512_cplx24x24s0.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft1024_cplx24x24s0.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft2048_cplx24x24s0.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft4096_cplx24x24s0.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft16_cplx24x24s1.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft32_cplx24x24s1.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft64_cplx24x24s1.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft128_cplx24x24s1.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft256_cplx24x24s1.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft512_cplx24x24s1.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft1024_cplx24x24s1.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft2048_cplx24x24s1.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft4096_cplx24x24s1.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft16_cplx24x24s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft32_cplx24x24s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft64_cplx24x24s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft128_cplx24x24s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft256_cplx24x24s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft512_cplx24x24s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft1024_cplx24x24s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft2048_cplx24x24s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft4096_cplx24x24s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft16_cplx24x24s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft32_cplx24x24s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft64_cplx24x24s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft128_cplx24x24s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft256_cplx24x24s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft512_cplx24x24s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft1024_cplx24x24s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft2048_cplx24x24s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/fft4096_cplx24x24s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },

    { 1, "cfft1/ifft16_cplx24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft32_cplx24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft64_cplx24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft128_cplx24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft256_cplx24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft512_cplx24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft1024_cplx24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft2048_cplx24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft4096_cplx24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft16_cplx24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft32_cplx24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft64_cplx24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft128_cplx24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft256_cplx24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft512_cplx24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft1024_cplx24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft2048_cplx24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft4096_cplx24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft16_cplx24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft32_cplx24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft64_cplx24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft128_cplx24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft256_cplx24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft512_cplx24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft1024_cplx24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft2048_cplx24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft4096_cplx24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft16_cplx24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft32_cplx24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft64_cplx24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft128_cplx24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft256_cplx24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft512_cplx24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft1024_cplx24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft2048_cplx24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
    { 1, "cfft1/ifft4096_cplx24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24,GetFFT_handle) },
#endif

#if 1
    { 1, "cfft1/fft16_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/fft32_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/fft64_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/fft128_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/fft256_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/fft512_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/fft1024_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/fft2048_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/fft4096_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/ifft16_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/ifft32_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/ifft64_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/ifft128_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/ifft256_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/ifft512_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/ifft1024_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/ifft2048_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
    { 1, "cfft1/ifft4096_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16,GetFFT_handle) },
#endif

#if 1
    { 1, "cfft1/fft16_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft32_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft64_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft128_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft256_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft512_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft1024_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft2048_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft4096_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },

    { 1, "cfft1/fft16_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft32_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft64_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft128_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft256_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft512_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft1024_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft2048_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },
    { 1, "cfft1/fft4096_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL,GetFFT_handle) },

    { 1, "cfft1/ifft16_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft32_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft64_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft128_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft256_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft512_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft1024_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft2048_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft4096_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },

    { 1, "cfft1/ifft16_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft32_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft64_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft128_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft256_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft512_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft1024_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft2048_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
    { 1, "cfft1/ifft4096_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32,GetFFT_handle) },
#endif
    { 0 } /* End of table */
}; // testTbl_cfft

/* Perform all tests for fft_cplx, ifft_cplx API functions. */
int main_cfft(int phaseNum, int isFull, int isVerbose, int breakOnError)
{
    int ix, res;
    for (ix = 0, res = 1; testTbl_cfft[ix].seqFilename && (res || !breakOnError); ix++)
    {
        if (phaseNum == 0 || phaseNum == testTbl_cfft[ix].phaseNum)
        {
            tTestEngTarget target = testTbl_cfft[ix].target;

            if (!IS_PRESENT(testTbl_cfft[ix].desc.frwTransFxn) &&
                !IS_PRESENT(testTbl_cfft[ix].desc.invTransFxn))
            {
                target = (tTestEngTarget)testTbl_cfft[ix].desc.frwTransFxn;
            }

            res &= (0 != TestEngRun(target,
                &testTbl_cfft[ix].desc.desc,
                testTbl_cfft[ix].seqFilename,
                isFull, isVerbose, breakOnError,0));
        }
    }
    return (res);

} /* main_cfft() */

int func_cfft1(int isFull, int isVerbose, int breakOnError)
{
    return main_cfft(1, isFull, isVerbose, breakOnError);
}
