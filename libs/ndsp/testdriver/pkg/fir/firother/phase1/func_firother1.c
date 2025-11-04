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
 * Test procedures for FIR
 */
#include "../common/test_firother.h"

/* API test definitions. */
static const struct 
{
  tTestEngTarget   funcList[MAX_FUNC_NUM];
  tTestEngDesc     testDesc;
}
testDefTbl[] =
{
    {  FUNC_LIST( (tTestEngTarget)&fir_convol32x32),  TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT32, 0                 , TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_convol16x16),    TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT16, 0                 , TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_convol32x16),    TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT16, TE_FIR_OTHER_32X16, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_convol32x16),  TEST_DESC_CONVOLVE( FMT_CPLX|FMT_FRACT16, TE_FIR_OTHER_32X16, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_convola16x16),   TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT16, 0                 , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_convola32x32),   TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT32, 0                 , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_convola32x16),   TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT16, TE_FIR_OTHER_32X16, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_convola32x16), TEST_DESC_CONVOLVE( FMT_CPLX|FMT_FRACT16, TE_FIR_OTHER_32X16, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_lconvola16x16),  TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT16, TE_FIR_LINEAR_API , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_lconvola32x32),  TEST_DESC_CONVOLVE( FMT_REAL|FMT_FRACT32, TE_FIR_LINEAR_API , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorr32x32),   TEST_DESC_AUTOCORR( FMT_REAL|FMT_FRACT32, 0                 , TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorr16x16),     TEST_DESC_AUTOCORR( FMT_REAL|FMT_FRACT16, 0                 , TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorra16x16),    TEST_DESC_AUTOCORR( FMT_REAL|FMT_FRACT16, 0                 , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorra32x32),    TEST_DESC_AUTOCORR( FMT_REAL|FMT_FRACT32, 0                 , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_lacorra16x16),   TEST_DESC_AUTOCORR( FMT_REAL|FMT_FRACT16, TE_FIR_LINEAR_API , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_lacorra32x32),   TEST_DESC_AUTOCORR( FMT_REAL|FMT_FRACT32, TE_FIR_LINEAR_API , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorr32x32),   TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT32, 0                 , TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorr16x16),     TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT16, 0                 , TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorr32x16),     TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT16, TE_FIR_OTHER_32X16, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorra16x16),    TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT16, 0                 , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorra32x32),    TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT32, 0                 , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorra32x16),    TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT16, TE_FIR_OTHER_32X16, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_lxcorra16x16),   TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT16, TE_FIR_LINEAR_API , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_lxcorra32x32),   TEST_DESC_CROSSCORR( FMT_REAL|FMT_FRACT32, TE_FIR_LINEAR_API , TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms32x32),    { FMT_REAL|FMT_FRACT32, 0                 ,NULL,TE_DIM_NUM_2, TE_ALIGN_YES, &te_create_lms,&te_destroy_lms, &te_loadFxn_lms, &te_processFxn_lms } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms16x16),      { FMT_REAL|FMT_FRACT16, 0                 ,NULL,TE_DIM_NUM_2, TE_ALIGN_YES, &te_create_lms,&te_destroy_lms, &te_loadFxn_lms, &te_processFxn_lms } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms16x32),      { FMT_REAL|FMT_FRACT16, TE_FIR_OTHER_32X16,NULL,TE_DIM_NUM_2, TE_ALIGN_YES, &te_create_lms,&te_destroy_lms, &te_loadFxn_lms, &te_processFxn_lms } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms16x16_convergence),      { FMT_REAL|FMT_FRACT16, 0                 ,(void *)fir_blms16x16  ,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms16x32_convergence),      { FMT_REAL|FMT_FRACT16, TE_FIR_OTHER_32X16,(void *)fir_blms16x32  ,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms32x32_convergence),      { FMT_REAL|FMT_FRACT32, 0,                 (void *)fir_blms32x32  ,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( NULL ),              TEST_DESC(  0, 0, 0, 0, NULL, NULL ) } /* End of table */
};

/* Perform all tests for FIR API functions. */
int func_firother1(int isFull, int isVerbose, int breakOnError)
{
    int res = 1;

    #define DO_TEST(fxn, seqFile)                                                                    \
        if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                               \
                                                           sizeof(testDefTbl)/sizeof(testDefTbl[0]), \
                                                           MAX_FUNC_NUM,                             \
                                                           (tTestEngTarget)(fxn), "firother1/" seqFile,	 \
                                                           isFull, isVerbose, breakOnError ) )
    

    DO_TEST( &fir_acorr16x16     ,  "fir_acorr16x16.seq"      );
    DO_TEST( &fir_acorr32x32     ,  "fir_acorr32x32.seq"      );
    DO_TEST( &fir_xcorr16x16     ,  "fir_xcorr16x16.seq"      );
    DO_TEST( &fir_xcorr32x16     ,  "fir_xcorr32x16.seq"      );

    DO_TEST( &fir_xcorr32x32     ,  "fir_xcorr32x32.seq"      );
    DO_TEST( &fir_convol16x16    ,  "fir_convol16x16.seq"     );
    DO_TEST( &fir_convol32x16    ,  "fir_convol32x16.seq"     );

    DO_TEST( &fir_convol32x32    ,  "fir_convol32x32.seq"     );
    DO_TEST( &cxfir_convol32x16  ,  "cxfir_convol32x16.seq"   );
    DO_TEST( &fir_acorra16x16    ,  "fir_acorra16x16.seq"     );

    DO_TEST( &fir_acorra32x32    ,  "fir_acorra32x32.seq"     );
    DO_TEST( &fir_lacorra16x16   ,  "fir_lacorra16x16.seq"    );
    DO_TEST( &fir_lacorra32x32   ,  "fir_lacorra32x32.seq"    );
    DO_TEST( &fir_xcorra16x16    ,  "fir_xcorra16x16.seq"     );
    DO_TEST( &fir_xcorra32x16    ,  "fir_xcorra32x16.seq"     );

    DO_TEST( &fir_xcorra32x32    ,  "fir_xcorra32x32.seq"     );
    DO_TEST( &fir_lxcorra16x16   ,  "fir_lxcorra16x16.seq"    );
    DO_TEST( &fir_lxcorra32x32   ,  "fir_lxcorra32x32.seq"    );
    DO_TEST( &fir_convola16x16   ,  "fir_convola16x16.seq"    );
    DO_TEST( &fir_convola32x16   ,  "fir_convola32x16.seq"    );
    DO_TEST( &fir_convola32x32   ,  "fir_convola32x32.seq"    );
    DO_TEST( &cxfir_convola32x16 ,  "cxfir_convola32x16.seq"  );
    DO_TEST( &fir_lconvola16x16  ,  "fir_lconvola16x16.seq"   );
    DO_TEST( &fir_lconvola32x32  ,  "fir_lconvola32x32.seq"   );

    DO_TEST( &fir_blms16x16      ,  "fir_blms16x16.seq"       );
    DO_TEST( &fir_blms32x32      ,  "fir_blms32x32.seq"       );
    DO_TEST( &fir_blms16x32      ,  "fir_blms16x32.seq"       );
    DO_TEST( &fir_blms16x16_convergence    ,  "fir_blms16x16_convergence.seq"    );
    DO_TEST( &fir_blms16x32_convergence    ,  "fir_blms16x32_convergence.seq"    );
    DO_TEST( &fir_blms32x32_convergence    ,  "fir_blms32x32_convergence.seq"    );

    return (res);
}
