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
    {  FUNC_LIST( (tTestEngTarget)&fir_convolf),        TEST_DESC_CONVOLVE ( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorrf),         TEST_DESC_AUTOCORR ( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorrf),         TEST_DESC_CROSSCORR( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_xcorrf),       TEST_DESC_CROSSCORR( FMT_CPLX|FMT_FLOAT32, 0, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_convolaf),       TEST_DESC_CONVOLVE ( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorraf),        TEST_DESC_AUTOCORR ( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorraf),        TEST_DESC_CROSSCORR( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_xcorraf),      TEST_DESC_CROSSCORR( FMT_CPLX|FMT_FLOAT32, 0, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_blmsf),          { FMT_REAL|FMT_FLOAT32, 0,NULL,TE_DIM_NUM_2, TE_ALIGN_YES, &te_create_lms,&te_destroy_lms, &te_loadFxn_lms, &te_processFxn_lms } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blmsf_convergence),          { FMT_REAL|FMT_FLOAT32, 0,                 (void *)fir_blmsf      ,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( NULL ),              TEST_DESC(  0, 0, 0, 0, NULL, NULL ) } /* End of table */
};
/* Perform all tests for FIR API functions. */
int func_firother2(int isFull, int isVerbose, int breakOnError)
{
    int res = 1;

    #define DO_TEST(fxn, seqFile)                                                                    \
        if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                               \
                                                           sizeof(testDefTbl)/sizeof(testDefTbl[0]), \
                                                           MAX_FUNC_NUM,                             \
                                                           (tTestEngTarget)(fxn), "firother2/" seqFile,    \
                                                           isFull, isVerbose, breakOnError ) )
    
    DO_TEST( &fir_acorrf       , "fir_acorrf.seq"       );
    DO_TEST( &fir_xcorrf       , "fir_xcorrf.seq"       );
    DO_TEST( &cxfir_xcorrf     , "cxfir_xcorrf.seq"     );
    DO_TEST( &fir_convolf      , "fir_convolf.seq"      );
    DO_TEST( &fir_acorraf      , "fir_acorraf.seq"      );
    DO_TEST( &fir_xcorraf      , "fir_xcorraf.seq"      );
    DO_TEST( &cxfir_xcorraf    , "cxfir_xcorraf.seq"    );
    DO_TEST( &fir_convolaf     , "fir_convolaf.seq"     );
    DO_TEST( &fir_blmsf        , "fir_blmsf.seq"        );
    DO_TEST( &fir_blmsf_convergence        , "fir_blmsf_convergence.seq"        );

    return (res);
}
