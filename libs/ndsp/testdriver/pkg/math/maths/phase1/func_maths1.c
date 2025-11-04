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
 * Test procedures for vector mathematics
 */

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environemnt convfiguration. */
#include "config.h"
#include "packages.h"
/* DSP Library API: arithmetic and logic functions on data vectors. */
#include LIBRARY_HEADER(math)
/* Test engine API. */
#include "../../common/testeng_math.h"

/* Test executive function. Performs the specified test on a brief or full or sanity version
 * of the designated SEQ-file. Return the test result (non-zero if passed). */
static int testExec( tTestEngTarget targetFxn, const char * seqName,
                     int errhExtendedTest, int isFull, int isVerbose, int breakOnError, int testBitexactness );


/* Perform all tests for math API functions. */
int func_maths1(int isFull, int isVerbose, int breakOnError)
{
  int res = 1;

  #define DO_TEST( fxn, seqFile, extraFlags )                                                         \
    if ( res || !breakOnError ) res &= ( 0 != testExec( (tTestEngTarget)(fxn), "maths1/" seqFile,       \
                                                        (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE, \
                                                        isFull, isVerbose, breakOnError, 0 ) )
  
  DO_TEST( &scl_recip16x16         , "scl_recip16x16.seq"       , 0 );
  DO_TEST( &scl_recip24x24         , "scl_recip24x24.seq"       , 0 );
  DO_TEST( &scl_recip32x32         , "scl_recip32x32.seq"       , 0 );
  DO_TEST( &scl_divide16x16        , "scl_divide16x16.seq"      , 0 );
  DO_TEST( &scl_divide32x32        , "scl_divide32x32.seq"      , 0 );
  DO_TEST( &scl_divide24x24        , "scl_divide24x24.seq"      , 0 );
  DO_TEST( &scl_divide64x32        , "scl_divide64x32.seq"      , 0 );
  DO_TEST( &scl_log2_32x32         , "scl_log2_32x32.seq"       , 0 );
  DO_TEST( &scl_logn_32x32         , "scl_logn_32x32.seq"       , 0 );
  DO_TEST( &scl_log10_32x32        , "scl_log10_32x32.seq"      , 0 );
  DO_TEST( &scl_log2_24x24         , "scl_log2_24x24.seq"       , 0 );
  DO_TEST( &scl_logn_24x24         , "scl_logn_24x24.seq"       , 0 );
  DO_TEST( &scl_log10_24x24        , "scl_log10_24x24.seq"      , 0 );
  DO_TEST( &scl_antilog2_32x32     , "scl_antilog2_32x32.seq"   , 0 );
  DO_TEST( &scl_antilogn_32x32     , "scl_antilogn_32x32.seq"   , 0 );
  DO_TEST( &scl_antilog10_32x32    , "scl_antilog10_32x32.seq"  , 0 );
  DO_TEST( &scl_antilog2_24x24     , "scl_antilog2_24x24.seq"   , 0 );
  DO_TEST( &scl_antilogn_24x24     , "scl_antilogn_24x24.seq"   , 0 );
  DO_TEST( &scl_antilog10_24x24    , "scl_antilog10_24x24.seq"  , 0 );
  DO_TEST( &scl_sqrt16x16          , "scl_sqrt16x16.seq"        , 0 );
  DO_TEST( &scl_sqrt24x24          , "scl_sqrt24x24.seq"        , 0 );
  DO_TEST( &scl_sqrt32x32          , "scl_sqrt32x32.seq"        , 0 );
  DO_TEST( &scl_sqrt64x32          , "scl_sqrt64x32.seq"        , 0 );
  DO_TEST( &scl_sine24x24          , "scl_sine24x24.seq"        , 0 );
  DO_TEST( &scl_sine32x32          , "scl_sine32x32.seq"        , 0 );
  DO_TEST( &scl_cosine24x24        , "scl_cosine24x24.seq"      , 0 );
  DO_TEST( &scl_cosine32x32        , "scl_cosine32x32.seq"      , 0 );
  DO_TEST( &scl_atan24x24          , "scl_atan24x24.seq"        , 0 );
  DO_TEST( &scl_atan32x32          , "scl_atan32x32.seq"        , 0 );
  DO_TEST( &scl_atan2_24x24        , "scl_atan2_24x24.seq"      , 0 );
  DO_TEST( &scl_tan24x24           , "scl_tan24x24.seq"         , 0 );
  DO_TEST( &scl_tan32x32           , "scl_tan32x32.seq"         , 0 );
  DO_TEST( &scl_rsqrt16x16         , "scl_rsqrt16x16.seq"       , 0 );
  DO_TEST( &scl_rsqrt32x32         , "scl_rsqrt32x32.seq"       , 0 );
  DO_TEST( &scl_tanh32x32          , "scl_tanh32x32.seq"        , 0 );
  DO_TEST( &scl_sigmoid32x32       , "scl_sigmoid32x32.seq"     , 0 );

  return (res);
}

/* Test executive function. Performs the specified test on a brief or full version
 * of the designated SEQ-file. Return the test result (non-zero if passed). */
static int testExec( tTestEngTarget   targetFxn, const char * seqName, 
              int errhExtendedTest, int isFull, int isVerbose, int breakOnError, int testBitexactness )
{
  #define MAX_FUNC_NUM   16
  /* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
  #define FUNC_LIST(...) { __VA_ARGS__, NULL }
  /* Initializer for a test description structure. */
  #define TEST_DESC( fmt, dimNum, align, loadFxn, procFxn ) { (fmt),0,NULL,(dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }

  /* vec API test definitions. */
  static const struct 
  {
    tTestEngTarget   funcList[MAX_FUNC_NUM];
    tTestEngDesc     testDesc;
  }
  testDefTbl[] =
  {
    /*
     * Stage 1
     */
  #if 1// HiFi3/3z API
    {
      FUNC_LIST( (tTestEngTarget)&scl_atan24x24),
      TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_math_processFxn_scl_vXvZ) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_sqrt24x24),
      TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_math_processFxn_scl_vXvZ) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_sine24x24,     (tTestEngTarget)&scl_cosine24x24,   (tTestEngTarget)&scl_tan24x24,
                 (tTestEngTarget)&scl_log2_24x24,    (tTestEngTarget)&scl_logn_24x24,    (tTestEngTarget)&scl_log10_24x24,
                 (tTestEngTarget)&scl_antilog2_24x24,(tTestEngTarget)&scl_antilogn_24x24,(tTestEngTarget)&scl_antilog10_24x24),
                 TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_math_processFxn_scl_vXvZ) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_recip24x24),
      TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZd, &te_math_processFxn_scl_recip32x32) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_divide24x24),
      TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZd, &te_math_processFxn_scl_divide32x32) },
    {
      FUNC_LIST((tTestEngTarget)&scl_atan2_24x24),
      TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_math_processFxn_scl_atan2) },
  #endif
    {
      FUNC_LIST( (tTestEngTarget)&scl_sine32x32,     (tTestEngTarget)&scl_cosine32x32,   (tTestEngTarget)&scl_tan32x32,
                 (tTestEngTarget)&scl_log2_32x32,    (tTestEngTarget)&scl_logn_32x32,    (tTestEngTarget)&scl_log10_32x32,
                 (tTestEngTarget)&scl_antilog2_32x32,(tTestEngTarget)&scl_antilogn_32x32,(tTestEngTarget)&scl_antilog10_32x32),
      TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_math_processFxn_scl_vXvZ ) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_atan32x32),
      TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_math_processFxn_scl_vXvZ) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_sqrt32x32),
      TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_math_processFxn_scl_vXvZ) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_recip16x16),
      TEST_DESC(FMT_REAL | FMT_FRACT16, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZf, &te_math_processFxn_scl_recip16x16) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_recip32x32),
      TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZd, &te_math_processFxn_scl_recip32x32) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_divide16x16),
      TEST_DESC(FMT_REAL | FMT_FRACT16, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZf, &te_math_processFxn_scl_divide16x16) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_divide32x32),
      TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZd, &te_math_processFxn_scl_divide32x32) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_divide64x32),
      TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vX64vYvZ, &te_math_processFxn_scl_divide64x32 ) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_rsqrt16x16),
      TEST_DESC(FMT_REAL | FMT_FRACT16, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZf, &te_math_processFxn_scl_recip16x16) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_rsqrt32x32),
      TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZd, &te_math_processFxn_scl_recip32x32) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_tanh32x32, (tTestEngTarget)&scl_sigmoid32x32),
      TEST_DESC( FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_math_processFxn_scl_vXvZ ) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_sqrt16x16),
      TEST_DESC( FMT_FRACT16, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_math_processFxn_scl_vXvZ ) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_sqrt64x32),
      TEST_DESC( FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vX64vZ, &te_math_processFxn_scl_vXcvZ ) },
    { 
      FUNC_LIST( NULL ), TEST_DESC(  0, 0, 0, NULL, NULL ) } /* End of table */
  };

  {
    int tblIx, funcIx;

    for ( tblIx=0; tblIx<(int)(sizeof(testDefTbl)/sizeof(testDefTbl[0])); tblIx++ )
    {
      for ( funcIx=0; funcIx<MAX_FUNC_NUM; funcIx++ )
      {
        if ( targetFxn == testDefTbl[tblIx].funcList[funcIx] )
        {
          tTestEngDesc testDesc = testDefTbl[tblIx].testDesc;
          testDesc.extraParam = (uint32_t)errhExtendedTest;

          return ( TestEngRun( targetFxn, &testDesc, 
                               seqName, isFull, 
                               isVerbose, breakOnError, testBitexactness ) );
        }
      }
    }

    ASSERT( !"Test not defined" );
    return (0);
  }
  return te_Exec(testDefTbl, sizeof(testDefTbl) / sizeof(testDefTbl[0]), MAX_FUNC_NUM, targetFxn, seqName, isFull, isVerbose, breakOnError);
}
