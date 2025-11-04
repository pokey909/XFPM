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

#include <math.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* Library API */
#include LIBRARY_HEADER(complex)
/* Test engine API. */
#include "testeng.h"

/* Test executive function. Performs the specified test on a brief or full or sanity version
 * of the designated SEQ-file. Return the test result (non-zero if passed). */
static int testExec( tTestEngTarget targetFxn, const char * seqName,
                     int errhExtendedTest, int isFull, int isVerbose, int breakOnError, int testBitexactness );

#define DO_TEST(fxn, seqFile, extraFlags)                                                \
{                                                                                                     \
    if ( res || !breakOnError ) res &= ( 0 != testExec( (tTestEngTarget)(fxn), "complexv1/" seqFile,       \
                                                        (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE, \
                                         isFull, isVerbose, breakOnError,0 ) ); }
int func_complexv1(int isFull, int isVerbose, int breakOnError)
{
  int res = 1;
  DO_TEST( &vec_cplxconj32x32   , "vec_cplxconj32x32.seq", 1 );
  DO_TEST( &vec_cplxconj16x16   , "vec_cplxconj16x16.seq",  1 );
  DO_TEST( &vec_cplx2cplx_mult32x32 , "vec_cplx2cplx_mult32x32.seq" ,1);
  DO_TEST( &vec_cplx2cplx_mult16x16 , "vec_cplx2cplx_mult16x16.seq", 1);
  DO_TEST( &vec_cplx2real_multv32x32, "vec_cplx2real_multv32x32.seq",1);
  DO_TEST( &vec_cplx2real_mults32x32, "vec_cplx2real_mults32x32.seq",1);
  DO_TEST( &vec_cplx2real_multv16x16, "vec_cplx2real_multv16x16.seq",1);
  DO_TEST( &vec_cplx2real_mults16x16, "vec_cplx2real_mults16x16.seq",1);
  DO_TEST( &vec_complex2mag32x32 , "vec_complex2mag32x32.seq",1 );
  DO_TEST( &vec_complex2mag16x16, "vec_complex2mag16x16.seq",1 );
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
  #define TEST_DESC( fmt,extraParam,dimNum, align, loadFxn, procFxn ) { (fmt),0,NULL,(dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }

  /* vec API test definitions. */
  static const struct 
  {
    tTestEngTarget   funcList[MAX_FUNC_NUM];
    tTestEngDesc     testDesc;
  }
  testDefTbl[] =
  {
    {
	  FUNC_LIST( (tTestEngTarget)&vec_cplxconj32x32),
	  TEST_DESC( FMT_CPLX|FMT_FRACT32, 0, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvZc, &te_processFxn_vZvX ) },
	{
	  FUNC_LIST( (tTestEngTarget)&vec_cplxconj16x16 ),
	  TEST_DESC( FMT_CPLX|FMT_FRACT16, 0, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvZc, &te_processFxn_vZvX ) },
	{
	  FUNC_LIST( (tTestEngTarget)&vec_cplx2cplx_mult32x32 ),
	  TEST_DESC( FMT_CPLX|FMT_FRACT32,0, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZc, &te_processFxn_vZvXvY ) },
	{
	  FUNC_LIST( (tTestEngTarget)&vec_cplx2cplx_mult16x16 ),
	  TEST_DESC( FMT_CPLX|FMT_FRACT16,0, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZc, &te_processFxn_vZvXvY ) },
	{
	  FUNC_LIST( (tTestEngTarget)&vec_cplx2real_multv32x32 ),
	  TEST_DESC( FMT_CPLX|FMT_FRACT32,0, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXcvYrvZc, &te_processFxn_vZvXvY ) },
	{
	  FUNC_LIST( (tTestEngTarget)&vec_cplx2real_mults32x32 ),
	  TEST_DESC( FMT_CPLX|FMT_FRACT32,0, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXcsYrvZc, &te_processFxn_vZcvXcsYr ) },
	{
	  FUNC_LIST( (tTestEngTarget)&vec_cplx2real_multv16x16 ),
	  TEST_DESC( FMT_CPLX|FMT_FRACT16,0, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXcvYrvZc, &te_processFxn_vZvXvY ) },
	{
	  FUNC_LIST( (tTestEngTarget)&vec_cplx2real_mults16x16 ),
	  TEST_DESC( FMT_CPLX|FMT_FRACT16,0, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXcsYrvZc, &te_processFxn_vZcvXcsYr ) },
	{
	  FUNC_LIST( (tTestEngTarget)&vec_cplx2cplx_mult16x16 ),
	  TEST_DESC( FMT_CPLX|FMT_FRACT16,0, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZc, &te_processFxn_vZvXvY ) },
	{
	  FUNC_LIST( (tTestEngTarget)&vec_complex2mag16x16),
	  TEST_DESC( FMT_REAL|FMT_FRACT16, 0, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXcvZ, &te_processFxn_vZvX ) },
	{
	  FUNC_LIST( (tTestEngTarget)&vec_complex2mag32x32 ),
	  TEST_DESC( FMT_REAL|FMT_FRACT32, 0, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXcvZ, &te_processFxn_vZvX ) },
    { 
      FUNC_LIST( NULL ), TEST_DESC(  0, 0,0, 0, NULL, NULL ) } /* End of table */
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
