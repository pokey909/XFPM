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
 * Test procedures for arithmetic and logic functions on data vectors.
 */

#include <string.h> 
/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* Library API. */
#include LIBRARY_HEADER(vector)
/* Test Engine add-on for Vector Mathematics functions */
#include "../common/testeng_vector.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

#define MAX_FUNC_NUM   10
/* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
#define FUNC_LIST(...) { __VA_ARGS__, NULL }
/* Initializer for a test description structure. */
#define TEST_DESC( fmt, dimNum, align, loadFxn, procFxn ) { (fmt),0,NULL,(dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }
#define TEST_DESC_EX( fmt, dimNum, align, loadFxn, procFxn,api ) { (fmt),0,(void*)api,(dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }

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
    FUNC_LIST( (tTestEngTarget)&vec_dot24x24 ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYsZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot24x24_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvYsZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_add24x24 ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_add24x24_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_power24x24 ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsY32sZ64, &te_processFxn_vXsY32sZ64 ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_power24x24_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsY32sZ64, &te_processFxn_vXsY32sZ64 ) },
  { 
    FUNC_LIST((tTestEngTarget)&vec_min24x24,
              (tTestEngTarget)&vec_max24x24),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_min24x24_fast,
               (tTestEngTarget)&vec_max24x24_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_scale24x24,
               (tTestEngTarget)&vec_shift24x24,
               (tTestEngTarget)&vec_scale32x24),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsYvZ, &te_processFxn_vZvXsY ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_scale32x24_fast, (tTestEngTarget)&vec_scale24x24_fast,
                (tTestEngTarget)&vec_shift24x24_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsYvZ, &te_processFxn_vZvXsY ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_bexp24 ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_bexp24_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  {
    FUNC_LIST( (tTestEngTarget)&scl_bexp24 ),
    TEST_DESC( FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ, &processFxn_scl_vXvZ32 ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_scale32x24_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsYvZ, &te_processFxn_vZvXsY ) },
#endif
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot16x16 ),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYsZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot32x16 ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvY16sZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot32x32 ),
    TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYsZ64, &te_processFxn_sZ64vXvY) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot32x32_fast ),
    TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvYsZ64, &te_processFxn_sZ64vXvY) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot32x16_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvY16sZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot16x16_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvYsZ32, &te_processFxn_sZ32vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot64x32),
    TEST_DESC( FMT_REAL|FMT_INT64, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvY32sZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot64x32_fast),
    TEST_DESC( FMT_REAL|FMT_INT64, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvY32sZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot64x64, 
                (tTestEngTarget)&vec_dot64x64i),
    TEST_DESC( FMT_REAL|FMT_INT64, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYsZ64, &te_processFxn_sZ64vXvY ) },
   {
    FUNC_LIST(  (tTestEngTarget)&vec_dot64x64_fast,
                (tTestEngTarget)&vec_dot64x64i_fast),
    TEST_DESC( FMT_REAL|FMT_INT64, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvYsZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_add32x32 ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_add32x32_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_add16x16 ),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_add16x16_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_power16x16 ),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsY32sZ64, &te_processFxn_vXsY32sZ64 ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_power16x16_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsY32sZ64, &te_processFxn_vXsY32sZ64 ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_power32x32 ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsY32sZ64, &te_processFxn_vXsY32sZ64 ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_power32x32_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsY32sZ64, &te_processFxn_vXsY32sZ64 ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_min16x16, (tTestEngTarget)&vec_max16x16),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_min16x16_fast, (tTestEngTarget)&vec_max16x16_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_min32x32,(tTestEngTarget)&vec_max32x32),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_min32x32_fast,(tTestEngTarget)&vec_max32x32_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_scale32x32, 
                (tTestEngTarget)&vec_shift32x32),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsYvZ, &te_processFxn_vZvXsY ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_scale32x32_fast,
                (tTestEngTarget)&vec_shift32x32_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsYvZ, &te_processFxn_vZvXsY ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_scale16x16), 
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsYvZ, &te_processFxn_vZvXsY ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_scale16x16_fast), 
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsYvZ, &te_processFxn_vZvXsY ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_shift16x16), 
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsY32vZ, &te_processFxn_vZvXsY32 ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_shift16x16_fast), 
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsY32vZ, &te_processFxn_vZvXsY32 ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_bexp32),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_bexp32_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_bexp16),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_bexp16_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  {
    FUNC_LIST( (tTestEngTarget)&scl_bexp16),
    TEST_DESC( FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ32, &processFxn_scl_vXvZ32 ) },
  {
    FUNC_LIST((tTestEngTarget)&scl_bexp32 ),
    TEST_DESC( FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ, &processFxn_scl_vXvZ32 ) },

	/* MCRL Functions - Phase1   */
  {
	FUNC_LIST( (tTestEngTarget)&vec_eleabs32x32),
	TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvZ, &te_processFxn_vXvZ ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_elemax32x32 ),
	TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_elemin32x32 ),
	TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_elesub32x32 ),
	TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_elemult32x32 ),
	TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_sum32x32 ),
	TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXsZ32, &te_processFxn_vXsZ32 ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_mean32x32 ),
	TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXsZ32, &te_processFxn_vXsZ32 ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_rms32x32 ),
	TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXsZ32, &te_processFxn_vXsZ32 ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_var32x32 ),
 	TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXsZ32, &te_processFxn_vXsZ32 ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_stddev32x32 ),
	TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXsZ32, &te_processFxn_vXsZ32 ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_eleabs16x16),
	TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvZ, &te_processFxn_vXvZ ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_elemax16x16 ),
	TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_elemin16x16 ),
	TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_elesub16x16 ),
	TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_elemult16x16 ),
	TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_sum16x16 ),
	TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXsZ, &te_processFxn_vXsZ ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_mean16x16 ),
	TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXsZ, &te_processFxn_vXsZ ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_rms16x16 ),
	TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXsZ, &te_processFxn_vXsZ ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_var16x16 ),
	TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXsZ, &te_processFxn_vXsZ ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_stddev16x16 ),
	TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXsZ, &te_processFxn_vXsZ ) },
  {
    FUNC_LIST( NULL ), TEST_DESC(  0, 0, 0, NULL, NULL ) } /* End of table */
};
/* Perform all functional tests for Vector Mathematics API functions. */
int func_vector1(int isFull, int isVerbose, int breakOnError)
{
	int res = 1;

#define DO_TEST(fxn, seqFile)                                                                   \
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                              \
                                                       sizeof(testDefTbl)/sizeof(testDefTbl[0]),\
                                                       MAX_FUNC_NUM,                            \
                                                       (tTestEngTarget)(fxn), "vec1/" seqFile,   \
                                                       isFull, isVerbose, breakOnError ) )

  DO_TEST( &vec_dot16x16        , "vec_dot16x16.seq"        );
  DO_TEST( &vec_dot24x24        , "vec_dot24x24.seq"        );
  DO_TEST( &vec_dot32x16        , "vec_dot32x16.seq"        );
  DO_TEST( &vec_dot32x32        , "vec_dot32x32.seq"        );
  DO_TEST( &vec_dot64x32        , "vec_dot64x32.seq"        );
  DO_TEST( &vec_dot64x64        , "vec_dot64x64.seq"        );
  DO_TEST( &vec_dot64x64i       , "vec_dot64x64i.seq"       );
  DO_TEST( &vec_dot16x16_fast   , "vec_dot16x16_fast.seq"   );
  DO_TEST( &vec_dot24x24_fast   , "vec_dot24x24_fast.seq"   );
  DO_TEST( &vec_dot32x16_fast   , "vec_dot32x16_fast.seq"   );
  DO_TEST( &vec_dot32x32_fast   , "vec_dot32x32_fast.seq"   );
  DO_TEST( &vec_dot64x32_fast   , "vec_dot64x32_fast.seq"   );
  DO_TEST( &vec_dot64x64_fast   , "vec_dot64x64_fast.seq"   );
  DO_TEST( &vec_dot64x64i_fast  , "vec_dot64x64i_fast.seq"  );
  DO_TEST( &vec_add32x32        , "vec_add32x32.seq"        );
  DO_TEST( &vec_add24x24        , "vec_add24x24.seq"        );
  DO_TEST( &vec_add16x16        , "vec_add16x16.seq"        );
  DO_TEST( &vec_add32x32_fast   , "vec_add32x32_fast.seq"   );
  DO_TEST( &vec_add24x24_fast   , "vec_add24x24_fast.seq"   );
  DO_TEST( &vec_add16x16_fast   , "vec_add16x16_fast.seq"   );
  DO_TEST( &vec_power32x32      , "vec_power32x32.seq"      );
  DO_TEST( &vec_power24x24      , "vec_power24x24.seq"      );
  DO_TEST( &vec_power16x16      , "vec_power16x16.seq"      );
  DO_TEST( &vec_power32x32_fast , "vec_power32x32_fast.seq" );
  DO_TEST( &vec_power24x24_fast , "vec_power24x24_fast.seq" );
  DO_TEST( &vec_power16x16_fast , "vec_power16x16_fast.seq" );
  DO_TEST( &vec_shift32x32      , "vec_shift32x32.seq"      );
  DO_TEST( &vec_shift24x24      , "vec_shift24x24.seq"      );
  DO_TEST( &vec_shift16x16      , "vec_shift16x16.seq"      );
  DO_TEST( &vec_scale32x24      , "vec_scale32x24.seq"      );
  DO_TEST( &vec_scale24x24      , "vec_scale24x24.seq"      );
  DO_TEST( &vec_scale16x16      , "vec_scale16x16.seq"      );
  DO_TEST( &vec_scale32x32      , "vec_scale32x32.seq"      );
  DO_TEST( &vec_shift32x32_fast , "vec_shift32x32_fast.seq" );
  DO_TEST( &vec_shift24x24_fast , "vec_shift24x24_fast.seq" );
  DO_TEST( &vec_shift16x16_fast , "vec_shift16x16_fast.seq" );
  DO_TEST( &vec_scale32x24_fast , "vec_scale32x24_fast.seq" );
  DO_TEST( &vec_scale24x24_fast , "vec_scale24x24_fast.seq" );
  DO_TEST( &vec_scale16x16_fast , "vec_scale16x16_fast.seq" );
  DO_TEST( &vec_scale32x32_fast , "vec_scale32x32_fast.seq" );
  DO_TEST( &vec_min32x32        , "vec_min32x32.seq"        );
  DO_TEST( &vec_min24x24        , "vec_min24x24.seq"        );
  DO_TEST( &vec_min16x16        , "vec_min16x16.seq"        );
  DO_TEST( &vec_max32x32        , "vec_max32x32.seq"        );
  DO_TEST( &vec_max24x24        , "vec_max24x24.seq"        );
  DO_TEST( &vec_max16x16        , "vec_max16x16.seq"        );
  DO_TEST( &vec_min32x32_fast   , "vec_min32x32_fast.seq"   );
  DO_TEST( &vec_min24x24_fast   , "vec_min24x24_fast.seq"   );
  DO_TEST( &vec_min16x16_fast   , "vec_min16x16_fast.seq"   );
  DO_TEST( &vec_max32x32_fast   , "vec_max32x32_fast.seq"   );
  DO_TEST( &vec_max24x24_fast   , "vec_max24x24_fast.seq"   );
  DO_TEST( &vec_max16x16_fast   , "vec_max16x16_fast.seq"   );
  DO_TEST( &scl_bexp32             , "scl_bexp32.seq"       );
  DO_TEST( &scl_bexp24             , "scl_bexp24.seq"       );
  DO_TEST( &scl_bexp16             , "scl_bexp16.seq"       );
  DO_TEST( &vec_bexp32             , "vec_bexp32.seq"       );
  DO_TEST( &vec_bexp24             , "vec_bexp24.seq"       );
  DO_TEST( &vec_bexp16             , "vec_bexp16.seq"       );
  DO_TEST( &vec_bexp32_fast        , "vec_bexp32_fast.seq"  );
  DO_TEST( &vec_bexp24_fast        , "vec_bexp24_fast.seq"  );
  DO_TEST( &vec_bexp16_fast        , "vec_bexp16_fast.seq"  );

  /* MCRL */
  DO_TEST( &vec_eleabs32x32     , "vec_eleabs32x32.seq"     );
  DO_TEST( &vec_elemin32x32     , "vec_elemin32x32.seq"     );
  DO_TEST( &vec_elemax32x32     , "vec_elemax32x32.seq"     );
  DO_TEST( &vec_elesub32x32     , "vec_elesub32x32.seq"     );
  DO_TEST( &vec_elemult32x32    , "vec_elemult32x32.seq"    );
  DO_TEST( &vec_sum32x32        , "vec_sum32x32.seq"     	  );
  DO_TEST( &vec_mean32x32       , "vec_mean32x32.seq"       );
  DO_TEST( &vec_rms32x32        , "vec_rms32x32.seq"        );
  DO_TEST( &vec_var32x32        , "vec_var32x32.seq"        );
  DO_TEST( &vec_stddev32x32     , "vec_stddev32x32.seq"     );
  DO_TEST( &vec_eleabs16x16     , "vec_eleabs16x16.seq"     );
  DO_TEST( &vec_elemin16x16     , "vec_elemin16x16.seq"     );
  DO_TEST( &vec_elemax16x16     , "vec_elemax16x16.seq"     );
  DO_TEST( &vec_elesub16x16     , "vec_elesub16x16.seq"     );
  DO_TEST( &vec_elemult16x16    , "vec_elemult16x16.seq"    );
  DO_TEST( &vec_sum16x16        , "vec_sum16x16.seq"     	);
  DO_TEST( &vec_mean16x16       , "vec_mean16x16.seq"       );
  DO_TEST( &vec_rms16x16        , "vec_rms16x16.seq"        );
  DO_TEST( &vec_var16x16        , "vec_var16x16.seq"        );
  DO_TEST( &vec_stddev16x16     , "vec_stddev16x16.seq"     );

  return (res);
}

