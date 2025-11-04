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
 * Test procedures for matrix functions
 */

#include <string.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API. */
#include LIBRARY_HEADER(matop)
/* Test engine API. */
#include "../common/testeng_matop.h"

#define MAX_FUNC_NUM   10
/* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
#define FUNC_LIST(...) { __VA_ARGS__, NULL }
/* Initializer for a test description structure. */
#define TEST_DESC( fmt, extra, api, dimNum, align, loadFxn, procFxn ) { (fmt),extra,api,(dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }

/* vec API test definitions. */
static const struct 
{
tTestEngTarget   funcList[MAX_FUNC_NUM];
tTestEngDesc     testDesc;
}
testDefTbl[] =
{
    { FUNC_LIST( (tTestEngTarget)&mtx_mpyf ),          TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,NULL  ,TE_DIM_NUM_4, TE_ALIGN_NO , &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpyf ),       TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,NULL  ,TE_DIM_NUM_4, TE_ALIGN_NO , &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpyf_fast),      TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,NULL  ,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpyf_fast ),  TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,NULL  ,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
    { FUNC_LIST( NULL ), TEST_DESC(  0, 0, 0, 0, 0, NULL, NULL ) } /* End of table */
};

/* Perform all tests for mat API functions. */
int func_matop2(int isFull, int isVerbose, int breakOnError)
{
  int res = 1;

  #define DO_TEST(fxn, seqFile)                                                                  \
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                               \
                                                       sizeof(testDefTbl)/sizeof(testDefTbl[0]), \
                                                       MAX_FUNC_NUM,                             \
                                                       (tTestEngTarget)(fxn), "mtx2/" seqFile,    \
                                                       isFull, isVerbose, breakOnError ) )

    DO_TEST( &mtx_mpyf           , "mtx_mpyf.seq"           );
    DO_TEST( &mtx_mpyf_fast      , "mtx_mpyf_fast.seq"      );
    DO_TEST( &mtx_vecmpyf        , "mtx_vecmpyf.seq"        );
    DO_TEST( &mtx_vecmpyf_fast   , "mtx_vecmpyf_fast.seq"   );

    return (res);
}
