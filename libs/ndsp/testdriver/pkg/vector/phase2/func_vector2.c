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
 * Test procedures for arithmetic and logic functions on data vectors
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
  * Stage 2
  */
  {
    FUNC_LIST( (tTestEngTarget)&vec_dotf ),
    TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYsZ, &te_processFxn_vXvYsZ ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_addf),
    TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_scalef ),
    TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsYvZ, &te_processFxn_vZvXsY ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_shiftf), 
    TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsY32vZ, &te_processFxn_vZvXsY32 ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_minf, (tTestEngTarget)&vec_maxf),
    TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_powerf ),
    TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ, &te_processFxn_vXsZ ) },
  {
    FUNC_LIST((tTestEngTarget)&vec_scale_sf),
    TEST_DESC(FMT_REAL | FMT_FLOAT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_vector_loadFxn_vXsF0sF1sYvZ, &te_vector_processFxn_vXsF0sF1sYvZ) },
  {
    FUNC_LIST( (tTestEngTarget)&scl_bexpf),
    TEST_DESC( FMT_FLOAT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ32, &processFxn_scl_vXvZ32 ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_bexpf),
    TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ32  , &te_processFxn_vXsZ32 ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_eleabsf),
	TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvZ, &te_processFxn_vXvZ ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_elemaxf ),
	TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_eleminf ),
	TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_elesubf ),
	TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_elemultf ),
	TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_sumf ),
	TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXsZ, &te_processFxn_vXsZ ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_meanf ),
	TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXsZ, &te_processFxn_vXsZ ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_rmsf ),
	TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXsZ, &te_processFxn_vXsZ ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_varf ),
	TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXsZ, &te_processFxn_vXsZ ) },
  {
	FUNC_LIST( (tTestEngTarget)&vec_stddevf ),
	TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXsZ, &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( NULL ), TEST_DESC(  0, 0, 0, NULL, NULL ) } /* End of table */
};

/* Perform all functional tests for Vector Mathematics API functions. */
int func_vector2(int isFull, int isVerbose, int breakOnError)
{
	int res = 1;

#define DO_TEST(fxn, seqFile)																	\
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                              \
                                                       sizeof(testDefTbl)/sizeof(testDefTbl[0]),\
                                                       MAX_FUNC_NUM,                            \
                                                       (tTestEngTarget)(fxn), "vec2/" seqFile,   \
                                                       isFull, isVerbose, breakOnError ) )
    DO_TEST(&vec_dotf			, "vec_dotf.seq");
    DO_TEST(&vec_addf			, "vec_addf.seq");
    DO_TEST(&vec_powerf			, "vec_powerf.seq");
    DO_TEST(&vec_shiftf			, "vec_shiftf.seq");
    DO_TEST(&vec_scalef			, "vec_scalef.seq");
    DO_TEST(&vec_scale_sf		, "vec_scale_sf.seq");
    DO_TEST(&vec_minf			, "vec_minf.seq");
    DO_TEST(&vec_maxf			, "vec_maxf.seq");
    DO_TEST(&scl_bexpf			, "scl_bexpf.seq");
    DO_TEST(&vec_bexpf			, "vec_bexpf.seq");

    DO_TEST(&vec_eleabsf          	, "vec_eleabsf.seq" 	);
	DO_TEST(&vec_elemaxf          	, "vec_elemaxf.seq"    	);
	DO_TEST(&vec_eleminf          	, "vec_eleminf.seq"    	);
	DO_TEST(&vec_elesubf          	, "vec_elesubf.seq"     );
	DO_TEST(&vec_elemultf         	, "vec_elemultf.seq"    );
	DO_TEST(&vec_sumf             	, "vec_sumf.seq"        );
	DO_TEST(&vec_meanf            	, "vec_meanf.seq"     	);
	DO_TEST(&vec_rmsf             	, "vec_rmsf.seq"       	);
	DO_TEST(&vec_varf             	, "vec_varf.seq"       	);
	DO_TEST(&vec_stddevf          	, "vec_stddevf.seq"    	);
	return (res);
}

