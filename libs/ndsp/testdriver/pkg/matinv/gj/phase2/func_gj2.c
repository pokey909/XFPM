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
 * Test procedures for matrix inversion and related functions
 */

#include <string.h>
#include <stdlib.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API. */
#include LIBRARY_HEADER(matinv)
/* Test engine API. */
#include "../../common/testeng_matinv.h"

static const tMtxInvApi mtx_inv2x2f_Api    =      {1,NULL };
static const tMtxInvApi mtx_inv3x3f_Api    =      {1,NULL };
static const tMtxInvApi mtx_inv4x4f_Api    =      {1,NULL };

static const tTestEngDesc descr_mtx_inv2x2f =  { FMT_REAL | FMT_FLOAT32, MTXINV_PLAIN,&mtx_inv2x2f_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxinv, &te_processFxn_matinv };
static const tTestEngDesc descr_mtx_inv3x3f =  { FMT_REAL | FMT_FLOAT32, MTXINV_PLAIN,&mtx_inv3x3f_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxinv, &te_processFxn_matinv };
static const tTestEngDesc descr_mtx_inv4x4f =  { FMT_REAL | FMT_FLOAT32, MTXINV_PLAIN,&mtx_inv4x4f_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxinv, &te_processFxn_matinv };
static const tTestMtxinv tests[] = 
{
  { &descr_mtx_inv2x2f,   (tTestEngTarget)&mtx_inv2x2f,   1,"gj2/mtx_inv2x2f_cond10.seq"      },
  { &descr_mtx_inv2x2f,   (tTestEngTarget)&mtx_inv2x2f,   1,"gj2/mtx_inv2x2f_cond100.seq"     },
  { &descr_mtx_inv2x2f,   (tTestEngTarget)&mtx_inv2x2f,   1,"gj2/mtx_inv2x2f_cond1000.seq"    },
  { &descr_mtx_inv3x3f,   (tTestEngTarget)&mtx_inv3x3f,   1,"gj2/mtx_inv3x3f_cond10.seq"      },
  { &descr_mtx_inv3x3f,   (tTestEngTarget)&mtx_inv3x3f,   1,"gj2/mtx_inv3x3f_cond100.seq"     },
  { &descr_mtx_inv3x3f,   (tTestEngTarget)&mtx_inv3x3f,   1,"gj2/mtx_inv3x3f_cond1000.seq"    },
  { &descr_mtx_inv4x4f,   (tTestEngTarget)&mtx_inv4x4f,   1,"gj2/mtx_inv4x4f_cond10.seq"      },
  { &descr_mtx_inv4x4f,   (tTestEngTarget)&mtx_inv4x4f,   1,"gj2/mtx_inv4x4f_cond100.seq"     },
  { &descr_mtx_inv4x4f,   (tTestEngTarget)&mtx_inv4x4f,   1,"gj2/mtx_inv4x4f_cond1000.seq"    },

  { NULL }  /* end of list */ 
};
/* Perform all tests for matrix inversion API functions. */
int func_gj2(int isFull, int isVerbose, int breakOnError)
{
  int res = 1;
  res &= (0!=te_ExecMtxInv(tests,isFull,isVerbose,breakOnError));
  return res;
}
