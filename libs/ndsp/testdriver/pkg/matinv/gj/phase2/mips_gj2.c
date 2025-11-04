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
* Test module for testing cycle performance (Matrix Decomposition 
* and Inversion)
*/

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API. */
#include LIBRARY_HEADER(matinv)
/* MIPS measurement means. */
#include "mips.h"
/* Utility functions and macros. */
#include "utils.h"

#define PROFILE_MATINV(cond,verb,fun,N,suffix)          PROFILE_NORMALIZED(cond,verb,fun,(mips.inp0.suffix   ),fout,"",prf_cyclesmtx,1);

void mips_gj2(int isFull, int isVerbose, FILE * fout)
{
    int i;
    /* fill with random floating point data */
    Rand_reset(12,273);
    for (i=0; i<200; i++)
    {
        mips.inp0.f32[i]=((int16_t)Rand())*(1.f/32768.f);
    }
    PROFILE_MATINV(     1, isVerbose, mtx_inv2x2f  , 2,f32);
    PROFILE_MATINV(     1, isVerbose, mtx_inv3x3f  , 3,f32);
    PROFILE_MATINV(     1, isVerbose, mtx_inv4x4f  , 4,f32);
}
