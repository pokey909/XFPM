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
 * NatureDSP_Signal Library API
 * Matrix Operations
 * Annotations
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Signal_matop.h"
#include "common.h"

ANNOTATE_FUN(mtx_mpy24x24,      "Matrix Multiply (24-bit data)");
ANNOTATE_FUN(mtx_mpy32x32,      "Matrix Multiply (32-bit data)");
ANNOTATE_FUN(mtx_mpy16x16,      "Matrix Multiply (16-bit data)");
ANNOTATE_FUN(mtx_mpyf,          "Matrix Multiply (floating point data)");
ANNOTATE_FUN(mtx_mpy24x24_fast, "Fast Matrix Multiply (24-bit data)");
ANNOTATE_FUN(mtx_mpy32x32_fast, "Fast Matrix Multiply (32-bit data)");
ANNOTATE_FUN(mtx_mpy16x16_fast, "Fast Matrix Multiply (16-bit data)");
ANNOTATE_FUN(mtx_mpyf_fast,     "Fast Matrix Multiply (floating point data)");
ANNOTATE_FUN(mtx_vecmpy32x32,       "Matrix by Vector Multiply (32-bit data)");
ANNOTATE_FUN(mtx_vecmpy24x24,       "Matrix by Vector Multiply (24-bit data)");
ANNOTATE_FUN(mtx_vecmpy16x16,       "Matrix by Vector Multiply (16-bit data)");
ANNOTATE_FUN(mtx_vecmpy32x32_fast,  "Fast Matrix by Vector Multiply (32-bit data)");
ANNOTATE_FUN(mtx_vecmpy24x24_fast,  "Fast Matrix by Vector Multiply (24-bit data)");
ANNOTATE_FUN(mtx_vecmpy16x16_fast,  "Fast Matrix by Vector Multiply (16-bit data)");
ANNOTATE_FUN(mtx_vecmpyf,           "Matrix by Vector Multiply (floating point data)");
ANNOTATE_FUN(mtx_vecmpyf_fast,      "Fast Matrix by Vector Multiply (floating point data)");
