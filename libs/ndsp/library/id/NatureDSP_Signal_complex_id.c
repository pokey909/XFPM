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
 * Complex Math Functions
 * Annotations
 */

#include "NatureDSP_types.h"
#include "NatureDSP_Signal_complex.h"
#include "common.h"


ANNOTATE_FUN(vec_complex2mag,    "Vector Complex Magnitude (floating point data)            ");
ANNOTATE_FUN(vec_complex2invmag, "Vector Reciprocal Complex Magnitude (floating point data) ");
ANNOTATE_FUN(scl_complex2mag,    "Scalar Complex Magnitude (floating point data)            ");
ANNOTATE_FUN(scl_complex2invmag, "Scalar Reciprocal Complex Magnitude (floating point data) ");
ANNOTATE_FUN(vec_cplx2cplx_multf,   "Vector Shift with Saturation (floating point data)");
ANNOTATE_FUN(vec_cplx2cplx_mult32x32,	"Vector Complex to Complex Multiplication (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_cplx2cplx_mult16x16,   "Vector Complex to Complex Multiplication (16-bit input, 16-bit output)");
ANNOTATE_FUN(vec_cplx2real_multvf,      "Vector Complex to Real Multiplication with Vector (floating point data)");
ANNOTATE_FUN(vec_cplx2real_multsf,      "Vector Complex to Real Multiplication with Scalar (floating point data)");
ANNOTATE_FUN(vec_cplx2real_multv32x32,  "Vector Complex to Real Multiplication with Vector (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_cplx2real_mults32x32,  "Vector Complex to Real Multiplication with Scalar (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_cplx2real_multv16x16,  "Vector Complex to Real Multiplication with Vector (16-bit input, 16-bit output)");
ANNOTATE_FUN(vec_cplx2real_mults16x16,  "Vector Complex to Real Multiplication with Scalar (16-bit input, 16-bit output)");
ANNOTATE_FUN(vec_cplx_Conjf,        	"Vector Complex Conjugate Operation (floating point data)");
ANNOTATE_FUN(vec_cplxconj32x32,     "Vector Complex Conjugate Operation (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_cplxconj16x16,     "Vector Complex Conjugate Operation (16-bit input, 16-bit output)");
ANNOTATE_FUN(vec_complex2mag32x32,  "Vector complex magnitude (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_complex2mag16x16,  "Vector complex magnitude (16-bit input, 16-bit output)");
