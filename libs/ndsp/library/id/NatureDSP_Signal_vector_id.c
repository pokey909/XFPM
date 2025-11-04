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
 * Vector Operations
 * Annotations
 */

#include "NatureDSP_types.h"
#include "NatureDSP_Signal_vector.h"
#include "common.h"

ANNOTATE_FUN(vec_dot24x24,        "Vector Dot product (24x24-bit data, 64-bit output)");
ANNOTATE_FUN(vec_dot32x16,        "Vector Dot product (32x16-bit data, 64-bit output)");
ANNOTATE_FUN(vec_dot32x32,        "Vector Dot product (32x32-bit data, 64-bit output)");
ANNOTATE_FUN(vec_dot16x16,        "Vector Dot product (16x16-bit data, 64-bit output)");
ANNOTATE_FUN(vec_dot64x32,        "Vector Dot product (64x32-bit data, 64-bit output)");
ANNOTATE_FUN(vec_dot64x64,        "Vector Dot product (64x64-bit data, 64-bit output)");
ANNOTATE_FUN(vec_dot64x64i,       "Vector Dot product (64x64-bit data, 64-bit output (low 64 bit of integer multiply))");
ANNOTATE_FUN(vec_dotf,            "Vector Dot product (floating point data)");
ANNOTATE_FUN(vec_dot24x24_fast,   "Fast Vector Dot product (24x24-bit data, 64-bit output)");
ANNOTATE_FUN(vec_dot32x16_fast,   "Fast Vector Dot product (32x16-bit data, 64-bit output)");
ANNOTATE_FUN(vec_dot32x32_fast,   "Fast Vector Dot product (32x32-bit data, 64-bit output)");
ANNOTATE_FUN(vec_dot16x16_fast,   "Fast Vector Dot product (16x16-bit data, 32-bit output)");
ANNOTATE_FUN(vec_dot64x32_fast,   "Fast Vector Dot product (64x32-bit data, 64-bit output)");
ANNOTATE_FUN(vec_dot64x64_fast,   "Fast Vector Dot product (64x64-bit data, 64-bit output)");
ANNOTATE_FUN(vec_dot64x64i_fast,  "Fast Vector Dot product (64x64-bit data, 64-bit output (low 64 bit of integer multiply))");
ANNOTATE_FUN(vec_add32x32,        "Vector Sum (32-bit data)");
ANNOTATE_FUN(vec_add32x32_fast,   "Fast Vector Sum (32-bit data)");
ANNOTATE_FUN(vec_add24x24,        "Vector Sum (24-bit data)");
ANNOTATE_FUN(vec_addf,            "Vector Sum (floating point data)");
ANNOTATE_FUN(vec_add24x24_fast,   "Fast Vector Sum (24-bit data)");
ANNOTATE_FUN(vec_add16x16,        "Vector Sum (16-bit data)");
ANNOTATE_FUN(vec_add16x16_fast,   "Fast Vector Sum (16-bit data)");
ANNOTATE_FUN(vec_power24x24,      "Power of a Vector (24x24-bit data, 64-bit output)");
ANNOTATE_FUN(vec_power32x32,      "Power of a Vector (32x32-bit data, 64-bit output)");
ANNOTATE_FUN(vec_power16x16,      "Power of a Vector (16x16-bit data, 64-bit output)");
ANNOTATE_FUN(vec_powerf,          "Power of a Vector (floating point data)");
ANNOTATE_FUN(vec_power24x24_fast, "Fast Power of a Vector (24x24-bit data, 64-bit output)");
ANNOTATE_FUN(vec_power32x32_fast, "Fast Power of a Vector (32x32-bit data, 64-bit output)");
ANNOTATE_FUN(vec_power16x16_fast, "Fast Power of a Vector (16x16-bit data, 64-bit output)");
ANNOTATE_FUN(vec_shift24x24,      "Vector Shift with Saturation (24-bit data)");
ANNOTATE_FUN(vec_shift32x32,      "Vector Shift with Saturation (32-bit data)");
ANNOTATE_FUN(vec_shift16x16,      "Vector Shift with Saturation (16-bit data)");
ANNOTATE_FUN(vec_shiftf,          "Vector Shift with Saturation (floating point data)");
ANNOTATE_FUN(vec_scale32x24,      "Vector Scaling with Saturation (32-bit input, 32-bit output, 24-bit scale factor)");
ANNOTATE_FUN(vec_scale32x32,      "Vector Scaling with Saturation (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_scale24x24,      "Vector Scaling with Saturation (24-bit input, 24-bit output)");
ANNOTATE_FUN(vec_scale16x16,      "Vector Scaling with Saturation (16-bit input, 16-bit output)");
ANNOTATE_FUN(vec_scalef,          "Vector Scaling with Saturation (floating point data)");
ANNOTATE_FUN(vec_scale_sf,        "Vector Scaling with Saturation (floating point data)");
ANNOTATE_FUN(vec_shift24x24_fast, "Fast Vector Shift with Saturation (24-bit data)");
ANNOTATE_FUN(vec_shift32x32_fast, "Fast Vector Shift with Saturation (32-bit data)");
ANNOTATE_FUN(vec_shift16x16_fast, "Fast Vector Shift with Saturation (16-bit data)");
ANNOTATE_FUN(vec_scale32x24_fast, "Fast Vector Scaling with Saturation (32-bit input, 32-bit output, 24-bit scale factor)");
ANNOTATE_FUN(vec_scale32x32_fast, "Fast Vector Scaling with Saturation (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_scale24x24_fast, "Fast Vector Scaling with Saturation (24-bit input, 24-bit output)");
ANNOTATE_FUN(vec_scale16x16_fast, "Fast Vector Scaling with Saturation (16-bit input, 16-bit output)");
ANNOTATE_FUN(vec_bexp32,          "Common Exponent (32-bit input data)");
ANNOTATE_FUN(vec_bexp24,          "Common Exponent (24-bit input data)");
ANNOTATE_FUN(vec_bexp16,          "Common Exponent (16-bit input data)");
ANNOTATE_FUN(vec_bexpf,           "Common Exponent (floating point input data)");
ANNOTATE_FUN(scl_bexp32,          "Exponent (32-bit input data)");
ANNOTATE_FUN(scl_bexp24,          "Exponent (24-bit input data)");
ANNOTATE_FUN(scl_bexp16,          "Exponent (16-bit input data)");
ANNOTATE_FUN(scl_bexpf,           "Exponent (floating point input data)");
ANNOTATE_FUN(vec_bexp32_fast,     "Fast Common Exponent (32-bit input data)");
ANNOTATE_FUN(vec_bexp24_fast,     "Fast Common Exponent (24-bit input data)");
ANNOTATE_FUN(vec_bexp16_fast,     "Fast Common Exponent (16-bit input data)");
ANNOTATE_FUN(vec_min32x32,        "Vector Minimum Value (32-bit data)");
ANNOTATE_FUN(vec_min24x24,        "Vector Minimum Value (24-bit data)");
ANNOTATE_FUN(vec_min16x16,        "Vector Minimum Value (16-bit data)");
ANNOTATE_FUN(vec_minf,            "Vector Minimum Value (floating point data)");
ANNOTATE_FUN(vec_max32x32,        "Vector Maximum Value (32-bit data)");
ANNOTATE_FUN(vec_max24x24,        "Vector Maximum Value (24-bit data)");
ANNOTATE_FUN(vec_max16x16,        "Vector Maximum Value (16-bit data)");
ANNOTATE_FUN(vec_maxf,            "Vector Maximum Value (floating point data)");
ANNOTATE_FUN(vec_min32x32_fast,   "Fast Vector Minimum Value (32-bit data)");
ANNOTATE_FUN(vec_min24x24_fast,   "Fast Vector Minimum Value (24-bit data)");
ANNOTATE_FUN(vec_min16x16_fast,   "Fast Vector Minimum Value (16-bit data)");
ANNOTATE_FUN(vec_max32x32_fast,   "Fast Vector Maximum Value (32-bit data)");
ANNOTATE_FUN(vec_max24x24_fast,   "Fast Vector Maximum Value (24-bit data)");
ANNOTATE_FUN(vec_max16x16_fast,   "Fast Vector Maximum Value (16-bit data)");

ANNOTATE_FUN(vec_eleabsf,         "Vector element-wise absolute (floating point data)");
ANNOTATE_FUN(vec_eleabs32x32,     "Vector element-wise absolute (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_eleabs16x16,     "Vector element-wise absolute (16-bit input, 16-bit output)");
ANNOTATE_FUN(vec_elemaxf,         "Vector element-wise maximum (floating point data)");
ANNOTATE_FUN(vec_elemax32x32,     "Vector element-wise maximum (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_elemax16x16,     "Vector element-wise maximum (16-bit input, 16-bit output)");
ANNOTATE_FUN(vec_eleminf,         "Vector element-wise minimum (floating point data)");
ANNOTATE_FUN(vec_elemin32x32,     "Vector element-wise minimum (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_elemin16x16,     "Vector element-wise minimum (16-bit input, 16-bit output)");
ANNOTATE_FUN(vec_elesubf,         "Vector element-wise subtraction (floating point data)");
ANNOTATE_FUN(vec_elesub32x32,     "Vector element-wise subtraction (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_elesub16x16,     "Vector element-wise subtraction (16-bit input, 16-bit output)");
ANNOTATE_FUN(vec_elemultf,        "Vector element-wise multiplication (floating point data)");
ANNOTATE_FUN(vec_elemult32x32,    "Vector element-wise multiplication  (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_elemult16x16,    "Vector element-wise multiplication  (16-bit input, 16-bit output)");

ANNOTATE_FUN(vec_sumf,            "Vector sum (floating point data)");
ANNOTATE_FUN(vec_sum32x32,        "Vector sum (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_sum16x16,        "Vector sum (16-bit input, 16-bit output)");
ANNOTATE_FUN(vec_meanf,        	  "Vector mean (floating point data)");
ANNOTATE_FUN(vec_mean32x32,       "Vector mean (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_mean16x16,       "Vector mean (16-bit input, 16-bit output)");
ANNOTATE_FUN(vec_rmsf,        	  "Vector rms (floating point data)");
ANNOTATE_FUN(vec_rms32x32,        "Vector rms  (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_rms16x16,        "Vector rms  (16-bit input, 16-bit output)");
ANNOTATE_FUN(vec_varf,            "Vector variance  (floating point data)");
ANNOTATE_FUN(vec_var32x32,        "Vector variance (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_var16x16,        "Vector variance (16-bit input, 16-bit output)");
ANNOTATE_FUN(vec_stddevf,         "Vector Standard Deviation (floating point data)");
ANNOTATE_FUN(vec_stddev32x32,     "Vector Standard Deviation (32-bit input, 32-bit output)");
ANNOTATE_FUN(vec_stddev16x16,     "Vector Standard Deviation (16-bit input, 16-bit output)");
