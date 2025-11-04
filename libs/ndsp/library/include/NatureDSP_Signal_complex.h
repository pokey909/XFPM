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
#ifndef __NATUREDSP_SIGNAL_COMPLEX_H__
#define __NATUREDSP_SIGNAL_COMPLEX_H__

#include "NatureDSP_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
  Complex Mathematics:
  vec_complex2mag, vec_complex2invmag   Complex Magnitude
===========================================================================*/

/*-------------------------------------------------------------------------
  Complex magnitude
  Routines compute complex magnitude or its reciprocal

  Precision: 
  f     single precision floating point

  Input:
  x[N]  input complex data
  N     length of vector
  Output:
  y[N]  output data

  Restriction:
  none
-------------------------------------------------------------------------*/
void       vec_complex2mag    (float32_t  * y, const complex_float  * x, int N);
void       vec_complex2invmag (float32_t  * y, const complex_float  * x, int N);
float32_t  scl_complex2mag    (complex_float x);
float32_t  scl_complex2invmag (complex_float x);

void vec_cplx_Conjf(complex_float * r, const complex_float * x, int N);
void vec_cplxconj32x32(complex32_t*  z, const complex32_t*  x, int N);
void vec_cplxconj16x16(complex16_t*  z, const complex16_t*  x, int N);

void vec_cplx2cplx_multf (complex_float*  z, complex_float*  x, complex_float*  y, int N);
void vec_cplx2cplx_mult32x32 (complex32_t*  z, complex32_t*  x, complex32_t*  y, int N);
void vec_cplx2cplx_mult16x16 (complex16_t*  z, complex16_t*  x, complex16_t*  y, int N);

void vec_cplx2real_multvf (complex_float*  z, complex_float*  x, float32_t*  y, int N);
void vec_cplx2real_multv32x32 (complex32_t*  z, complex32_t*  x, int32_t*  y, int N);
void vec_cplx2real_multv16x16 (complex16_t*  z, complex16_t*  x, int16_t*  y, int N);

void vec_cplx2real_multsf (complex_float* restrict z, complex_float* restrict x, float32_t y, int N);
void vec_cplx2real_mults32x32 (complex32_t*  z, complex32_t*  x, int32_t  y, int N);
void vec_cplx2real_mults16x16 (complex16_t*  z, complex16_t*  x, int16_t  y, int N);

void vec_complex2mag32x32 (int32_t* z, const complex32_t* x, int N);
void vec_complex2mag16x16 (int16_t* z, const complex16_t* x, int N);

#ifdef __cplusplus
}
#endif

#endif/* __NATUREDSP_SIGNAL_COMPLEX_H__ */
