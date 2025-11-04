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
* Test module for testing cycle performance (Vector Operations)
*/
#include "types.h"
#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(vector)
#include "mips.h"


void mips_vector2(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_NORMALIZED(    1, isVerbose, vec_dotf,       (mips.inp0.f32,mips.inp1.f32,200           ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(    1, isVerbose, vec_addf,       (mips.out0.f32,mips.inp0.f32,mips.inp1.f32,200  ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(    1, isVerbose, vec_powerf,     (mips.inp0.f32,200                    ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(    1, isVerbose, vec_shiftf,     (mips.out0.f32, mips.inp0.f32, 1, 200      ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(    1, isVerbose, vec_scalef,     (mips.out0.f32, mips.inp0.f32, 1.   , 200  ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(    1, isVerbose, vec_scale_sf,   (mips.out0.f32, mips.inp0.f32, 32767 ,- 1000, 1000, 200), fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(    1, isVerbose, vec_minf,       (mips.inp0.f32,200                    ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(    1, isVerbose, vec_maxf,       (mips.inp0.f32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(    1, isVerbose, vec_bexpf,      (mips.inp0.f32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_SIMPLE(     1, isVerbose, scl_bexpf,(9621325.f),      fout,""     ,prf_cycle);
    PROFILE_NORMALIZED(1, isVerbose, vec_eleabsf			, (mips.out0.f32, mips.inp0.f32, 				200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_elemaxf			, (mips.out0.f32, mips.inp0.f32, mips.inp1.f32, 200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_eleminf			, (mips.out0.f32, mips.inp0.f32, mips.inp1.f32, 200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_elesubf			, (mips.out0.f32, mips.inp0.f32, mips.inp1.f32, 200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_sumf				, (mips.inp0.f32,                				200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_elemultf			, (mips.out0.f32, mips.inp0.f32, mips.inp1.f32, 200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_meanf				, (mips.inp0.f32, 								200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_rmsf				, (mips.inp0.f32, 								200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_varf				, (mips.inp0.f32, 								200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_stddevf			, (mips.inp0.f32, 								200     ), fout,"N=200", prf_cyclespts, 200);
}

