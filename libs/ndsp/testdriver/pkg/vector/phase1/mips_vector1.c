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

typedef void tFxn(void   *z, const void* x, const void* y, int rsh, int N, int M);

void mips_vector1(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot16x16,       (mips.inp0.i16,mips.inp1.i16,200           ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot24x24,       (mips.inp0.i32, mips.inp1.i32 , 200        ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot32x16,       (mips.inp0.i32,mips.inp1.i16  ,200         ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot32x16,       (mips.inp0.i32,mips.inp1.i16+1,200         ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot32x32,       (mips.inp0.i32, mips.inp1.i32, 200         ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot32x32,       (mips.inp0.i32+1, mips.inp1.i32, 200       ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot64x32,       (mips.inp0.i64, mips.inp1.i32, 200         ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot64x64,       (mips.inp0.i64, mips.inp1.i64, 200         ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot64x64i,      (mips.inp0.i64, mips.inp1.i64, 200         ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_dot16x16_fast , (mips.inp0.i16,mips.inp1.i16,200           ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_dot24x24_fast,  (mips.inp0.i32, mips.inp1.i32, 200         ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_dot32x16_fast , (mips.inp0.i32,mips.inp1.i16,200           ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_dot32x32_fast,  (mips.inp0.i32, mips.inp1.i32, 200         ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_dot64x32_fast,  (mips.inp0.i64, mips.inp1.i32, 200         ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_dot64x64_fast,  (mips.inp0.i64, mips.inp1.i64, 200         ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_dot64x64i_fast, (mips.inp0.i64, mips.inp1.i64, 200         ),fout,"N=200",prf_cyclespts,200);

    PROFILE_NORMALIZED(isFull, isVerbose, vec_add16x16     ,  (mips.out0.i16,mips.inp0.i16,mips.inp1.i16,200  ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_add16x16     ,  (mips.out0.i16,mips.inp0.i16+1,mips.inp1.i16,200  ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_add24x24,       (mips.out0.i32,mips.inp0.i32,mips.inp1.i32,200  ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_add32x32     ,  (mips.out0.i32,mips.inp0.i32,mips.inp1.i32,200  ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_add32x32     ,  (mips.out0.i32,mips.inp0.i32+1,mips.inp1.i32,200  ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_add16x16_fast,  (mips.out0.i16,mips.inp0.i16,mips.inp1.i16,200  ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_add24x24_fast,  (mips.out0.i32,mips.inp0.i32,mips.inp1.i32, 200 ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_add32x32_fast,  (mips.out0.i32,mips.inp0.i32,mips.inp1.i32,200  ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_power16x16,     (mips.inp0.i16,24,200                 ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_power16x16,     (mips.inp0.i16+1,24,200               ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_power24x24,     (mips.inp0.i32, 44, 200               ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_power32x32,     (mips.inp0.i32,44,200                 ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_power32x32,     (mips.inp0.i32+1,44,200               ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_power16x16_fast,(mips.inp0.i16,24,200                 ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_power24x24_fast,(mips.inp0.i32, 44, 200               ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_power32x32_fast,(mips.inp0.i32,44,200                 ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_shift16x16,     (mips.out0.i16,mips.inp0.i16, 1,200        ),fout,"shift>0, x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_shift16x16,     (mips.out0.i16+1,mips.inp0.i16+1, 1,200    ),fout,"shift>0, x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_shift16x16,     (mips.out0.i16,mips.inp0.i16,-1,200        ),fout,"shift<0, x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_shift16x16,     (mips.out0.i16+1,mips.inp0.i16+1, -1,200    ),fout,"shift<0, x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_shift24x24,     (mips.out0.i32, mips.inp0.i32, 1, 200      ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_shift32x32,     (mips.out0.i32, mips.inp0.i32, 1, 200      ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_shift32x32,     (mips.out0.i32+1, mips.inp0.i32+1, 1, 200      ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_scale16x16,     (mips.out0.i16,mips.inp0.i16,32767,200     ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_scale16x16,     (mips.out0.i16+1,mips.inp0.i16+1,32767,200     ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_scale24x24,     (mips.out0.i32, mips.inp0.i32, 32767, 200  ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_scale32x24,     (mips.out0.i32, mips.inp0.i32, 32767, 200  ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_scale32x32,     (mips.out0.i32, mips.inp0.i32, 32767, 200  ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_scale32x32,     (mips.out0.i32, mips.inp0.i32+1, 32767, 200),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_shift16x16_fast,(mips.out0.i16,mips.inp0.i16, 1,200        ),fout,"shift>0, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_shift16x16_fast,(mips.out0.i16,mips.inp0.i16,-1,200        ),fout,"shift<0, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_shift32x32_fast,(mips.out0.i32, mips.inp0.i32, 1, 200      ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_scale16x16_fast,(mips.out0.i16,mips.inp0.i16,32767,200     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_scale24x24_fast,(mips.out0.i32, mips.inp0.i32, 32767, 200  ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_scale32x24_fast,(mips.out0.i32, mips.inp0.i32, 32767, 200  ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_scale32x32_fast,(mips.out0.i32, mips.inp0.i32, 32767, 200  ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_shift24x24_fast,(mips.out0.i32, mips.inp0.i32, 1, 200      ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_max16x16,       (mips.inp0.i16,200                    ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_max16x16,       (mips.inp0.i16+1,200                    ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_min16x16,       (mips.inp0.i16,200                    ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_min16x16,       (mips.inp0.i16+1,200                    ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_max24x24,       (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_min24x24,       (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_max32x32,       (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_min32x32,       (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_max16x16_fast,  (mips.inp0.i16, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_min16x16_fast,  (mips.inp0.i16, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_max24x24_fast,  (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_min24x24_fast,  (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_max32x32_fast,  (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_min32x32_fast,  (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_bexp16,  (mips.inp0.i16, 200                          ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_bexp24,  (mips.inp0.i32, 200                          ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_bexp32,  (mips.inp0.i32, 200                          ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_bexp16_fast,  (mips.inp0.i16, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_bexp24_fast,  (mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_bexp32_fast,  (mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_SIMPLE(     1, isVerbose, scl_bexp16,(3917),                           fout,"",prf_cycle);
    PROFILE_SIMPLE(     1, isVerbose, scl_bexp24,(9621325),                        fout,"",prf_cycle);
    PROFILE_SIMPLE(     1, isVerbose, scl_bexp32,(9621325),                        fout,"",prf_cycle);
    /* MCRL */
    PROFILE_NORMALIZED(1, isVerbose, vec_eleabs16x16		, (mips.out0.i16, mips.inp0.i16, 				200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_eleabs32x32		, (mips.out0.i32, mips.inp0.i32, 				200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_elemax16x16		, (mips.out0.i16, mips.inp0.i16, mips.inp1.i16, 200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_elemax32x32		, (mips.out0.i32, mips.inp0.i32, mips.inp1.i32, 200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_elemin16x16		, (mips.out0.i16, mips.inp0.i16, mips.inp1.i16, 200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_elemin32x32		, (mips.out0.i32, mips.inp0.i32, mips.inp1.i32, 200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_elesub16x16		, (mips.out0.i16, mips.inp0.i16, mips.inp1.i16, 200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_elesub32x32		, (mips.out0.i32, mips.inp0.i32, mips.inp1.i32, 200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_elemult16x16		, (mips.out0.i16, mips.inp0.i16, mips.inp1.i16, 200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_elemult32x32		, (mips.out0.i32, mips.inp0.i32, mips.inp1.i32, 200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_sum16x16			, (mips.inp0.i16,                				200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_sum32x32			, (mips.inp0.i32,                				200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_mean16x16			, (mips.inp0.i16, 								200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_mean32x32			, (mips.inp0.i32, 								200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_rms16x16			, (mips.inp0.i16, 								200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_rms32x32			, (mips.inp0.i32, 								200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_var16x16			, (mips.inp0.i16, 								200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_var32x32			, (mips.inp0.i32, 								200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_stddev16x16		, (mips.inp0.i16, 								200     ), fout,"N=200", prf_cyclespts, 200);
    PROFILE_NORMALIZED(1, isVerbose, vec_stddev32x32		, (mips.inp0.i32, 								200     ), fout,"N=200", prf_cyclespts, 200);
}
