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
* Test module for testing cycle performance (Matrix Operations)
*/
#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(matop)
#include "mips.h"

#define PROFILE_mpy16x16(cond,verb,M,N,P)      PROFILE_INVERTED(cond,verb,mtx_mpy16x16,(pScr,mips.out0.i16, mips.inp0.i16,mips.inp1.i16, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpy16x16_fast(cond,verb,M,N,P) PROFILE_INVERTED(cond,verb,mtx_mpy16x16_fast,(mips.out0.i16, mips.inp0.i16, mips.inp1.i16, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpy24x24(cond,verb,M,N,P)      PROFILE_INVERTED(cond,verb,mtx_mpy24x24,(pScr,mips.out0.i32, mips.inp0.i32, mips.inp1.i32, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpy24x24_fast(cond,verb,M,N,P) PROFILE_INVERTED(cond,verb,mtx_mpy24x24_fast,(mips.out0.i32, mips.inp0.i32, mips.inp1.i32, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpy32x32(cond,verb,M,N,P)      PROFILE_INVERTED(cond,verb,mtx_mpy32x32,(pScr,mips.out0.i32, mips.inp0.i32, mips.inp1.i32, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpy32x32_fast(cond,verb,M,N,P) PROFILE_INVERTED(cond,verb,mtx_mpy32x32_fast,(mips.out0.i32, mips.inp0.i32, mips.inp1.i32, M, N, P,0),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));

#define PROFILE_mtx_vecmpy16x16(cond,verb,fun,M,N) PROFILE_INVERTED(cond,verb,fun,(mips.out0.i16, mips.inp0.i16, mips.inp1.i16, M, N,0),fout,#M"x"#N" x "#N"x1",prf_maccycle, (N*M));
#define PROFILE_mtx_vecmpy24x24(cond,verb,fun,M,N) PROFILE_INVERTED(cond,verb,fun,(mips.out0.i32, mips.inp0.i32, mips.inp1.i32, M, N,0),fout,#M"x"#N" x "#N"x1",prf_maccycle, (N*M));
#define PROFILE_mtx_vecmpy32x32(cond,verb,fun,M,N) PROFILE_INVERTED(cond,verb,fun,(mips.out0.i32, mips.inp0.i32, mips.inp1.i32, M, N,0),fout,#M"x"#N" x "#N"x1",prf_maccycle, (N*M));

static void mips_matop_mpy(int isFull, int isVerbose, FILE * fout)
{
  void* pScr = (void*)mips.scratch0.i32;
  PROFILE_mpy16x16(     1, isVerbose, 16,16,16);
  PROFILE_mpy16x16(     1, isVerbose, 32,32,32);
  PROFILE_mpy16x16(isFull, isVerbose, 40,80,8);
  PROFILE_mpy16x16(isFull, isVerbose, 40,81,8);
  PROFILE_mpy16x16(isFull, isVerbose, 40,82,8);
  PROFILE_mpy16x16(isFull, isVerbose, 40,83,8);
  PROFILE_mpy16x16(isFull, isVerbose, 40,84,8);
  PROFILE_mpy16x16(isFull, isVerbose, 40,85,8);
  PROFILE_mpy16x16(isFull, isVerbose, 40,86,8);
  PROFILE_mpy16x16(isFull, isVerbose, 40,87,8);
  PROFILE_mpy16x16(isFull, isVerbose, 40,88,8);
  PROFILE_mpy16x16(isFull, isVerbose, 2,100,8);
  PROFILE_mpy16x16(isFull, isVerbose, 8,80 ,2);
  PROFILE_mpy16x16(isFull, isVerbose, 8,4  ,2);
  PROFILE_mpy16x16(isFull, isVerbose, 8,16 ,2);
  PROFILE_mpy16x16(isFull, isVerbose, 8,32 ,2);
  PROFILE_mpy16x16_fast(     1, isVerbose, 16,16,16);
  PROFILE_mpy16x16_fast(     1, isVerbose, 32,32,32);
  PROFILE_mpy16x16_fast(     1, isVerbose, 8,80,4);
  PROFILE_mpy16x16_fast(isFull, isVerbose, 8,84,4);
  PROFILE_mpy16x16_fast(isFull, isVerbose, 8,4 ,4);
  PROFILE_mpy16x16_fast(isFull, isVerbose, 8,16,4);
  PROFILE_mpy16x16_fast(isFull, isVerbose, 8,32,4);
  PROFILE_mpy24x24(isFull, isVerbose, 40,80,8);
  PROFILE_mpy24x24(isFull, isVerbose, 40,81,8);
  PROFILE_mpy24x24(isFull, isVerbose, 40,82,8);
  PROFILE_mpy24x24(isFull, isVerbose, 40,83,8);
  PROFILE_mpy24x24(isFull, isVerbose, 2,100,8);
  PROFILE_mpy24x24(isFull, isVerbose, 8,80,2);
  PROFILE_mpy24x24(isFull, isVerbose, 8,4,2);
  PROFILE_mpy24x24(isFull, isVerbose, 8,16,2);
  PROFILE_mpy24x24(isFull, isVerbose, 8,32,2);
  PROFILE_mpy24x24_fast(     1, isVerbose,8,80,4);
  PROFILE_mpy24x24_fast(isFull, isVerbose,8,84,4);
  PROFILE_mpy24x24_fast(isFull, isVerbose,8,4,4);
  PROFILE_mpy24x24_fast(isFull, isVerbose,8,16,4);
  PROFILE_mpy24x24_fast(isFull, isVerbose,8,32,4);

  PROFILE_mpy32x32(     1, isVerbose, 16,16,16);
  PROFILE_mpy32x32(     1, isVerbose, 32,32,32);
  PROFILE_mpy32x32(isFull, isVerbose,40,80,8);
  PROFILE_mpy32x32(isFull, isVerbose,40,81,8);
  PROFILE_mpy32x32(isFull, isVerbose,40,82,8);
  PROFILE_mpy32x32(isFull, isVerbose,40,83,8);
  PROFILE_mpy32x32(isFull, isVerbose,2,100,8);
  PROFILE_mpy32x32(isFull, isVerbose,8,80 ,2);
  PROFILE_mpy32x32(isFull, isVerbose,8,4  ,2);
  PROFILE_mpy32x32(isFull, isVerbose,8,16 ,2);
  PROFILE_mpy32x32(isFull, isVerbose,8,32 ,2);
  PROFILE_mpy32x32_fast(     1, isVerbose, 16,16,16);
  PROFILE_mpy32x32_fast(     1, isVerbose, 32,32,32);
  PROFILE_mpy32x32_fast(     1, isVerbose,8,80,4);
  PROFILE_mpy32x32_fast(isFull, isVerbose,8,84,4);
  PROFILE_mpy32x32_fast(isFull, isVerbose,8,4 ,4);
  PROFILE_mpy32x32_fast(isFull, isVerbose,8,16,4);
  PROFILE_mpy32x32_fast(isFull, isVerbose,8,32,4);
}

static void mips_matop_vecmpy(int isFull, int isVerbose, FILE * fout)
{
  PROFILE_mtx_vecmpy16x16(isFull, isVerbose, mtx_vecmpy16x16,     16,100);
  PROFILE_mtx_vecmpy16x16(isFull, isVerbose, mtx_vecmpy16x16,     16,104);
  PROFILE_mtx_vecmpy16x16(isFull, isVerbose, mtx_vecmpy16x16,     40,40);
  PROFILE_mtx_vecmpy16x16(isFull, isVerbose, mtx_vecmpy16x16_fast,16,100);
  PROFILE_mtx_vecmpy16x16(     1, isVerbose, mtx_vecmpy16x16_fast,16,104);
  PROFILE_mtx_vecmpy16x16(isFull, isVerbose, mtx_vecmpy16x16_fast,40,40);
  PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24,     16,100);
  PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24,     16,101);
  PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24,     16,102);
  PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24,     16,103);
  PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24,     16,104);
  PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24,     40,40);
  PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24_fast,16,100);
  PROFILE_mtx_vecmpy24x24(     1, isVerbose, mtx_vecmpy24x24_fast,16,104);
  PROFILE_mtx_vecmpy24x24(isFull, isVerbose, mtx_vecmpy24x24_fast,40,40);
  PROFILE_mtx_vecmpy32x32(isFull, isVerbose,mtx_vecmpy32x32,     16,100);
  PROFILE_mtx_vecmpy32x32(isFull, isVerbose,mtx_vecmpy32x32,     16,101);
  PROFILE_mtx_vecmpy32x32(isFull, isVerbose,mtx_vecmpy32x32,     16,102);
  PROFILE_mtx_vecmpy32x32(isFull, isVerbose,mtx_vecmpy32x32,     16,103);
  PROFILE_mtx_vecmpy32x32(isFull, isVerbose,mtx_vecmpy32x32,     16,104);
  PROFILE_mtx_vecmpy32x32(isFull, isVerbose,mtx_vecmpy32x32,     40,40);
  PROFILE_mtx_vecmpy32x32(isFull, isVerbose,mtx_vecmpy32x32_fast,16,100);
  PROFILE_mtx_vecmpy32x32(     1, isVerbose,mtx_vecmpy32x32_fast,16,104);
  PROFILE_mtx_vecmpy32x32(isFull, isVerbose,mtx_vecmpy32x32_fast,40,40);
}


void mips_matop1(int isFull, int isVerbose, FILE * fout)
{
  mips_matop_mpy(isFull, isVerbose, fout);
  mips_matop_vecmpy(isFull, isVerbose, fout);
}
