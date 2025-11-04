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
* Test module for testing cycle performance (Correlation, Convolution, 
* Dispreading, LMS)
*/

#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(fir)
#include "mips.h"

void mips_convol(int isFull, int isVerbose, FILE * fout)
{
    void* pScr = (void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose,fir_convol16x16,(mips.out0.i16, mips.inp2.i16, mips.inp1.i16,    80, 56),fout,"N: 80; M: 56" ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convol16x16,(mips.out0.i16, mips.inp2.i16, mips.inp1.i16,    80, 60),fout,"N: 80; M: 60" ,prf_maccycle,    80*60 );
    PROFILE_INVERTED(     1, isVerbose,fir_convol16x16,(mips.out0.i16, mips.inp2.i16, mips.inp1.i16,   256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convol16x16,(mips.out0.i16, mips.inp2.i16, mips.inp1.i16,   256, 84),fout,"N: 256; M: 84",prf_maccycle,   84*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convol32x16,(mips.out0.i32, mips.inp2.i32, mips.inp1.i16+8,  80, 56),fout,"N: 80; M: 56" ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convol32x16,(mips.out0.i32, mips.inp2.i32, mips.inp1.i16  ,  80, 60),fout,"N: 80; M: 60" ,prf_maccycle,    80*60 );
    PROFILE_INVERTED(     1, isVerbose,fir_convol32x16,(mips.out0.i32, mips.inp2.i32, mips.inp1.i16+8, 256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convol32x16,(mips.out0.i32, mips.inp2.i32, mips.inp1.i16  , 256, 84),fout,"N: 256; M: 84",prf_maccycle,   84*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convol24x24,(mips.out0.i32, mips.inp2.i32, mips.inp1.i32,    80, 56),fout,"N: 80; M: 56" ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convol24x24,(mips.out0.i32, mips.inp2.i32, mips.inp1.i32,   256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );

    PROFILE_INVERTED(isFull, isVerbose,fir_convol32x32,(mips.out0.i32, mips.inp2.i32, mips.inp1.i32,    80, 56),fout,"N: 80; M: 56" ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convol32x32,(mips.out0.i32, mips.inp2.i32, mips.inp1.i32,   256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );

    PROFILE_INVERTED(isFull, isVerbose,fir_convola16x16,(pScr, mips.out0.i16, mips.inp2.i16, mips.inp1.i16,    80, 56),fout,"N=80; M=56"   ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convola16x16,(pScr, mips.out0.i16, mips.inp2.i16, mips.inp1.i16,   256, 80),fout,"N=256; M=80"  ,prf_maccycle,   80*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convola32x16,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16+4,  80, 56),fout,"N: 80; M: 56" ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convola32x16,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16+4,  80, 60),fout,"N: 80; M: 60" ,prf_maccycle,    80*60 );
    PROFILE_INVERTED(     1, isVerbose,fir_convola32x16,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16+4, 256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convola24x24,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32,    80, 56),fout,"N: 80; M: 56" ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convola24x24,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32,   256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_convola32x16,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16+4, 256, 84),fout,"N: 256; M: 84",prf_maccycle,   84*256 );

    PROFILE_INVERTED(isFull, isVerbose,fir_convola32x32,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32,    80, 56),fout,"N=80; M=56"   ,prf_maccycle,    80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convola32x32,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32,   256, 80),fout,"N=256; M=80"  ,prf_maccycle,   80*256 );

    PROFILE_INVERTED(isFull, isVerbose, cxfir_convol32x16, (mips.out0.ci32, mips.inp2.ci32, mips.inp1.ci16 + 4, 80, 56), fout, "N: 80; M: 56", prf_maccycle, 4 * 80 * 56);
    PROFILE_INVERTED(     1, isVerbose,cxfir_convol32x16, (mips.out0.ci32, mips.inp2.ci32, mips.inp1.ci16 + 4, 256, 80), fout, "N: 256; M: 80",prf_maccycle, 4*80 * 256);
    PROFILE_INVERTED(isFull, isVerbose, cxfir_convola32x16, (pScr, mips.out0.ci32, mips.inp2.ci32, mips.inp1.ci16 + 2, 80, 56), fout, "N: 80; M: 56", prf_maccycle, 4 * 80 * 56);
    PROFILE_INVERTED(     1, isVerbose,cxfir_convola32x16, (pScr, mips.out0.ci32, mips.inp2.ci32, mips.inp1.ci16 + 2, 256, 80), fout, "N: 256; M: 80",prf_maccycle, 4*80 * 256);
}

void mips_lconvol(int isFull, int isVerbose, FILE * fout)
{
    void* pScr=(void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose,fir_lconvola16x16,(pScr, mips.out0.i16, mips.inp2.i16, mips.inp1.i16,    80, 56),fout,"N=80; M=56",prf_maccycle,     80*56);
    PROFILE_INVERTED(     1, isVerbose,fir_lconvola16x16,(pScr, mips.out0.i16, mips.inp2.i16, mips.inp1.i16,   256, 80),fout,"N=256; M=80",prf_maccycle,   80*256);
    PROFILE_INVERTED(isFull, isVerbose,fir_lconvola32x32,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32,    80, 56),fout,"N=80; M=56",prf_maccycle,     80*56);
    PROFILE_INVERTED(     1, isVerbose,fir_lconvola32x32,(pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32,   256, 80),fout,"N=256; M=80",prf_maccycle,   80*256);
}

void mips_xcorr(int isFull, int isVerbose, FILE * fout)
{
    void* pScr = (void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose, fir_xcorr16x16, (mips.out0.i16, mips.inp2.i16, mips.inp1.i16, 80, 56  ),fout, "N: 80; M: 56" , prf_maccycle, 80 * 56);
    PROFILE_INVERTED(     1, isVerbose, fir_xcorr16x16, (mips.out0.i16, mips.inp2.i16, mips.inp1.i16, 256, 80 ),fout, "N: 256; M: 80", prf_maccycle, 80 * 256);
    PROFILE_INVERTED(isFull, isVerbose, fir_xcorr32x16, (mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 80, 56  ),fout, "N: 80; M: 56" , prf_maccycle, 80 * 56);
    PROFILE_INVERTED(isFull, isVerbose, fir_xcorr32x16, (mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 80, 60  ),fout, "N: 80; M: 60" , prf_maccycle, 80 * 60);
    PROFILE_INVERTED(     1, isVerbose, fir_xcorr32x16, (mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 256, 80 ),fout, "N: 256; M: 80", prf_maccycle, 80 * 256);
    PROFILE_INVERTED(isFull, isVerbose, fir_xcorr32x16, (mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 256, 84 ),fout, "N: 256; M: 84", prf_maccycle, 84 * 256);
    PROFILE_INVERTED(isFull, isVerbose, fir_xcorr24x24, (mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 80, 56  ),fout, "N: 80; M: 56" , prf_maccycle, 80 * 56);
    PROFILE_INVERTED(     1, isVerbose, fir_xcorr24x24, (mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 256, 80 ),fout, "N: 256; M: 80", prf_maccycle, 80 * 256);

    PROFILE_INVERTED(isFull, isVerbose,fir_xcorr32x32,  (mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 80, 56  ),fout, "N: 80; M: 56" , prf_maccycle, 80 * 56);
    PROFILE_INVERTED(     1, isVerbose,fir_xcorr32x32,  (mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 256, 80 ),fout, "N: 256; M: 80", prf_maccycle, 80 * 256);

    PROFILE_INVERTED(isFull, isVerbose,fir_xcorra16x16,  (pScr, mips.out0.i16, mips.inp2.i16, mips.inp1.i16, 80, 56  ),fout,"N: 80; M: 56" ,prf_maccycle, 80 * 56);
    PROFILE_INVERTED(     1, isVerbose,fir_xcorra16x16,  (pScr, mips.out0.i16, mips.inp2.i16, mips.inp1.i16, 256, 80 ),fout,"N: 256; M: 80",prf_maccycle, 80 * 256);
    PROFILE_INVERTED(isFull, isVerbose, fir_xcorra32x16, (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 80, 56  ),fout,"N: 80; M: 56" ,prf_maccycle, 80 * 56);
    PROFILE_INVERTED(isFull, isVerbose, fir_xcorra32x16, (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 80, 60  ),fout,"N: 80; M: 60" ,prf_maccycle, 80 * 60);
    PROFILE_INVERTED(     1, isVerbose, fir_xcorra32x16, (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 256, 80 ),fout,"N: 256; M: 80",prf_maccycle, 80 * 256);
    PROFILE_INVERTED(isFull, isVerbose, fir_xcorra32x16, (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i16, 256, 84 ),fout,"N: 256; M: 84",prf_maccycle, 84 * 256);
    PROFILE_INVERTED(isFull, isVerbose, fir_xcorra24x24, (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 80, 56  ),fout,"N: 80; M: 56" ,prf_maccycle, 80 * 56);
    PROFILE_INVERTED(     1, isVerbose, fir_xcorra24x24, (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 256, 80 ),fout,"N: 256; M: 80",prf_maccycle, 80 * 256);

    PROFILE_INVERTED(isFull, isVerbose,fir_xcorra32x32,  (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 80, 56  ),fout,"N: 80; M: 56" ,prf_maccycle, 80 * 56);
    PROFILE_INVERTED(     1, isVerbose,fir_xcorra32x32,  (pScr, mips.out0.i32, mips.inp2.i32, mips.inp1.i32, 256, 80 ),fout,"N: 256; M: 80",prf_maccycle, 80 * 256);
}

void mips_lxcorr(int isFull, int isVerbose, FILE * fout)
{
    void* pScr=(void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose, fir_lxcorra16x16, (pScr, mips.out0.i16,  mips.inp2.i16,  mips.inp1.i16, 80, 56), fout, "N=80; M=56", prf_maccycle, 80 * 56);
    PROFILE_INVERTED(     1, isVerbose,fir_lxcorra16x16,  (pScr, mips.out0.i16,  mips.inp2.i16,  mips.inp1.i16, 256, 80   ),fout,"N=256; M=80",prf_maccycle,   256 * 80);
    PROFILE_INVERTED(isFull, isVerbose,fir_lxcorra32x32,  (pScr, mips.out0.i32,  mips.inp2.i32,  mips.inp1.i32, 80, 56  ),fout,"N=80; M=56" ,prf_maccycle,    80*56);
    PROFILE_INVERTED(     1, isVerbose,fir_lxcorra32x32,  (pScr, mips.out0.i32,  mips.inp2.i32,  mips.inp1.i32, 256, 80 ),fout,"N=256; M=80",prf_maccycle,   80*256);
}

void mips_acorr(int isFull, int isVerbose, FILE * fout)
{
    void* pScr = (void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose,fir_acorr16x16,(mips.out0.i16, mips.inp2.i16, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorr16x16,(mips.out0.i16, mips.inp2.i16, 256),fout,"N: 256",prf_maccycle,  256*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_acorr24x24,(mips.out0.i32, mips.inp2.i32, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorr24x24,(mips.out0.i32, mips.inp2.i32, 256),fout,"N: 256",prf_maccycle,  256*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_acorr32x32,(mips.out0.i32, mips.inp2.i32, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorr32x32,(mips.out0.i32, mips.inp2.i32, 256),fout,"N: 256",prf_maccycle,  256*256 );

    PROFILE_INVERTED(isFull, isVerbose,fir_acorra16x16,(pScr, mips.out0.i16, mips.inp2.i16, 80),fout,"N=80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorra16x16,(pScr, mips.out0.i16, mips.inp2.i16, 256),fout,"N=256",prf_maccycle,  256*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_acorra24x24,(pScr, mips.out0.i32, mips.inp2.i32, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorra24x24,(pScr, mips.out0.i32, mips.inp2.i32, 256),fout,"N: 256",prf_maccycle,  256*256 );

    PROFILE_INVERTED(isFull, isVerbose,fir_acorra32x32,(pScr, mips.out0.i32, mips.inp2.i32, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorra32x32,(pScr, mips.out0.i32, mips.inp2.i32, 256),fout,"N: 256",prf_maccycle,  256*256 );
}

void mips_lacorr(int isFull, int isVerbose, FILE * fout)
{
    void* pScr=(void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose,fir_lacorra16x16,(pScr, mips.out0.i16, mips.inp2.i16, 80),fout,"N=80",prf_maccycle,    80*80  /2 );
    PROFILE_INVERTED(     1, isVerbose,fir_lacorra16x16,(pScr, mips.out0.i16, mips.inp2.i16, 256),fout,"N=256",prf_maccycle,  256*256/2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_lacorra32x32,(pScr, mips.out0.i32, mips.inp2.i32, 80),fout,"N=80",prf_maccycle,    80*80   /2);
    PROFILE_INVERTED(     1, isVerbose,fir_lacorra32x32,(pScr, mips.out0.i32, mips.inp2.i32, 256),fout,"N=256",prf_maccycle,  256*256 /2);
}

void mips_blms(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x16,(mips.out1.i16, mips.out2.i16, mips.inp1.i16, mips.inp0.i16, 0x111, 111, 80, 16),fout,"N: 80; M: 16",prf_maccycle,    80*16*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x16,(mips.out1.i16, mips.out2.i16, mips.inp1.i16, mips.inp0.i16, 0x111, 111, 64, 16),fout,"N: 64; M: 16",prf_maccycle,    64*16*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x16,(mips.out1.i16, mips.out2.i16, mips.inp1.i16, mips.inp0.i16, 0x111, 111, 64, 64),fout,"N: 64; M: 64",prf_maccycle,    64*64*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x16,(mips.out1.i16, mips.out2.i16, mips.inp1.i16, mips.inp0.i16, 0x111, 111, 80, 64),fout,"N: 80; M: 64",prf_maccycle,    80*64*2 );
    PROFILE_INVERTED(     1, isVerbose,fir_blms16x16,(mips.out1.i16, mips.out2.i16, mips.inp2.i16, mips.inp0.i16, 0x111, 111, 80, 128),fout,"N: 80; M: 128",prf_maccycle, 80*128*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x16,(mips.out1.i16, mips.out2.i16, mips.inp2.i16, mips.inp0.i16, 0x111, 111, 64, 128),fout,"N: 64; M: 128",prf_maccycle, 64*128*2 );

    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x32,(mips.out1.i32, mips.out2.i32, mips.inp1.i16, mips.inp0.i16, 0x111111, 111, 80, 16), fout,"N: 80; M: 16", prf_maccycle, 80*16 *2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x32,(mips.out1.i32, mips.out2.i32, mips.inp1.i16, mips.inp0.i16, 0x111111, 111, 64, 16), fout,"N: 64; M: 16", prf_maccycle, 64*16 *2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x32,(mips.out1.i32, mips.out2.i32, mips.inp1.i16, mips.inp0.i16, 0x111111, 111, 64, 64), fout,"N: 64; M: 64", prf_maccycle, 64*64 *2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x32,(mips.out1.i32, mips.out2.i32, mips.inp1.i16, mips.inp0.i16, 0x111111, 111, 80, 64), fout,"N: 80; M: 64", prf_maccycle, 80*64 *2 );
    PROFILE_INVERTED(     1, isVerbose,fir_blms16x32,(mips.out1.i32, mips.out2.i32, mips.inp2.i16, mips.inp0.i16, 0x111111, 111, 80, 128),fout,"N: 80; M: 128",prf_maccycle, 80*128*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms16x32,(mips.out1.i32, mips.out2.i32, mips.inp2.i16, mips.inp0.i16, 0x111111, 111, 64, 128),fout,"N: 64; M: 128",prf_maccycle, 64*128*2 );

    PROFILE_INVERTED(isFull, isVerbose,fir_blms24x24,(mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 80, 16),fout,"N: 80; M: 16",prf_maccycle,    80*16*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms24x24,(mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 64, 16),fout,"N: 64; M: 16",prf_maccycle,    64*16*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms24x24,(mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 64, 64),fout,"N: 64; M: 64",prf_maccycle,    64*64*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms24x24,(mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 80, 64),fout,"N: 80; M: 64",prf_maccycle,    80*64*2 );
    PROFILE_INVERTED(     1, isVerbose,fir_blms24x24,(mips.out1.i32, mips.out2.i32, mips.inp2.i32, mips.inp0.i32, 0x111111, 111, 80, 128),fout,"N: 80; M: 128",prf_maccycle, 80*128*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blms24x24,(mips.out1.i32, mips.out2.i32, mips.inp2.i32, mips.inp0.i32, 0x111111, 111, 64, 128),fout,"N: 64; M: 128",prf_maccycle, 64*128*2 );

    PROFILE_INVERTED(isFull, isVerbose,fir_blms32x32, (mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 80, 16), fout, "N: 80; M: 16", prf_maccycle, 80 * 16 * 2);
    PROFILE_INVERTED(isFull, isVerbose,fir_blms32x32, (mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 64, 16), fout, "N: 64; M: 16", prf_maccycle, 64 * 16 * 2);
    PROFILE_INVERTED(isFull, isVerbose,fir_blms32x32, (mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 64, 64), fout, "N: 64; M: 64", prf_maccycle, 64 * 64 * 2);
    PROFILE_INVERTED(isFull, isVerbose,fir_blms32x32, (mips.out1.i32, mips.out2.i32, mips.inp1.i32, mips.inp0.i32, 0x111111, 111, 80, 64), fout, "N: 80; M: 64", prf_maccycle, 80 * 64 * 2);
    PROFILE_INVERTED(     1, isVerbose,fir_blms32x32, (mips.out1.i32, mips.out2.i32, mips.inp2.i32, mips.inp0.i32, 0x111111, 111, 80, 128), fout, "N: 80; M: 128", prf_maccycle, 80 * 128 * 2);
    PROFILE_INVERTED(isFull, isVerbose,fir_blms32x32, (mips.out1.i32, mips.out2.i32, mips.inp2.i32, mips.inp0.i32, 0x111111, 111, 64, 128), fout, "N: 64; M: 128", prf_maccycle, 64 * 128 * 2);
}

void mips_firother1(int isFull, int isVerbose, FILE * fout)
{
    mips_convol(isFull, isVerbose, fout);
    mips_lconvol(isFull, isVerbose, fout);
    mips_xcorr(isFull, isVerbose, fout);
    mips_lxcorr(isFull, isVerbose, fout);
    mips_acorr(isFull, isVerbose, fout);
    mips_lacorr(isFull, isVerbose, fout);
    mips_blms(isFull, isVerbose, fout);
}
