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
plug and play package API
list of all (not necessary available
package processing functions
*/
#ifndef __PACKAGES_H__
#define __PACKAGES_H__

#include <stdio.h>

// order of packages
#define PACKAGE_FIR     1
#define PACKAGE_IIR     2
#define PACKAGE_MATH    3
#define PACKAGE_COMPLEX 4
#define PACKAGE_VECTOR	5
#define PACKAGE_MATOP   7
#define PACKAGE_MATINV  8
#define PACKAGE_FIT     9
#define PACKAGE_FFT     10

#ifdef __cplusplus
extern "C"
{
#endif
    /*========================== FIR ============================*/
    int  func_firblk1(int isFull, int isVerbose, int breakOnError); 
    void mips_firblk1(int isFull, int isVerbose, FILE* f);
    int  func_firblk2(int isFull, int isVerbose, int breakOnError); 
    void mips_firblk2(int isFull, int isVerbose, FILE* f);

    int  func_firdec1(int isFull, int isVerbose, int breakOnError); 
    void mips_firdec1(int isFull, int isVerbose, FILE* f);
    int  func_firdec2(int isFull, int isVerbose, int breakOnError); 
    void mips_firdec2(int isFull, int isVerbose, FILE* f);

    int  func_firint1(int isFull, int isVerbose, int breakOnError); 
    void mips_firint1(int isFull, int isVerbose, FILE* f);
    int  func_firint2(int isFull, int isVerbose, int breakOnError); 
    void mips_firint2(int isFull, int isVerbose, FILE* f);
	
    int  func_firother1(int isFull, int isVerbose, int breakOnError); 
    void mips_firother1(int isFull, int isVerbose, FILE* f);
    int  func_firother2(int isFull, int isVerbose, int breakOnError); 
    void mips_firother2(int isFull, int isVerbose, FILE* f);

    /*========================== IRR ===========================*/
    int  func_iirbq1(int isFull, int isVerbose, int breakOnError); 
    void mips_iirbq1(int isFull, int isVerbose, FILE* f);
    int  func_iirbq2(int isFull, int isVerbose, int breakOnError); 
    void mips_iirbq2(int isFull, int isVerbose, FILE* f);

    int  func_iirlt1(int isFull, int isVerbose, int breakOnError); 
    void mips_iirlt1(int isFull, int isVerbose, FILE* f);
    int  func_iirlt2(int isFull, int isVerbose, int breakOnError); 
    void mips_iirlt2(int isFull, int isVerbose, FILE* f);

    /*========================= MATH ===========================*/
    int  func_mathv1(int isFull, int isVerbose, int breakOnError); 
    void mips_mathv1(int isFull, int isVerbose, FILE* f);
    int  func_mathv2(int isFull, int isVerbose, int breakOnError); 
    void mips_mathv2(int isFull, int isVerbose, FILE* f);
    int  func_mathvf1(int isFull, int isVerbose, int breakOnError); 
    void mips_mathvf1(int isFull, int isVerbose, FILE* f);

    int  func_maths1(int isFull, int isVerbose, int breakOnError);
    void mips_maths1(int isFull, int isVerbose, FILE* f);
    int  func_maths2(int isFull, int isVerbose, int breakOnError);
    void mips_maths2(int isFull, int isVerbose, FILE* f);
	int   acc_maths2(int isFull, int isVerbose, int breakOnError, int optAccuracy, int optException);

    /*========================== COMPLEX ==========================*/
    int  func_complexv1(int isFull, int isVerbose, int breakOnError);
    int  func_complexv2(int isFull, int isVerbose, int breakOnError);

    void mips_complexv1(int isFull, int isVerbose, FILE* f);
    void mips_complexv2(int isFull, int isVerbose, FILE* f);

    int  func_complexs2(int isFull, int isVerbose, int breakOnError);
    void mips_complexs2(int isFull, int isVerbose, FILE* f);

    /*========================= VECTOR ==========================*/
    int  func_vector1(int isFull, int isVerbose, int breakOnError);
    void mips_vector1(int isFull, int isVerbose, FILE* f);
    int  func_vector2(int isFull, int isVerbose, int breakOnError);
    void mips_vector2(int isFull, int isVerbose, FILE* f);

    /*========================= MATOP ==========================*/
    int  func_matop1(int isFull, int isVerbose, int breakOnError);
    void mips_matop1(int isFull, int isVerbose, FILE* f);
    int  func_matop2(int isFull, int isVerbose, int breakOnError);
    void mips_matop2(int isFull, int isVerbose, FILE* f);

    /*======================= MATINV ========================*/
    int  func_gj2(int isFull, int isVerbose, int breakOnError);
    void mips_gj2(int isFull, int isVerbose, FILE* f);

    /*========================= FIT ===========================*/
    int  func_pfit1(int isFull, int isVerbose, int breakOnError);
    void mips_pfit1(int isFull, int isVerbose, FILE* f);
    int  func_pfit2(int isFull, int isVerbose, int breakOnError);
    void mips_pfit2(int isFull, int isVerbose, FILE* f);

    /*========================= FFT ===========================*/
    int  func_cfft1(int isFull, int isVerbose, int breakOnError);
    void mips_cfft1(int isFull, int isVerbose, FILE* f);
    int  func_cfftie1(int isFull, int isVerbose, int breakOnError);
    void mips_cfftie1(int isFull, int isVerbose, FILE* f);
    int  func_cfftie2(int isFull, int isVerbose, int breakOnError);
    void mips_cfftie2(int isFull, int isVerbose, FILE* f);

    int  func_cnfft1(int isFull, int isVerbose, int breakOnError);
    void mips_cnfft1(int isFull, int isVerbose, FILE* f);

    int  func_rfft1(int isFull, int isVerbose, int breakOnError);
    void mips_rfft1(int isFull, int isVerbose, FILE* f);

    int  func_rfftie1(int isFull, int isVerbose, int breakOnError);
    void mips_rfftie1(int isFull, int isVerbose, FILE* f);
    int  func_rfftie2(int isFull, int isVerbose, int breakOnError);
    void mips_rfftie2(int isFull, int isVerbose, FILE* f);

    int  func_rnfft1(int isFull, int isVerbose, int breakOnError);
    void mips_rnfft1(int isFull, int isVerbose, FILE* f);

    int  func_dct1(int isFull, int isVerbose, int breakOnError);
    void mips_dct1(int isFull, int isVerbose, FILE* f);
    int  func_dct2(int isFull, int isVerbose, int breakOnError);
    void mips_dct2(int isFull, int isVerbose, FILE* f);

#ifdef __cplusplus
}
#endif

#endif /* __PACKAGES_H__ */
