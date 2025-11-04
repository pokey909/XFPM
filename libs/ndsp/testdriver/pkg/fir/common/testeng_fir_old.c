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
 * Test-engine add-on for older FIR API
 */

#include <string.h>
#include <stdlib.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"
#include "testeng_fir.h"
#include "testeng_fir_old.h"
/* Aligning memory allocator. */
#include "malloc16.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

typedef struct
{
    tVec ir;          /* impulse response */
    int  L;           /* number of streams */
    int  D;           /* interpolation/decimation ratio */
    void *firMem;     /* memory allocated for FIR */ 
    void *handle;     /* filter handle */
    void *extIrMem;   /* memory allocated for external IR */
}
tFirOldContext;


/* function reads impulse response from file and creates FIR structure. returns 0 if failed */
int te_create_fir_old(tTestEngContext * context)
{
    int extIR=(context->desc->extraParam & TE_FIR_EXTIR)?1:0;
    void *pExtIr;
    int res;
    size_t szObj=0,szExtIr;
    int M, L, D, fmtIr=0, irLen=0;
    tFirOldContext *firContext;
    if (seqFileScanf(context->seqFile, "%d %d %d ", &M, &L, &D) != 3)
    {
        printf("bad SEQ-file format\n");
        return 0;
    }
    firContext = (tFirOldContext *)malloc(sizeof(tFirOldContext));
    context->target.handle = (void*)firContext;
    if (firContext==NULL) return 0;
    memset(firContext, 0, sizeof(*firContext));
    firContext->D=D;
    firContext->L=L;
    switch(context->desc->extraParam & TE_FIR_FILTER_TYPE_MASK)
    {
        case TE_FIR_FIR: irLen=M;   fmtIr=context->desc->fmt;                                 break;
        case TE_FIR_DN:  irLen=M;   fmtIr=context->desc->fmt;                                 break;
        case TE_FIR_UP:  irLen=M*D; fmtIr=context->desc->fmt;                                 break;
        case TE_FIR_PP:  irLen=M*D; fmtIr = FMT_REAL | (context->desc->fmt & FMT_DTYPE_MASK); break;
        default: ASSERT("wrong extraParam");
    }
    if ( context->desc->extraParam & TE_FIR_FILTER_32X16 )
    {
        fmtIr = (0==(context->desc->fmt & FMT_CPLX)) ? FMT_REAL|FMT_FRACT16 : FMT_CPLX|FMT_FRACT16;
    }
    if (!vecAlloc(&firContext->ir, irLen, context->desc->isAligned, fmtIr, NULL)) return 0;
    if (!IS_PRESENT(((tFirOldDescr*)context->target.fut)->g.alloc) ||
        !IS_PRESENT((((tFirOldDescr*)context->target.fut)->g.init)))
    {
        firContext->firMem=malloc(1);    // FUT is not defined
        firContext->extIrMem = mallocAlign(1,0);
        return -1; 
    }
    szExtIr = (extIR) ? ((tFirOldDescr*)context->target.fut)->e.allocExtIr(M) : 0;
    switch(context->desc->extraParam & TE_FIR_FILTER_API_MASK)
    {
        case TE_FIR_OLDGENERIC:   /* generic filters */
            szObj=((tFirOldDescr*)context->target.fut)->g.alloc(M); break;
        case TE_FIR_OLDEXTIR:     /* filters with external IR */
            szObj=((tFirOldDescr*)context->target.fut)->e.alloc(M,extIR); 
            break;
        case TE_FIR_OLDDECIMATOR: /* decimator/interpolator */
        case TE_FIR_POLYPHASE:    /* polyphase filter */
            szObj=((tFirOldDescr*)context->target.fut)->d.alloc(D,M); break;
        default: 
            ASSERT("wrong extraParam");
    }

    firContext->firMem=malloc(szObj);
    res=seqFileReadVec(context->seqFile, &firContext->ir);
    /* allocate memory for extIR if any */
    pExtIr = firContext->extIrMem = mallocAlign(szExtIr, 0);
    if (extIR)
    {
        ((tFirOldDescr*)context->target.fut)->e.copyExtIr(pExtIr,vecGetElem(&firContext->ir,0),M);
    }
    switch(context->desc->extraParam & TE_FIR_FILTER_API_MASK)
    {
        case TE_FIR_OLDGENERIC:   /* generic filters */
            firContext->handle=((tFirOldDescr*)context->target.fut)->g.init(firContext->firMem,M,vecGetElem(&firContext->ir, 0));
            break;
        case TE_FIR_OLDEXTIR:     /* filters with external IR */
            firContext->handle=((tFirOldDescr*)context->target.fut)->e.init(firContext->firMem,M,extIR,
                extIR ? pExtIr : vecGetElem(&firContext->ir, 0));
            break;
        case TE_FIR_OLDDECIMATOR: /* decimator/interpolator */
        case TE_FIR_POLYPHASE:    /* polyphase filter */
            firContext->handle=(((tFirOldDescr*)context->target.fut)->d.init)(firContext->firMem,D,M,vecGetElem(&firContext->ir, 0));
            break;
    }

    return res;
}

/* function destroys FIR structure, returns 0 if failed */
int te_destroy_fir_old(tTestEngContext * context)
{
    tFirOldContext *firContext;
    firContext = (tFirOldContext *)context->target.handle;
    if (firContext)
    {
        ASSERT(firContext->firMem);
        vecFree(&firContext->ir);
        free(firContext->firMem);
        freeAlign(firContext->extIrMem);
        free(firContext);
    }
    return 1;
}

/* 
   Allocate vectors and load the data set for FIR:
*  vector X (in), vector Z (out) */
int te_loadFxn_fir_old(tTestEngContext * context)
{
    int fmtX;
    int M, N, L,D;
    int nElemIn=0,nElemOut=0, res;
    tFirOldContext *firContext;
    firContext = (tFirOldContext *)context->target.handle;

    ASSERT(context && context->seqFile);

    M = context->args.dim[0];
    N = context->args.dim[1];
    L = firContext->L;
    D = firContext->D;

    nElemIn = MAX(0, M*N*L);
    switch(context->desc->extraParam & TE_FIR_FILTER_TYPE_MASK)
    {
        case TE_FIR_FIR: nElemIn = nElemOut = MAX(0, M*N*L)          ; break;
        case TE_FIR_DN:  nElemIn = MAX(0, M*N*L); nElemOut =nElemIn/D; break;
        case TE_FIR_UP:  nElemIn = MAX(0, M*N*L); nElemOut =nElemIn*D; break;
        case TE_FIR_PP:  nElemIn = nElemOut = MAX(0, M*L);             break;
        default: ASSERT("wrong extraParam");
    }

    memset(&context->dataSet, 0, sizeof(context->dataSet));
    fmtX=context->desc->fmt;

    /* Allocate data vectors memory. */
    res = vecAlloc(&context->dataSet.X, nElemIn, context->desc->isAligned, fmtX, NULL);
    if (context->desc->extraParam & TE_FIR_LSH)
    {
        res &= vecAlloc(&context->dataSet.Y, 1, TE_ALIGN_NO, FMT_INT32, NULL);
    }
    res &= (3 == vecsAlloc(0, fmtX,
                            &context->dataSet.Z, nElemOut,
                            &context->dataSet.Zlo, nElemOut,
                            &context->dataSet.Zhi, nElemOut, 0));
    if (res)
    {
        /* Load vectors data from the SEQ-file. */
        if (context->desc->extraParam & TE_FIR_LSH)
        {
            res = seqFileReadVecs(context->seqFile,
                                  &context->dataSet.X,
                                  &context->dataSet.Y,
                                  &context->dataSet.Zlo,
                                  &context->dataSet.Zhi, 0);
        }
        else
        {
            res = seqFileReadVecs(context->seqFile,
                                  &context->dataSet.X,
                                  &context->dataSet.Zlo,
                                  &context->dataSet.Zhi, 0);
        }

        if (!res)
        {
            printf("te_loadFxn_fir_old(): failed to read vectors data; "
                    "fmt = 0x%02x, nElemIn = %d, nElemOut = %d\n",
                    (unsigned)context->desc->fmt, nElemIn, nElemOut);
        }
    }
    else
    {
        printf("te_loadFxn_fir_old(): failed to allocate vectors; "
                "fmt = 0x%02x, nElemIn = %d, nElemOut = %d\n",
                (unsigned)context->desc->fmt, nElemIn, nElemOut);
    }

    /* Free vectors data if failed. */
    if (!res) te_freeVectors(context);
    return (res);

} /* te_loadFxn_fir_old() */


/* Apply FIR function to the test case data set.
*  vector X (in), vector Z (out) */
void te_processFxn_fir_old(tTestEngContext * context)
{
    tFirOldFxnProcess *fxn;
    void *X, *Z;
    tFirOldContext* firContext;
    int D,N;

    ASSERT(context && context->target.fut);

    firContext = (tFirOldContext *)context->target.handle;
    X = vecGetElem(&context->dataSet.X, 0);
    Z = vecGetElem(&context->dataSet.Z, 0);

    N = context->args.dim[0];
    D = firContext->D;
    switch(context->desc->extraParam & TE_FIR_FILTER_TYPE_MASK)
    {
        case TE_FIR_FIR: N=N*1;   break;
        case TE_FIR_DN:  N=N/D;   break;
        case TE_FIR_UP:  N=N*1;   break;
        case TE_FIR_PP:  N=N*1;   break;
        default: ASSERT("wrong extraParam");
    }

    {   /* logging */
        const tFirOldDescr* descr=(const tFirOldDescr*)context->target.fut;
        tReportFUT fut[3];
        fut[0]=(tReportFUT)descr->g.process;
        fut[1]=(tReportFUT)descr->g.alloc;
        fut[2]=(tReportFUT)descr->g.init;
        vReportAdd( fut, 3, "", context->seqFile->filename, context->args.caseType, te_vGetDataSize(context));
    }

    fxn = ((tFirOldDescr*)context->target.fut)->g.process;
    if ((context->desc->extraParam & TE_FIR_FILTER_API_MASK) == TE_FIR_POLYPHASE)
    {
        if (context->desc->extraParam & TE_FIR_LSH)
        {
            int lsh = *vecGetElem_i32(&context->dataSet.Y, 0);
            ((tFirOldPolyphaseFxnProcess_lsh*)fxn)(firContext->handle, Z, X, lsh);
        }
        else
        {
            ((tFirOldPolyphaseFxnProcess*)fxn)(firContext->handle, Z, X);
        }
    }
    else
    {
        fxn(firContext->handle, Z, X, N);
    }
} /* te_processFxn_fir_old() */


/* management of external IR */
size_t bkfir16x16_allocExtIr(int M)
{
    int coefLen = M+4;
    return sizeof(int16_t)*(coefLen);
}
size_t bkfir32x16_allocExtIr(int M)
{
    int coefLen = (M>32) ? M+8:M+4;
    return sizeof(int16_t)*(coefLen);
}
size_t bkfir24x24_allocExtIr(int M)
{
    int coefLen = M+4;
    return sizeof(f24)*(coefLen);
}
size_t bkfir32x32_allocExtIr(int M)
{
    int coefLen = M+4;
    return sizeof(int32_t)*(coefLen);
}
size_t bkfir24x24p_allocExtIr(int M)
{
    int coefLen = M+(-M&4)+8;
    return sizeof(uint8_t)*(coefLen*3);
}
size_t cxfir16x16_allocExtIr(int M)
{
    int coefLen;
    if (!NatureDSP_Signal_get_isa_opt(NATUREDSP_ISA_OPT_HIFI3Z))
    {
        coefLen = M+2;
    }
    else
    {
        coefLen = 2 * (M + 4);
    }
    return sizeof(complex_fract16)*(coefLen);
}
size_t cxfir24x24_allocExtIr(int M)
{
    int coefLen = M+2;
    return 2*sizeof(f24)*(coefLen);
}
size_t cxfir32x16_allocExtIr(int M)
{
    int coefLen = M+2;
    return sizeof(complex_fract16)*(coefLen);
}
size_t cxfir32x32_allocExtIr(int M)
{
    int coefLen = M;
    return sizeof(complex_fract32)*(coefLen);
}
size_t bkfirf_allocExtIr(int M)
{
    int coefLen = M;
    return sizeof(float32_t)*(coefLen);
}
size_t cxfirf_allocExtIr(int M)
{
    int coefLen = M;
    return sizeof(complex_float)*(coefLen);
}


/* left padding 2 bytes, right padding 6 bytes, inverted order */
void bkfir16x16_copyExtIr(int16_t* pExtIr, const int16_t* ir, int M)
{   
    int n;
    *pExtIr++=0;
    for (n=0; n<M;  n++) *pExtIr++=ir[M-1-n];
    *pExtIr++=0;
    *pExtIr++=0;
    *pExtIr++=0;
}

/* left padding 2 or 10 bytes, right padding 6 bytes, inverted order */
void bkfir32x16_copyExtIr(int16_t* pExtIr, const int16_t* ir, int M)
{   
    int n;
    for (n=0; n<(M>32 ? 5:1); n++) *pExtIr++=0;
    for (n=0; n<M;  n++) *pExtIr++=ir[M-1-n];
    for (n=0; n<3; n++) *pExtIr++=0;
}

/* left padding 4 bytes, right padding 12 bytes, inverted order */
void bkfir24x24_copyExtIr(f24* pExtIr, const f24* ir, int M)
{   
    int n;
    *pExtIr++=0;
    for (n=0; n<M;  n++) *pExtIr++=ir[M-1-n];
    *pExtIr++=0;
    *pExtIr++=0;
    *pExtIr++=0;
}

/* left padding 4 bytes, right padding 12 bytes, inverted order */
void bkfir32x32_copyExtIr(int32_t* pExtIr, const int32_t* ir, int M)
{   
    int n;
    *pExtIr++=0;
    for (n=0; n<M;  n++) *pExtIr++=ir[M-1-n];
    *pExtIr++=0;
    *pExtIr++=0;
    *pExtIr++=0;
}
/* left padding ((-M&4)+1)*3 bytes, right padding 7, inverted order */
void bkfir24x24p_copyExtIr(uint8_t* pExtIr, const uint8_t* ir, int M)
{   
    int n,left,right;
    left = ((-M&4)+1);
    right= 7;
    for (n=0; n<left*3; n++)  *pExtIr++=0;
    for (n=0; n<M;  n++) 
    {
        *pExtIr++=ir[(M-1-n)*4+1];  /* note: ir in int32_t (4 bytes), pExtIr is packed f24 (3 bytes, least significant byte is omitted) */
        *pExtIr++=ir[(M-1-n)*4+2];
        *pExtIr++=ir[(M-1-n)*4+3];
    }
    for (n=0; n<right*3; n++)  *pExtIr++=0;
}

void cxfir16x16_copyExtIr(int16_t* pExtIr, const int16_t* ir, int M)
{   
    int n;
    /* code for no Quad MAC option */
    if (!NatureDSP_Signal_get_isa_opt(NATUREDSP_ISA_OPT_HIFI3Z))
    {
        *pExtIr++=0; *pExtIr++=0;
        for (n=0; n<M;  n++) 
        {
            *pExtIr++=ir[(M-1-n)*2+0];
            *pExtIr++=ir[(M-1-n)*2+1];
        }
        *pExtIr++=0; *pExtIr++=0;
    }
    else /* code for Quad MAC option (HiFi3z only) */
    {
        int16_t* pExtIrB;
        pExtIrB = pExtIr + 2 * (M + 4);
        *pExtIr++=0; *pExtIr++=0;
        *pExtIrB++=0; *pExtIrB++=0;
    
        for (n = 1; n<M+1; n++)
        {
        *pExtIr++= ir[(M-n)*2+0];
        *pExtIr++=-ir[(M-n)*2+1];

        *pExtIrB++=ir[(M-n)*2+1];
        *pExtIrB++=ir[(M-n)*2+0];
        }

        for (; n<M + 4; n++)
        {
            *pExtIr++=0; *pExtIr++=0;
            *pExtIrB++=0; *pExtIrB++=0;
        }
    }
}

void cxfir32x16_copyExtIr(int16_t* pExtIr, const int16_t* ir, int M)
{   
    int n;
    *pExtIr++=0; *pExtIr++=0;
    for (n=0; n<M;  n++) 
    {
        *pExtIr++=ir[(M-1-n)*2+0];
        *pExtIr++=ir[(M-1-n)*2+1];
    }
    *pExtIr++=0; *pExtIr++=0;
}

/* no padding, inverted order, complex coefficients */
void cxfir24x24_copyExtIr(f24* pExtIr, const f24* ir, int M)
{   
    int n;
    for (n=0; n<M;  n++) 
    {
        *pExtIr++=ir[(M-1-n)*2+0];
        *pExtIr++=ir[(M-1-n)*2+1];
    }
}

/* left padding 8 bytes, inverted order, complex coefficients */
void cxfir32x32_copyExtIr(int32_t* pExtIr, const int32_t* ir, int M)
{   
    int n;
    for (n=0; n<M;  n++) 
    {
        *pExtIr++= ir[(M-1-n)*2+0];
        *pExtIr++=-ir[(M-1-n)*2+1];
    }
}
/* no padding, inverted order*/
void cxfir32x32ep_copyExtIr(int32_t* pExtIr, const int32_t* ir, int M)
{   
    int n;
    for (n=0; n<M;  n++) 
    {
        int32_t t;
        t=(ir[(M-1-n)*2+1]);
        t=t==MIN_INT32 ? MAX_INT32:-t;
        *pExtIr++= ir[(M-1-n)*2+0];
        *pExtIr++= t;
    }
}

/* no padding, direct order */
void bkfirf_copyExtIr(float32_t* pExtIr, const float32_t* ir, int M)
{   
    int n;
    for (n=0; n<M;  n++) *pExtIr++=ir[n];

}

/* no padding, direct order */
void cxfirf_copyExtIr(complex_float* pExtIr, const complex_float* ir, int M)
{   
    int n;
    for (n=0; n<M;  n++) 
    {
        *pExtIr++= ir[n];
    }
}
