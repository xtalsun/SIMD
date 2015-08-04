/***********************************************************************
* @file SIMD.h
* @author Shiheng Liang
* @date 2015-07-31
* @brief SIMD instructions used in vector arithmetic
***********************************************************************/
#ifndef SIMD_H
#define SIMD_H

#include <stdlib.h>     //size_t

namespace SIMD
{
    /**********************Notice:************************
     Assume buf[]=buf1[]={a1,a2,...,an}, buf2[]={b1,b2,...,bn},
     result[]={c1,c2,...,cn}, buflen = n
    *****************************************************/

    ///ci = ai + bi, i = 1 to n
    void _add(float buf1[], float buf2[], float result[], size_t buflen);
    void _add(double buf1[], double buf2[], double result[], size_t buflen);

    ///ci = ai - bi, i = 1 to n
    void _sub(float buf1[], float buf2[], float result[], size_t buflen);
    void _sub(double buf1[], double buf2[], double result[], size_t buflen);

    ///ci = ai * bi, i = 1 to n
    void _mul(float buf1[], float buf2[], float result[], size_t buflen);
    void _mul(double buf1[], double buf2[], double result[], size_t buflen);

    void _divide(float buf1[], float buf2[], float result[], size_t buflen);
    void _divide(double buf1[], double buf2[], double result[], size_t buflen);

    ///Round up, such as a1 = 5.4, after _ceil(buf), a1 = 6.
    void _ceil(float buf[], float result[], size_t buflen);
    void _ceil(double buf[], double result[], size_t buflen);
    ///Round down, such as a1 = 5.6, after _floor(buf), a1 = 5.
    void _floor(float buf[], float result[], size_t buflen);
    void _floor(double buf[], double result[], size_t buflen);

    ///Return max(a1,a2,...,an)
    float _max(float buf[], size_t buflen);
    double _max(double buf[], size_t buflen);
    ///Return min(a1,a2,...,an)
    float _min(float buf[], size_t buflen);
    double _min(double buf[], size_t buflen);

    ///ci = 1 / ai, i = 1 to n
    void _rcp(float buf[], float result[], size_t buflen);
    ///ci = 1 / sqrt(ai), i = 1 to n
    void _rsqrt(float buf[], float result[], size_t buflen);
    ///ci = sqrt(ai), i = 1 to n
    void _sqrt(float buf[], float result[], size_t buflen);
    void _sqrt(double buf[], double result[], size_t buflen);

    ///ci = ln(ai), i = 1 to n
    void _ln(float buf[], float result[], double ln_lookup[], size_t buflen);
    void _ln(double buf[], double result[], double ln_lookup[], size_t buflen);
}

#endif
