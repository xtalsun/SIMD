/***********************************************************************
* @file SIMD.cpp
* @author XtalSun
* @date 2015-07-31
* @brief SIMD instructions used in vector arithmetic
***********************************************************************/

#include <immintrin.h>      //AVX
#include <tmmintrin.h>      //SSSE3
#include <math.h>
#include "SIMD.h"
#include "fastlog.h"

#define BlockWidth_4 size_t BlockWidth = 4;\
                    size_t nBlocks = buflen / BlockWidth;\
                    size_t nRemain = buflen % BlockWidth;

#define BlockWidth_8 size_t BlockWidth = 8;\
                    size_t nBlocks = buflen / BlockWidth;\
                    size_t nRemain = buflen % BlockWidth;

#define PointerFloat_1 float *src = buf;\
                    float *dst = result;

#define PointerDouble_1 double *src = buf;\
                    double *dst = result;

#define PointerFloat_2 float *src1 = buf1;\
                    float *src2 = buf2;\
                    float *dst = result;

#define PointerDouble_2 double *src1 = buf1;\
                    double *src2 = buf2;\
                    double *dst = result;

#define PointerMove_1 src += BlockWidth;\
                    dst += BlockWidth;

#define PointerMove_2 src1 += BlockWidth;\
                    src2 += BlockWidth;\
                    dst += BlockWidth;

 void SIMD::_add(float buf1[], float buf2[], float result[], size_t buflen)
{
    BlockWidth_8;     //Declaration of BlockWidth, nBlocks, nRemain
    PointerFloat_2;     //Declaration of src1, src2, dst
    __m256 load1, load2;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load1 = _mm256_loadu_ps(src1);
        load2 = _mm256_loadu_ps(src2);
        load1 = _mm256_add_ps(load1, load2);
        _mm256_storeu_ps(dst, load1);
        PointerMove_2;
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        dst[i] = src1[i] + src2[i];
    }
}

void SIMD::_add(double buf1[], double buf2[], double result[], size_t buflen)
{
    BlockWidth_4;     //Declaration of BlockWidth, nBlocks, nRemain
    PointerDouble_2;    //Declaration of src1, src2, dst
    __m256d load1, load2;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load1 = _mm256_loadu_pd(src1);
        load2 = _mm256_loadu_pd(src2);
        load1 = _mm256_add_pd(load1, load2);
        _mm256_storeu_pd(dst, load1);
        PointerMove_2;
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        dst[i] = src1[i] + src2[i];
    }
}

void SIMD::_ceil(float buf[], float result[], size_t buflen)
{
    BlockWidth_8;       //Declaration of BlockWidth, nBlocks, nRemain
    PointerFloat_1;     //Declaration of src, dst
    __m256 load;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load = _mm256_loadu_ps(src);
        load = _mm256_ceil_ps(load);
        _mm256_storeu_ps(dst, load);
        PointerMove_1;
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        dst[i] = ceil(src[i]);
    }
}

void SIMD::_ceil(double buf[], double result[], size_t buflen)
{
    BlockWidth_4;       //Declaration of BlockWidth, nBlocks, nRemain
    PointerDouble_1;     //Declaration of src, dst
    __m256d load;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load = _mm256_loadu_pd(src);
        load = _mm256_ceil_pd(load);
        _mm256_storeu_pd(dst, load);
        PointerMove_1;
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        dst[i] = ceil(src[i]);
    }
}

void SIMD::_divide(float buf1[], float buf2[], float result[], size_t buflen)
{
    BlockWidth_8;
    PointerFloat_2;
    __m256 load1, load2;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load1 = _mm256_loadu_ps(src1);
        load2 = _mm256_loadu_ps(src2);
        load1 = _mm256_div_ps(load1, load2);
        _mm256_storeu_ps(dst, load1);
        PointerMove_2;
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        dst[i] = src1[i] / src2[i];
    }
}

void SIMD::_divide(double buf1[], double buf2[], double result[], size_t buflen)
{
    BlockWidth_4;
    PointerDouble_2;
    __m256d load1, load2;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load1 = _mm256_loadu_pd(src1);
        load2 = _mm256_loadu_pd(src2);
        load1 = _mm256_div_pd(load1, load2);
        _mm256_storeu_pd(dst, load1);
        PointerMove_2;
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        dst[i] = src1[i] / src2[i];
    }
}

void SIMD::_floor(float buf[], float result[], size_t buflen)
{
    BlockWidth_8;
    PointerFloat_1;
    __m256 load;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load = _mm256_loadu_ps(src);
        load = _mm256_floor_ps(load);
        _mm256_storeu_ps(dst, load);
        src += BlockWidth;
        dst += BlockWidth;
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        dst[i] = floor(src[i]);
    }
}

void SIMD::_floor(double buf[], double result[], size_t buflen)
{
    BlockWidth_4;
    PointerDouble_1;
    __m256d load;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load = _mm256_loadu_pd(src);
        load = _mm256_floor_pd(load);
        _mm256_storeu_pd(dst, load);
        src += BlockWidth;
        dst += BlockWidth;
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        dst[i] = floor(src[i]);
    }
}

float SIMD::_max(float buf[], size_t buflen)
{
    BlockWidth_8;
    float MAX, *tmp, *src = buf;
    __m256 load1, load2;
    load1 = _mm256_loadu_ps(src);
    load2 = _mm256_loadu_ps(src + BlockWidth);
    load1 = _mm256_max_ps(load1, load2);
    src += BlockWidth*2;
    for(size_t i = 2; i < nBlocks; i++)
    {
        load2 = _mm256_loadu_ps(src);
        load1 = _mm256_max_ps(load1, load2);
        src += BlockWidth;
    }
    tmp = (float *)&load1;
    MAX = tmp[0];
    for(size_t i = 1; i < BlockWidth; i++)
    {
        if(tmp[i] > MAX)
        {
            MAX = tmp[i];
        }
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        if(src[i] > MAX)
        {
            MAX = src[i];
        }
    }

    return MAX;
}

double SIMD::_max(double buf[], size_t buflen)
{
    BlockWidth_4;
    double MAX, *tmp, *src = buf;
    __m256d load1, load2;
    load1 = _mm256_loadu_pd(src);
    load2 = _mm256_loadu_pd(src + BlockWidth);
    load1 = _mm256_max_pd(load1, load2);
    src += BlockWidth*2;
    for(size_t i = 2; i < nBlocks; i++)
    {
        load2 = _mm256_loadu_pd(src);
        load1 = _mm256_max_pd(load1, load2);
        src += BlockWidth;
    }
    tmp = (double *)&load1;
    MAX = tmp[0];
    for(size_t i = 1; i < BlockWidth; i++)
    {
        if(tmp[i] > MAX)
        {
            MAX = tmp[i];
        }
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        if(src[i] > MAX)
        {
            MAX = src[i];
        }
    }

    return MAX;
}

__m256 recursive_max(float buf[], size_t BlockWidth, size_t head, size_t tail)
{
    if(tail - head <= 1)
    {
        __m256 load1, load2;
        float *src = buf;
        load1 = _mm256_loadu_ps(src + head*BlockWidth);
        load2 = _mm256_loadu_ps(src + tail*BlockWidth);
        load1 = _mm256_max_ps(load1, load2);
        return load1;
    }

    return _mm256_max_ps(recursive_max(buf, BlockWidth, head, (head+tail)/2),
                          recursive_max(buf, BlockWidth, (head+tail)/2 +1, tail));
}

float recur_max(float buf[], size_t buflen)
{
    size_t BlockWidth = 8;
    size_t nBlocks = buflen / BlockWidth;
    float MAX, *tmp, *src = buf;
    __m256 load;
    load = recursive_max(src, BlockWidth, 0, nBlocks-1);
    tmp = (float *)&load;
    MAX = tmp[0];
    for(size_t i = 1; i < BlockWidth; i++)
    {
        if(tmp[i] > MAX)
        {
            MAX = tmp[i];
        }
    }

    for(size_t i = BlockWidth *nBlocks; i < buflen; i++)
    {
        if(src[i] > MAX)
        {
            MAX = src[i];
        }
    }
    return MAX;
}

float nor_max(float buf[], size_t buflen)
{
    float MAX = buf[0];
    for(size_t i = 0; i < buflen; i++)
    {
        if(buf[i] > MAX)
        {
            MAX = buf[i];
        }
    }
    return MAX;
}

float SIMD::_min(float buf[], size_t buflen)
{
    BlockWidth_8;
    float MIN, *tmp, *src = buf;
    __m256 load1, load2;
    load1 = _mm256_loadu_ps(src);
    load2 = _mm256_loadu_ps(src + BlockWidth);
    load1 = _mm256_min_ps(load1, load2);
    src += BlockWidth*2;
    for(size_t i = 2; i < nBlocks; i++)
    {
        load2 = _mm256_loadu_ps(src);
        load1 = _mm256_min_ps(load1, load2);
        src += BlockWidth;
    }
    tmp = (float *)&load1;
    MIN = tmp[0];
    for(size_t i = 1; i < BlockWidth; i++)
    {
        if(tmp[i] < MIN)
        {
            MIN = tmp[i];
        }
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        if(src[i] < MIN)
        {
            MIN = src[i];
        }
    }

    return MIN;
}

double SIMD::_min(double buf[], size_t buflen)
{
    BlockWidth_4;
    double MIN, *tmp, *src = buf;
    __m256d load1, load2;
    load1 = _mm256_loadu_pd(src);
    load2 = _mm256_loadu_pd(src + BlockWidth);
    load1 = _mm256_min_pd(load1, load2);
    src += BlockWidth*2;
    for(size_t i = 2; i < nBlocks; i++)
    {
        load2 = _mm256_loadu_pd(src);
        load1 = _mm256_min_pd(load1, load2);
        src += BlockWidth;
    }
    tmp = (double *)&load1;
    MIN = tmp[0];
    for(size_t i = 1; i < BlockWidth; i++)
    {
        if(tmp[i] < MIN)
        {
            MIN = tmp[i];
        }
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        if(src[i] < MIN)
        {
            MIN = src[i];
        }
    }

    return MIN;
}

void SIMD::_mul(float buf1[], float buf2[], float result[], size_t buflen)
{
    BlockWidth_8;
    PointerFloat_2;
    __m256 load1, load2;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load1 = _mm256_loadu_ps(src1);
        load2 = _mm256_loadu_ps(src2);
        load1 = _mm256_mul_ps(load1, load2);
        _mm256_storeu_ps(dst, load1);
        PointerMove_2;
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        dst[i] = src1[i] * src2[i];
    }
}

void SIMD::_mul(double buf1[], double buf2[], double result[], size_t buflen)
{
    BlockWidth_4;
    PointerDouble_2;
    __m256d load1, load2;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load1 = _mm256_loadu_pd(src1);
        load2 = _mm256_loadu_pd(src2);
        load1 = _mm256_mul_pd(load1, load2);
        _mm256_storeu_pd(dst, load1);
        PointerMove_2;
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        dst[i] = src1[i] * src2[i];
    }
}

void SIMD::_rcp(float buf[], float result[], size_t buflen)
{
    BlockWidth_8;
    PointerFloat_1;
    __m256 load;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load = _mm256_loadu_ps(src);
        load = _mm256_rcp_ps(load);
        _mm256_storeu_ps(dst, load);
        PointerMove_1;
    }

    for(size_t i = 0; i < nRemain && src[i] != 0; i++)
    {
        dst[i] = 1 / src[i];
    }
}

void SIMD::_rsqrt(float buf[], float result[], size_t buflen)
{
    BlockWidth_8;
    PointerFloat_1;
    __m256 load;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load = _mm256_loadu_ps(src);
        load = _mm256_rsqrt_ps(load);
        _mm256_storeu_ps(dst, load);
        PointerMove_1;
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        if(src[i] > 0)
        {
            dst[i] = 1 / sqrt(src[i]);
        }
        else
        {
            dst[i] = 0;
        }
    }
}

void SIMD::_sqrt(float buf[], float result[], size_t buflen)
{
    BlockWidth_8;
    PointerFloat_1;
    __m256 load;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load = _mm256_loadu_ps(src);
        load = _mm256_sqrt_ps(load);
        _mm256_storeu_ps(dst, load);
        PointerMove_1;
    }

    for(size_t i = 0; i < nRemain && src[i] != 0; i++)
    {
        dst[i] = sqrt(src[i]);
    }
}

void SIMD::_sqrt(double buf[], double result[], size_t buflen)
{
    BlockWidth_4;
    PointerDouble_1;
    __m256d load;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load = _mm256_loadu_pd(src);
        load = _mm256_sqrt_pd(load);
        _mm256_storeu_pd(dst, load);
        PointerMove_1;
    }

    for(size_t i = 0; i < nRemain && src[i] > 0; i++)
    {
        dst[i] = sqrt(src[i]);
    }
}

void SIMD::_sub(float buf1[], float buf2[], float result[], size_t buflen)
{
    BlockWidth_8;
    PointerFloat_2;
    __m256 load1, load2;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load1 = _mm256_loadu_ps(src1);
        load2 = _mm256_loadu_ps(src2);
        load1 = _mm256_sub_ps(load1, load2);
        _mm256_storeu_ps(dst, load1);
        PointerMove_2;
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        dst[i] = src1[i] - src2[i];
    }
}

void SIMD::_sub(double buf1[], double buf2[], double result[], size_t buflen)
{
    BlockWidth_4;
    PointerDouble_2;
    __m256d load1, load2;

    for(size_t i = 0; i < nBlocks; i++)
    {
        load1 = _mm256_loadu_pd(src1);
        load2 = _mm256_loadu_pd(src2);
        load1 = _mm256_sub_pd(load1, load2);
        _mm256_storeu_pd(dst, load1);
        PointerMove_2;
    }

    for(size_t i = 0; i < nRemain; i++)
    {
        dst[i] = src1[i] - src2[i];
    }
}

void SIMD::_ln(float buf[], float result[], double ln_lookup[], size_t buflen)
{
    for(size_t i = 0; i < buflen; i++)
    {
        result[i] = (float)fastlog((double)buf[i], ln_lookup);
    }
}

void SIMD::_ln(double buf[], double result[], double ln_lookup[], size_t buflen)
{
    for(size_t i = 0; i < buflen; i++)
    {
        result[i] = fastlog(buf[i], ln_lookup);
    }
}
