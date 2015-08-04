#include <stdio.h>
#include <sys/time.h>	//timeval, gettimeofday()
#include <time.h>		//time()
#include <stdlib.h>		//size_t, srand(), rand()
#include <math.h>		//ceil(), floor(), log()
#include "SIMD.h"

using namespace SIMD;

#define BUFLEN 4096
#define _TIME (stop.tv_sec-start.tv_sec+(double)stop.tv_usec/1000000-(double)start.tv_usec/1000000)

void normal_add(double buf1[], double buf2[], double result[], size_t buflen)
{
    for(size_t i = 0; i < buflen; i++)
    {
        result[i] = buf1[i] + buf2[i];
    }
}
void normal_sub(double buf1[], double buf2[], double result[], size_t buflen)
{
    for(size_t i = 0; i < buflen; i++)
    {
        result[i] = buf1[i] - buf2[i];
    }
}
void normal_mul(double buf1[], double buf2[], double result[], size_t buflen)
{
    for(size_t i = 0; i < buflen; i++)
    {
        result[i] = buf1[i] * buf2[i];
    }
}
void normal_div(double buf1[], double buf2[], double result[], size_t buflen)
{
    for(size_t i = 0; i < buflen; i++)
    {
	    if(buf2[i] != 0)
		{
			result[i] = buf1[i] / buf2[i];
		}
    }
}
//void normal_log()

typedef void (*FUNC)(double *, double *, double *, size_t);

void RunTest(double buf1[], double buf2[], double result[], size_t buflen, FUNC testfunc, const char *funcName)
{
	timeval start, stop;	
	gettimeofday(&start, NULL);
    for(int loop = 0; loop < 5000000; loop++)
    {
        testfunc(buf1, buf2, result, buflen);
    }
    gettimeofday(&stop, NULL);
	printf("\n%s time cost: %f s.\n", funcName, _TIME);
	printf("Result check: %f, %f.\n", result[0], result[buflen-1]);
}

int main(int argc, char **argv)
{
    double a[BUFLEN], b[BUFLEN], result[BUFLEN];
    srand(time(NULL));
    for(int i = 0; i < BUFLEN; i++)
    {
        a[i] = (rand()&0x3f) * 1.1;
        b[i] = (rand()&0x3f) * 0.9;

    }
	RunTest(a, b, result, BUFLEN, SIMD::_add, "SIMD_add()");
	RunTest(a, b, result, BUFLEN, normal_add, "Normal_add()");

	RunTest(a, b, result, BUFLEN, SIMD::_sub, "SIMD_sub()");
	RunTest(a, b, result, BUFLEN, normal_sub, "Normal_sub()");

	RunTest(a, b, result, BUFLEN, SIMD::_mul, "SIMD_mul()");
	RunTest(a, b, result, BUFLEN, normal_mul, "Normal_mul()");

	RunTest(a, b, result, BUFLEN, SIMD::_divide, "SIMD_div()");
	RunTest(a, b, result, BUFLEN, normal_div, "Normal_div()");

	return 0;
}
