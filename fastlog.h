#ifndef FASTLOG_H
#define FASTLOG_H

#include <math.h>
#include <stdint.h>

double *fastlog_init(int prec);
void fastlog_free(double *fastlog_lookup);
double fastlog(double x, double *fastlog_lookup);

#endif


