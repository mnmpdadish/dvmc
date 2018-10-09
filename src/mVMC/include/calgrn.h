#ifndef _CALGRN
#define _CALGRN
#include <complex.h>

void CalculateGreenFunc(const double w, const double complex ip, int *eleIdx, int *eleCfg,
                         int *eleNum, int *eleProjCnt);

void CalculateGreenFuncMoments(const double w, const double complex ip, int *eleIdx, int *eleCfg,
                         int *eleNum, int *eleProjCnt);

void CalculateGreenFuncBF(const double w, const double ip, int *eleIdx, int *eleCfg,
                          int *eleNum, int *eleProjCnt, const int *eleProjBFCnt);
#endif
