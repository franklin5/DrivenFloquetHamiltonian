//
//  lgwt.h
//  tBdG
//
//  Created by Dong Lin on 1/29/14.
//  Copyright (c) 2014 Dong Lin. All rights reserved.
//

#ifndef LGWT_H_
#define LGWT_H_
#include "stdcpp.h"
double maxfunabs(const int N, const double *y, const double *y0);
void gauss_lgwt(int N, const double a,
		const double b, double x[], double w[]);
#endif
