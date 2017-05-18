/*
 *  coaltrees.h
 *  HKAdirect
 *
 *  Created by Sebastian E. Ramos Onsins on 12/03/2013.
 *  Copyright 2013 CRAG. All rights reserved.
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

float li_mean(long int *S, long int len);
float li_var(long int *S, long int len);
int sample(long int nrep, long int len, int *data, int replace, float *prob, int *result);
int unique(int *data, int lendata, int *result);
double ran1(void);