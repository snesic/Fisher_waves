/*
 *
 *  Created by Svetozar Nesic on 26.2.10..
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *  Numerical recepies codes
 */
#ifndef __RANDOM_NO_GENERATORS_H_INCLUDED__
#define __RANDOM_NO_GENERATORS_H_INCLUDED__

#include <iostream>
#include <stdlib.h>
#include <math.h>

using namespace std;


//-----Random generators - Numerical Recepies----Poisson and Gamma------------//



double ran1(long *);    // ovo je ran1 al da ne bi menjao u svim procedurama..

double gamdev(int, long *);

double gammln(double);					//Returns the value ln[Γ(xx)] for xx > 0.

double poidev(double, long *);

double gauss(double, double s, long *);			/* normal random variate generator */
				     				   /* mean m, standard deviation s */
double ran3(long *);


#endif