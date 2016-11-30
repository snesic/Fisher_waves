/*
 *
 FISHER WAVES  -   2D      ROUGHNESS FUNCTION, PSD AND KPZ BEHAVIOR, CORRELATION BETWEEN EDGE AND THE MIDDLE POINT OF THE FRONT, MAYBE, INSTABILITIES APPEAR

 Za pocetni uslov    step funkcija!!!        druga opcija:       f(x-ct + epsilon*sin(q*y)),  gde je q=2*PI*i/ly
 
 *
 *  Created by Svetozar Nesic on 26.2.10..
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include "splitting_routines.h"
#include "initial_conditions.h"

#define PI 3.141592653589793

using namespace std;



void stochastic(double **, int, int, long int, int, double, double, double, double, double *, double **, double **, double **, double *, double *, double *, double **, double **, double **, long int * );
