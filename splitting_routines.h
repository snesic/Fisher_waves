/*
 * Routines used by simulation
 *
 shift_arr_l shifts the waves back in order not to leave the box
 *
 fftw_line Fourier transforms the line
 *
 integrate_noise: Samples random numbers as a part of the predictor step
 *
 integrate_diffusion: Corrector step: solves the deterministic equation (diffusion term + reaction term)
 *
 check_wave_positions: for faster calculations (always keep track of rho=1/2) 
 *
 calculate_positions: Calculates the equipotential lines (rho= 1/2, 1/N, 0)
 *
 *
 *  Created by Svetozar Nesic on 26.2.10..
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <time.h>
#include <fftw3.h>
#include "random_no_generators.h"

#define PI 3.141592653589793

using namespace std;

void shift_arr_l(double **, int, int, double, int&, int&, int&, int&, int&, int&, double&, double&, double&); //SHIFTS SEQUENCE TO LEFT FOR 2l/3//


void fftw_line(double **, double *, const int&, const long int&); //---------Fourier transform----FRONT------------------//

void integrate_noise(double **, const int&, const int&, const double&, const double&, const double&, long int *);

void integrate_diffusion(double **, double *, double *, double *,  const int&, const int&, const int&, const double&, const double&, const double&);

void check_wave_positions(double **, int&, int&, int&, const double&);

void calculate_positions(double **, const int&, const int&, const long int&, double, double **, double *, double *, int&, double&, const double&);


