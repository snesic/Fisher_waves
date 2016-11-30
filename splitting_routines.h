/*
 *  splitting.cpp
 *  
 
 kompajlirati: g++ split_klasa.cpp -O4 -lfftw3 -lm -o proba
 
 
 FISHER WAVES  -   2D      ROUGHNESS FUNCTION, PSD AND KPZ BEHAVIOR, CORRELATION BETWEEN EDGE AND THE MIDDLE POINT OF THE FRONT, MAYBE, INSTABILITIES APPEAR
 
 
 Za pocetni uslov    step funkcija!!!        druga opcija:       f(x-ct + epsilon*sin(q*y)),  gde je q=2*PI*i/ly
 
 Sve velicine su normirane!!!
 
 racuna: 
 
 1. roughness function     ->	 roughness1.txt
 2. Sq in time             ->    matsq.txt
 3. rougness front mat     ->    mat_front.txt
 4. rougness edge mat      ->    mat_edge.txt
 
 
 5. roughness over max sq  ->    omega_in_sqmax.txt
 jos ne funkcionise 4. max sq over time       ->    sqmax_in_t.txt
 
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
//#include <complex.h>
//#include "/opt/local/include/fftw3.h"
#include <fftw3.h>
#include "random_no_generators.h"

#define PI 3.141592653589793

using namespace std;

void shift_arr_l(double **, int, int, double, int *, int *, int *, int *, int *, int *, double *, double *, double *); //SHIFTS SEQUENCE TO LEFT FOR 2l/3//


void fftw_line(double **, double *, int, int); //---------Fourier transform----FRONT------------------//

void integrate_noise(double **, int, int, double, double, double, long int *);

void integrate_diffusion(double **, double *, double *, double *,  int, int, int, double, double, double);

void check_wave_positions(double **, int, int, int, double);

void calculate_positions(double **, int, int, int, double **, double *, double *, int *, double *, double);


