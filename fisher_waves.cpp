/*
 *  fisher_waves.cpp
 *
 Simulation of stochastic FISHER WAVES in 2D.  ROUGHNESS FUNCTION, PSD AND KPZ BEHAVIOR, CORRELATION BETWEEN the EDGE AND FRONT..
 
 The algorithm employs the splitting step scheme to solve the stochastic partial differential equation.
 *
 Initial condition:    Step function!!!        other options:       f(x-ct + epsilon*sin(q*y)),  where q=2*PI*i/ly...
 *
 Traveling direction: x
 *
 All variables are normalized!!!
 *
 Variables:
 *
 1. roughness function w(t)  ->    roughness.txt
 2. Sq(t)                    ->    matsq.txt
 3. front(t)                 ->    mat_front.txt -> equipotential line where rho(x,y) = 1/2
 4. edge (t)                 ->    mat_edge.txt  -> equipotential line where rho(x,y) = 1/N
 5. end (t)                 ->    mat_edge.txt   -> equipotential line where rho(x,y) = 0
 *
 *  Created by Svetozar Nesic on 26.2.10..
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <time.h>
#include "read_write_msg.h"
// #include "initial_conditions.h"
#include "simulation.h"

using namespace std;


int main(int argc, char * const argv[])
{   double **s;
    int lx, ly, i, iy, brojac,no_int,i_no;
	long int t;
	double dx, dt, D, sigma;
    time_t calc_time;
    long int *inic = new long int;
    
    double *x_t, **y_t, **y_cut_t, **y_end_t;
    double *w, *w_cut, *w_end;
    double **Sq, **Sq_edge, **Sq_end;
    int no_points;
    

					// D,sigma:  	 	     - system constants. calc_time: calculation time
					// dx, dy=dx, dt 	     - discretization of space and time. t,lx, ly: system size. *s: system function.
					// w, w_cut 	 - averaged rougness of the front and edge. 
					// w_fourier, Sq 	     - a. roughness in fourier space,  structure factor averaged
					// y_t, y_cut_t, x_t     - coordinates of the front line, edge line. 
					// iy , i, i_no          - counters lateral size, usually time counter, realizations

calc_time = time(NULL);

no_int=1;

//system parameters:

    
    if (argc==2) // When the program runs in parallel (by script), only one realization is performed (no_int=1) 
    {   ly =  atoi(argv[1]);
        no_int = 1;
    }
    else ly = 2048;
    if (ly>2048) {ly=2048; cout << "ly too big, 2048 is used"<< endl;}

	sigma = 1e-3;
	lx  = 900;
	t   = 40;
	dt  = 0.05;
	dx  = 1;
	D   = 1;
    no_points=2; //t/2000;
    *inic=45;
    //*inic=time(NULL)*500;

    

	w  		= new double [no_points];
	w_cut  	= new double [no_points];
	w_end  	= new double [no_points];


	s  = new double *[lx];
	for(i=0;i<lx;i++)
		s[i] = new double [ly]; 


    initial_message(lx, ly, t, dt, dx, no_int, sigma);

    
// our parameters:


    x_t	= new double  [no_points];
        for(i=0;i<no_points;i++)
            x_t[i] = t/no_points*i*dt;
    
    y_t	= new double *[no_points];
            for(i=0;i<no_points;i++)
				y_t[i] = new double [ly];

	y_cut_t	= new double *[no_points];
			for(i=0;i<no_points;i++)
				y_cut_t[i] = new double [ly];

	y_end_t	= new double *[no_points];
			for(i=0;i<no_points;i++)
				y_end_t[i] = new double [ly];

	Sq 	= new double *[ly/2];
			for(i=0;i<ly/2;i++)
				Sq[i] = new double [no_points];

	Sq_edge = new double *[ly/2];
			for(i=0;i<ly/2;i++)
				Sq_edge[i] = new double [no_points];


	Sq_end = new double *[ly/2];
			for(i=0;i<ly/2;i++)
				Sq_end[i] = new double [no_points];

	
    for(i=0;i<no_points;i++)
    {
        w[i]=0;
        w_cut[i]=0;
        w_end[i]=0;
        for(iy=0;iy<ly;iy++)
        {
            y_t[i][iy]=0;
            y_cut_t[i][iy]=0;
            y_end_t[i][iy]=0;
        }
    }
    
	
    for(i_no=0;i_no<no_int;i_no++)  //----no_int simulations-------
    {
        stochastic(s, lx, ly, t, no_points, dx, dt, D, sigma, x_t, y_t, y_cut_t, y_end_t, w, w_cut, w_end, Sq, Sq_edge, Sq_end, inic );
        cout << "Calculation time: " << time(NULL)-calc_time << "     " << i_no << "\n" <<endl;
        calc_time=time(NULL);
    }
    
    
    for(i=0;i<no_points;i++)
    {
        w[i]     = sqrt(w[i]/no_int);
        w_cut[i] = sqrt(w_cut[i]/no_int);
        w_end[i] = sqrt(w_end[i]/no_int);
    }
	
    cout<< y_t[1][120] << " " << y_cut_t[1][120] << endl;
    
    write1d_data(x_t, w    , no_points, "roughness_front.txt");
    write1d_data(x_t, w_cut, no_points, "roughness_edge.txt" );
//    write1d_data(x_t, w_end, no_points, "roughness_end.txt"  );

    write2d_data_sq(Sq, ly, no_points, no_int, "matsq.txt");
    write2d_data_sq(Sq_edge, ly, no_points, no_int, "matsq_edge.txt");
//    write2d_data_sq(Sq_end, ly, no_points, no_int, "matsq_end.txt");

    write2d_data(y_t, ly, no_points, no_int, "mat_front.txt");
//    write2d_data(y_cut_t, ly, no_points, no_int, "mat_edge.txt");
	
return 0;

}
