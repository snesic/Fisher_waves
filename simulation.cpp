#include "simulation.h"


void stochastic(double **s, int lx, int ly, long int t, int no_points, double dx, double dt, double D, double sigma, double *x_t, double **y_t, double **y_cut_t, double **y_end_t, double *w, double *w_cut, double *w_end, double **Sq, double **Sq_edge, double **Sq_end, long int *inic )
{
	int ix, iy, brojac, brojac_cut, brojac_end, pozicija, pozicija_cut,pozicija_end;
	long int j , j1 , vcounter; 
	double path, path_cut, path_end, *poz_front, *poz_edge,*poz_end, *lastx1, *lastx2, *firstx,h,N;

	double lambda;
	lambda=2/(sigma*sigma*dt);				// samo za inicijalicaziju
	N=pow(1/sigma,2);	
	
			// samo za inicijalicaziju
	
	//---------INITIAL CONDITIONS-------------//
	
	inicon(s,lx,ly);
	pozicija=lx/3;
	pozicija_cut=lx/3;
	pozicija_end=lx/3;
	
	// speed parameters:	

	path=0;
	brojac=0;
	path_cut=0;
	brojac_cut=0;
	path_end=0;
	brojac_end=0;

	poz_front 	= new double [ly];       // positions of the front of a current simulation
	poz_edge 	= new double [ly];
	poz_end 	= new double [ly];

	lastx1 = new double[lx];
	lastx2 = new double[lx];
	firstx = new double[lx];
    
	for(ix=0;ix<lx;ix++) 
	{
        lastx1[ix]=s[ix][1];
        firstx[ix]=s[ix][1];
    }
	
	//-------------------SIMULATING EQUATION---------------------------//

	for(j1=0;j1<no_points;j1++)
	{
		for(j=j1*t/no_points;j<(j1+1)*t/no_points;j++)
		{  
			
            integrate_noise(s, pozicija, ly, dt, sigma, lambda, inic);  // Splitting step method: integrating stochastic part first:
            integrate_diffusion(s, firstx, lastx1, lastx2, pozicija, lx, ly, dx, dt, D); // Integrating deterministic part:

            check_wave_positions(s, pozicija, pozicija_cut, pozicija_end, N);
			
			if(s[2*lx/3-1][0]>s[0][0]/2 && s[2*lx/3][0]<s[0][0]/2) // SHIFTING IF FRONT HAS REACHED l/3
				shift_arr_l(s, lx, ly, dx, &pozicija, &pozicija_cut, &pozicija_end, &brojac, &brojac_cut, &brojac_end, &path, &path_cut, &path_end);
		}
		
		
		// Coordinate of the (FRONT, EDGE, END) in time: First calculates the position of the first wave(y=0) and then comapres it with the positions of the other waves.
		
        calculate_positions(s, lx, ly, j1, y_t    , poz_front, w    , &brojac    , &path    , dx);
        calculate_positions(s, lx, ly, j1, y_cut_t, poz_edge , w_cut, &brojac_cut, &path_cut, dx);
        calculate_positions(s, lx, ly, j1, y_end_t, poz_end  , w_end, &brojac_end, &path_end, dx);
        
        fftw_line(Sq, poz_front, ly, j1);     // Fourier transform----FRONT------
        fftw_line(Sq_edge, poz_edge, ly, j1); // Edge
        fftw_line(Sq_end, poz_end, ly, j1);   // End
				
		
	}	// end of the j1 counter, time counter
    
}

