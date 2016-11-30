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
#include "splitting_routines.h"



//--------SHIFTS SEQUENCE TO LEFT FOR 2l/3-------------//

void shift_arr_l(double **s, int lx, int ly, double dx, int *pozicija, int *pozicija_cut, int *pozicija_end, int *brojac, int *brojac_cut, int *brojac_end, double *path, double *path_cut, double *path_end)
{  int i,j; 

	for(j=0;j<ly;j++)
	{	
		for(i=0;i<2*lx/3;i++)
			s[i][j]=s[i+lx/3][j];

		for(i=0;i<lx/3;i++)
			s[i+2*lx/3][j]=0;
	}
    
    
    *path = *path + (lx/3-(*brojac))*dx;
    *brojac = 0;
    *pozicija = lx/3;
				
    *path_cut = *path_cut + (lx/3-(*brojac_cut))*dx;
    *brojac_cut = 0;
    *pozicija_cut = lx/3;
    
    *path_end = *path_end + (lx/3-(*brojac_end))*dx;
    *brojac_end = 0;
    *pozicija_end = lx/3;

    
}


//---------Fourier transform----FRONT------------------//

void fftw_line(double ** sq, double *poz_line, int ly, int j1)
{ 	fftw_plan p;
    fftw_complex *izlaz;

    
izlaz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ly);  // def complex sequence
p= fftw_plan_dft_r2c_1d(ly,poz_line,izlaz,FFTW_ESTIMATE);
fftw_execute(p);


sq[0][j1]=sq[0][j1] + (pow(izlaz[0][0],2) + pow(izlaz[0][1],2));
for(int iy=1;iy<ly/2;iy++)
    sq[iy][j1]=sq[iy][j1] + (pow(izlaz[iy][0],2)  + pow(izlaz[iy][1],2))/(ly);              // Structure factor

fftw_destroy_plan(p);
free(izlaz);
}

void integrate_noise(double **s, int pozicija, int ly, double dt, double sigma, double lambda, long int *inic)
{   double prev_value, new_value;
    int inic_noise;

    for(int iy=0;iy<ly;iy++)
    {
    
        inic_noise=pozicija-50;
        do inic_noise++;
        while (s[inic_noise][iy]>0.5);    // inic_noise, implements noise for values of rho<0.5
    
        for(int ix=inic_noise;ix<pozicija+250;ix++)
        {
            if(s[ix][iy]*lambda>1000)
                s[ix][iy]=s[ix][iy]+sigma*sqrt(s[ix][iy]*dt)*gauss(0,1,inic);
        
            else
                if(s[ix][iy]*lambda>0.0001)
                {   prev_value=s[ix][iy];
                    new_value=gamdev(poidev(lambda*s[ix][iy],inic),inic);
                    s[ix][iy]=new_value/lambda;
            
                    if(!isfinite(s[ix][iy]) || s[ix][iy]!=s[ix][iy])
                    {   s[ix][iy]=prev_value;
                        cout<< s[ix][iy] << "   2.deo p.!!  " << prev_value << "   " << ix << "  " << iy << "  "<< new_value << endl;
                    }
                }
        
            else s[ix][iy]=0;
        }

    }
}

void integrate_diffusion(double **s, double *firstx, double *lastx1, double *lastx2,  int pozicija, int lx, int ly, double dx, double dt, double D)
{   double last1, last2;
    int ix,iy;

    last1=s[0][0];
    firstx[0]=1;
    firstx[lx-1]=0;
    last1=s[pozicija-101][0];
    
    for(ix=pozicija-100;ix<pozicija+252;ix++)                      //lx-1   pozicija+50
    {
        last2=s[ix][0];							// FIRST WAVE, ON THE LEFT
        firstx[ix]=s[ix][0];
        lastx1[ix]=s[ix][0];
        s[ix][0]=s[ix][0] + D*dt*(s[ix+1][0] + last1 + s[ix][1] + s[ix][ly-1] - 4*s[ix][0])/(dx*dx) + s[ix][0]*dt - s[ix][0]*s[ix][0]*dt;       // Deterministic term b=1
        last1=last2;
    }
    
    for(iy=1;iy<ly-1;iy++)
    {
        last1=1; last1=s[pozicija-101][iy];
        for(ix=pozicija-100;ix<pozicija+252;ix++)                        // WAVES INSIDE THE REGION
        {
            last2=s[ix][iy];
            lastx2[ix]= s[ix][iy];
            s[ix][iy] = s[ix][iy] + D*dt*(s[ix+1][iy] + last1 + lastx1[ix] + s[ix][iy+1]  - 4*s[ix][iy])/(dx*dx) + s[ix][iy]*dt - s[ix][iy]*s[ix][iy]*dt;      // Deterministic term b=1
            last1=last2;
            lastx1[ix]=lastx2[ix];
        }
    }
    last1=1; last1=s[pozicija-101][iy];
    for(ix=pozicija-100;ix<pozicija+252;ix++)                      // LAST WAVE
    {
        last2=s[ix][ly-1];
        s[ix][ly-1]=s[ix][ly-1] + D*dt*(s[ix+1][ly-1] + last1 + firstx[ix] + lastx1[ix] - 4*s[ix][ly-1])/(dx*dx) + s[ix][ly-1]*dt - s[ix][ly-1]*s[ix][ly-1]*dt;      // Deterministic term b=1
        last1=last2;
    }
}


void check_wave_positions(double **s, int pozicija, int pozicija_cut, int pozicija_end, double N)
{
    if(s[pozicija][0]>0.5)
        {do pozicija=pozicija+1; while (s[pozicija][0]>0.5);}               //locating the front
    
    if(s[pozicija_cut][0]>1/N)
        {do pozicija_cut=pozicija_cut+1; while (s[pozicija_cut][0]>1/N);}	//locating the edge
    
    if(s[pozicija_end][0]>0)
        {do pozicija_end=pozicija_end+1; while (s[pozicija_end][0]>0);}     //locating the edge
}


void calculate_positions(double **s, int lx, int ly, int j1, double **y_fce_t, double *poz, double *omega, int *brojac, double *path_fce, double dx)
{   int iy, ix = *brojac-1;
    double path,h=0;
    do  ix++;
    while(s[ix+lx/3][0]>0.5);
    
    path = path + (ix-(*brojac))*dx;
    poz[0] = (path*(s[ix+lx/3-1][0]-s[ix+lx/3][0])+s[ix+lx/3][0]-0.5)/(s[ix+lx/3-1][0]-s[ix+lx/3][0]);
    y_fce_t[j1][0] = y_fce_t[j1][0] + poz[0];
    
    *brojac=ix;
    
    // Positions of other waves:
    
    for(iy=1;iy<ly;iy++)
    {
        ix = *brojac;
        if(s[ix+lx/3][iy]>0.5)
        	do  ix++;
            while(s[ix+lx/3][iy]>0.5);
        
        else
            do  ix--;
            while(s[ix+lx/3][iy]<0.5); ix++;
        
        
        poz[iy] = ((path+(ix-(*brojac))*dx)*(s[ix+lx/3-1][iy]-s[ix+lx/3][iy])+s[ix+lx/3][iy]-0.5)/(s[ix+lx/3-1][iy]-s[ix+lx/3][iy]);
        y_fce_t[j1][iy] = y_fce_t[j1][iy] + poz[iy];
    }
    
    
    for(iy=0;iy<ly;iy++)
        h = h + poz[iy];
    h=h/ly;
    
    for(iy=0;iy<ly;iy++)
        omega[j1]=omega[j1]+(poz[iy]-h)*(poz[iy]-h)/ly;
    
    *path_fce=path;
}



