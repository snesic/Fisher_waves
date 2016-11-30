/*
 *
 *  Created by Svetozar Nesic on 26.2.10..
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include "random_no_generators.h"

#define PI 3.141592653589793
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

//-----Random generators - Numerical Recepies----Poisson and Gamma------------//



double ran1(long *idum)    // ovo je ran1 al da ne bi menjao u svim procedurama..
{
int j;
long k;
static long iy=0;
static long iv[NTAB];
float temp;
if (*idum <= 0 || !iy) {
			if (-(*idum) < 1) *idum=1; 
			else *idum = -(*idum);
			for (j=NTAB+7;j>=0;j--) { 
							k=(*idum)/IQ;
							*idum=IA*(*idum-k*IQ)-IR*k;
							if (*idum < 0) *idum += IM;
							if (j < NTAB) iv[j] = *idum;
						}
						iy=iv[0];
			}
k=(*idum)/IQ; 
*idum=IA*(*idum-k*IQ)-IR*k; 
if (*idum < 0) *idum += IM; 
j=iy/NDIV; 
iy=iv[j];
iv[j] = *idum; 
if ((temp=AM*iy) > RNMX) return RNMX; 
else return temp;
}


double ran3(long *idum)
{
static int inext,inextp;
static long ma[56]; 
static int iff=0;
long mj,mk;
int i,ii,k;


	if (*idum < 0 || iff == 0) { 

			iff=1;
			mj=labs(MSEED-labs(*idum)); 
			mj %= MBIG; 
			ma[55]=mj;
			mk=1;
			for (i=1;i<=54;i++) { 
					ii=(21*i) % 55; 
					ma[ii]=mk; 
					mk=mj-mk;
					if (mk < MZ) mk += MBIG;
				        mj=ma[ii];
			}

			for (k=1;k<=4;k++) 
				for(i=1;i<=55;i++) {
							ma[i] -= ma[1+(i+30) % 55];
							if (ma[i] < MZ) ma[i] += MBIG;
				    		}
	inext=0; 
	inextp=31; 
	*idum=1;
	}

	if (++inext == 56) inext=1; 
	if (++inextp == 56) inextp=1; 
	mj=ma[inext]-ma[inextp]; 
	if (mj < MZ) mj += MBIG; 
	ma[inext]=mj;
 
return mj*FAC;
}


double gamdev(int ia, long *idum)						 						
{

double ran3(long *idum);
//void nrerror(char error_text[]);
int j;
double am,e,s,v1,v2,x,y;
if (ia < 1) return 0;
if (ia < 6) { 
x=1.0; 
for (j=1;j<=ia;j++) x *= ran3(idum);
x = -log(x);
} else { 
	do {
		do {
			do { 	v1=ran3(idum);
				v2=2.0*ran3(idum)-1.0;
			} while (v1*v1+v2*v2 > 1.0);
			y=v2/v1;
			am=ia-1;
			s=sqrt(2.0*am+1.0);
			x=s*y+am;
		} while (x <= 0.0); 
		e=(1.0+y*y)*exp(am*log(x/am)-s*y); 
	} while (ran3(idum) > e); 
       } 
return x;
}


double gammln(double xx)					//Returns the value ln[Γ(xx)] for xx > 0.
{
//Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
//accuracy is good enough.
double x,y,tmp,ser;
static double cof[6]={76.18009172947146,-86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;

return -tmp+log(2.5066282746310005*ser/x);
}



double poidev(double xm, long *idum)
{
double gammln(double xx);
double ran3(long *idum);
static double sq,alxm,g,oldm=(-1.0); 
double em,t,y; 

	if (xm < 12.0) { 
		if (xm != oldm) {
				oldm=xm;
				g=exp(-xm); 
				}
		em = -1;
		t=1.0;
		do { 
			++em;
			t *= ran3(idum);
		} while (t > g);	
	} else { 
		if (xm != oldm) { 
				sq=sqrt(2.0*xm);
				alxm=log(xm);
				g=xm*alxm-gammln(xm+1.0);
				}
		do {
			do { y=tan(PI*ran3(idum));
			     em=sq*y+xm; 
			   } while (em < 0.0); 
		em=floor(em); 
		t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);

		} while (ran3(idum) > t);
	}

return em;
}


double gauss(double m, double s, long *idum)			/* normal random variate generator */
{				     				   /* mean m, standard deviation s */
	
double ran3(long *idum);
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * ran3(idum) - 1.0;
			x2 = 2.0 * ran3(idum) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}
	if(!isfinite(m + y1 * s))  cout<< "Generator!!"<<endl;
	return( m + y1 * s );
}
