
#include "initial_conditions.h"


void inicon(double **s, int lx, int ly)
{ int i,j;
	for(j=0;j<ly;j++)
	{	
		for(i=0;i<lx/3;i++)
			s[i][j]=1;
		for(i=lx/3;i<lx;i++)
			s[i][j]=0;
	}	
}

