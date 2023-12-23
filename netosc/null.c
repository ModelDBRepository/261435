#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <limits.h>
#include "excon.h"
#include "exvc.c"
#define V_MIN -100
#define OFFSET 50
#define V_RANGE 70  /* actual range divided by increment */
#define V_INC 0.01
#define H_INC 0.01
#define H 100

 
int sgn(x)
double x;
{if (x>=0) return 1;
 if(x<0) return -1;
}

/* global variable declarations */

double current[C],vset1;
 
void main()
{
int stop,start,flag,istop,istart,iflag,kstart,kstop;
double n,state[N],h1v1[2*OFFSET+1],v2v1[2*OFFSET+1];
double time_=0,vnull,hnull;
int i,j,k,counter,index;
FILE  *vp,*hp;
double dv,dh;
double olddv,olddh;
double m1,m2,b1,b2;
vp= fopen("vn.data","w");
hp= fopen("hn.data","w");

for(k=0;k<10000;k++)
{
state[V_1] = V_MIN +k*V_INC ;
fprintf(hp,"%f %f\n",state[V_1],hbar(state[V_1]));
state[H_1] = 0.0;
state[H_2] = 0.0;
state[V_2] = -10.0;
current_(state,time_);
olddv= -1000*current[I_SOMA_1];
for(j=1;j<=H;j++)
{
state[H_1]= H_INC*j;
current_(state,time_);
dv= -1000*current[I_SOMA_1];
if(sgn(dv)!=sgn(olddv))
     {
 vnull = state[H_1] - H_INC*fabs(dv/(dv - olddv));
 fprintf(vp,"%f %f\n",state[V_1],vnull);
      }
   olddv = dv;
}
}
}
