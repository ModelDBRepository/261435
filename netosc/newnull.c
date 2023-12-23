#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <limits.h>
#include "netcon.h"
#include "vc.c"
#define V_MIN -70
#define OFFSET 70
#define V_RANGE 70  /* actual range divided by increment */
#define V_INC 1.0
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
double n,state[N],v1h1[2*OFFSET+1],v2h1[2*OFFSET+1],h2h1[2*OFFSET+1];
double v2,vnull,hnull;
double time_=0;
int i,j,k,counter,index;
FILE  *vp,*hp;
double dv;
double olddv;
double m1,m2,b1,b2;
vp= fopen("vn.data","w");
hp= fopen("hn.data","w");

for(k=0;k<100;k++)
{
state[V_1] = V_MIN +k*V_INC ;
state[V_2] = state[V_1] - OFFSET;
state[H_2] = hbar(state[V_2]);
current_(state,time_);
olddv= -current[I_SOMA_2]/(CM);
for(j=1;j<=2*OFFSET;j++)
{
state[V_2] = state[V_1] - OFFSET +j*V_INC;
state[H_2] = hbar(state[V_2]);
current_(state,time_);
dv= -current[I_SOMA_2]/(CM);
if(sgn(dv)!=sgn(olddv))
     {
 v2 = state[V_2] - V_INC*fabs(dv/(dv - olddv));
      }
   olddv = dv;
}
  state[H_1] = 0;
  state[V_2]=v2;
  current_(state,time_);
  olddv= -current[I_SOMA_1]/(CM);
  for(j=0;j<=H;j++)
   {
     state[H_1]= j*H_INC;
     current_(state,time_);
     dv= -current[I_SOMA_1]/(CM);
     if(sgn(olddv/dv)<0) 
       {
        vnull = state[H_1] - H_INC*fabs(dv/(dv - olddv));
        fprintf(vp,"%f %f\n", state[V_1],vnull);
        fprintf(hp,"%f %f\n", state[V_1],hbar(state[V_1]));
       }
     olddv=dv;
   }

}
}
