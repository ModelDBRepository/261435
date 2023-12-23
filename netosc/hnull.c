#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <limits.h>
#include "netcon.h"
#include "vc.c"
#define V_MIN -50
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
double hnull;
double time_=0;
int i,j,k,counter,index;
FILE  *vp,*hp;
double dv;
double olddv;
double m1,m2,b1,b2;
vp= fopen("vn.data","w");
hp= fopen("hn.data","w");

for(k=0;k<=H;k++)
{
state[H_1]= k*H_INC;
for(i=0;i<=2*OFFSET;i++) v2h1[i] = 1000;
for(i=0;i<=2*OFFSET;i++) v1h1[i] = 1000;
for(i=0;i<=2*OFFSET;i++) h2h1[i] = 1000;
for(i=0;i<=2*OFFSET;i++)
{
state[V_1] = V_MIN - OFFSET +i*V_INC;
state[V_2] = state[V_1] - OFFSET;
state[H_2] = hbar(state[V_2]);
current_(state,time_);
olddv= -current[I_SOMA_1]/(CM);
for(j=1;j<=2*OFFSET;j++)
{
state[V_2] = state[V_1] - OFFSET +j*V_INC;
state[H_2] = hbar(state[V_2]);
current_(state,time_);
dv= -current[I_SOMA_1]/(CM);
if(sgn(dv)!=sgn(olddv))
     {
 v2h1[i] = state[V_2] - V_INC*fabs(dv/(dv - olddv));
 v1h1[i] = state[V_1]; 
      }
   olddv = dv;
}
}
state[V_1]=v1h1[0];
state[V_2]=v2h1[0];
state[H_2]=hbar(state[V_2]);
current_(state,time_);
olddv= -current[I_SOMA_2]/(CM);
for(i=0;i<=2*OFFSET;i++)
  {
state[V_1]=v1h1[i];
state[V_2]=v2h1[i];
state[H_2]=hbar(state[V_2]);
current_(state,time_);
dv= -current[I_SOMA_2]/(CM);
if(sgn(dv)!=sgn(olddv))
     {
printf("%f %f %f %f\n", state[H_1],v1h1[i],v2h1[i],state[H_2]);
      }
   olddv = dv;
  }
}
}
