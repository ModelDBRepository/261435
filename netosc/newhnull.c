/* assumes fast variables are at equilibrium */
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <limits.h>
#include "netcon.h"
#include "vc.c"
#define V_MIN -70
#define OFFSET 70
#define V_RANGE 70  /* actual range divided by increment */
#define V_INC 0.1
#define O_INC 1.0
#define H_INC 0.01
#define H 100
#define K 1000

 
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
double n,state[N],v2v1[K+1],h1v1[K+1],h2v1[K+1];
double v2,vnull,hnull;
double time_=0;
int i,j,k,counter,index;
FILE  *vp,*hp,*v1p,*v2p;
double dv;
double olddv;
double m1,m2,b1,b2;
v1p= fopen("v1.data","w");
v2p= fopen("v2.data","w");
vp= fopen("vn.data","w");
hp= fopen("hn.data","w");

for(k=0;k<K;k++)
{
state[V_1] = V_MIN +k*V_INC ;
state[V_2] = state[V_1] - OFFSET;
state[H_2] = hbar(state[V_2]);
current_(state,time_);
olddv= -current[I_SOMA_2]/(CM);
for(j=1;j<=2*OFFSET;j++)
{
state[V_2] = state[V_1] - OFFSET +j*O_INC;
state[H_2] = hbar(state[V_2]);
current_(state,time_);
dv= -current[I_SOMA_2]/(CM);
if(sgn(dv)!=sgn(olddv))
     {
 v2 = state[V_2] - O_INC*fabs(dv/(dv - olddv));
      }
   olddv = dv;
}
  state[H_1] = 0;
  state[V_2]=v2;
  v2v1[k]=v2;
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
        h1v1[k]=vnull;
        h2v1[k]=hbar(v2);
        fprintf(vp,"%f %f\n", hbar(v2),vnull);
        fprintf(hp,"%f %f\n", vnull, hbar(v2));  
        fprintf(v1p,"%f %f\n",state[V_1],v2);
        fprintf(v2p,"%f %f\n",v2,state[V_1]);
       }
     olddv=dv;
   }

}
/*for(k=0;k<K;k++)
{
state[V_1] = V_MIN +k*V_INC ;
state[V_2] = v2v1[k] ;
state[H_1] = h1v1[k] ;
state[H_2] = h2v1[k] ;
printf("%f %f\n",state[V_1],state[H_1]);
} */
fclose(vp);
fclose(hp);
fclose(v1p);
fclose(v2p);
}
