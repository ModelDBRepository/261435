#include <math.h>
#include <stdio.h>
#include "netcon.h"


double boltz(v,half,slope)
double v,half,slope;
{
double arg;
arg = -(v-half)/slope;
/* if(fabs(arg<0.01)) return(1.0/(2.0 + expm1(arg)));
else*/
return(1.0/(1.0 + exp(arg)));}

double mbar(v)
double v;
{return( boltz(v,M_HALF,M_SLOPE));}

double hbar(v)
double v;
{return( boltz(v,H_HALF,H_SLOPE));}
 
double sbar(v)
double v;
{return( boltz(v,S_HALF,S_SLOPE));}
 
/* CURRENTS */

void current_(sv,time_)
double sv[N];   /* state vector */
double time_;
{
extern double current[C],vset1; 
current[I_L_1]=  G_L*(sv[V_1] - E_L);
current[I_L_2]=  G_L*(sv[V_2] - E_L);
current[I_CA_T_1] = G_CA_T*pow(mbar(sv[V_1]),3.0)*sv[H_1]*(sv[V_1] - E_CA);
current[I_CA_T_2] = G_CA_T*pow(mbar(sv[V_2]),3.0)*sv[H_2]*(sv[V_2] - E_CA);

current[I_SYN_1] = G_SYN*sbar(sv[V_2])*(sv[V_1] - E_SYN);
current[I_SYN_2] = G_SYN*sbar(sv[V_1])*(sv[V_2] - E_SYN); 
current[I_SOMA_1] = current[I_L_1] + current[I_CA_T_1] + current[I_SYN_1] ;
current[I_SOMA_2] = current[I_L_2] + current[I_CA_T_2] + current[I_SYN_2] ;
}


int deriv_(np,xp,Y,F)
double *F,*Y;
int *np;
double *xp;
{
extern double current[C]; 
int el;
double time_;
time_ = *xp;

el = *np;
current_(Y,time_);

F[V_1] =    ( -current[I_SOMA_1])/(CM); /*to get units in mV vs V*/
F[V_2] =    ( -current[I_SOMA_2])/(CM); /*to get units in mV vs V*/
F[H_1] = (-Y[H_1] + boltz(Y[V_1],H_HALF,H_SLOPE))/(TAU_H1*hbar(Y[V_1]));
F[H_2] = (-Y[H_2] + boltz(Y[V_2],H_HALF,H_SLOPE))/(TAU_H2*hbar(Y[V_2]));
return 0;
   }


void scan_(Y) 
double Y[N];
{FILE *fopen(),*sp;
int i;
sp = fopen("state.data","r");
for(i=0;i<N;i++) fscanf(sp,"%lf\n",&Y[i]);
fclose(sp);}

void dump_(Y) 
double Y[N];
{FILE *fopen(),*sp;
int i;
sp = fopen("end.data","w");
for(i=0;i<N;i++) fprintf(sp,"%.16f\n",Y[i]);
fclose(sp);}

int mas(n,amas,l)
        int *n;
        double *amas;
        int *l;
{return 0;}

int dummy(n,t,y,ydot)
        int *n;
        double *t;
        double *y;
        double *ydot;
{return 0;}
