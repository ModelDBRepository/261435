#include <stdio.h>
#include <math.h>
#include "newcon.h"
#include "hvc.c"
/*#include "solout.c"*/
#include "cblock.h"
/* for sag use newcon.h hvc.c, for pir use netcon.h vc.c */


/* global variable declarations */

/* double F[N];
double Y[N];
long int *np;
double *xp; */
 
double current[C],vset1;
 
void main()
{
double state[N];
double rtol[N],atol[N];

int itol,ijac,mljac,mujac;
int imas,mlmas,mumas;
int iout,lwork,liwork,lrcont;
int iwork[LIWORK];
double work[LWORK];
extern int out_();
int n,idid;

double x,xend,h;
FILE *fp;
int i,j;
n=N;
for(i=0;i<7;i++) {iwork[i]=0;
                  work[i]=0.0;}
for(i=0;i<N;i++) {rtol[i]=1e-13;
                  atol[i]=0.000001;}
itol=1;
ijac=0;
mljac=n;
mujac=0;
imas=0;
mlmas=n;
mumas=0;
lwork=LWORK;
liwork=LIWORK;
lrcont=LRCONT;
iout=1;

h=1e-4;
x=START_TIME;
scan_(state);
current_(state,x);
 
if(x<ENDTIME)
{xend=ENDTIME;

radau5_(&n,deriv_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        out_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);}

dump_(state);
}

