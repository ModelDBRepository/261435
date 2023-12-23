/* solout.c */
#include <math.h>
#include <stdio.h>
#include "netcon.h"

#define T_RES 100
#define V_LOW_RES 2.0
#define V_HIGH_RES 0.05

int out_(nr,xold,x,state,n,irtrn)
        int *nr;
        double *xold;
        double *x;
        double *state;
        int *n;
        int *irtrn;
{
static double t_old,yold,yolder;
static double dyold,dyolder;
extern double current[C];
int i;


if((*nr==1) || (fabs(yold -state[0])>V_LOW_RES) ||
(fabs(dyold - state[1])>V_LOW_RES) ||
(yold>state[0] && yold>yolder)  ||
(yold<state[0] && yold<yolder)  ||
(dyold>state[1] && dyold>dyolder)  ||
(dyold<state[1] && dyold<dyolder)  ||
 (fabs(*x - t_old) > T_RES))
{
printf("%f ",*x);
for(i=0;i<N;i++) printf("%f ", state[i]);
printf("\n");
fflush(stdout);
yolder=yold;
yold=state[0];
dyolder=dyold;
dyold=state[N];
t_old=*x;}
}
