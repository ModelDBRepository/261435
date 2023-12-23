#include <stdio.h>
#include <math.h>
#include "newcon.h"
#include "hvc.c"
double current[C];
void main()
{
double v;
int i;
v=-50;
for(i=0;i<50;i++)
 {
 printf("%f %f\n", v + 0.1*i, sbar(v + 0.1*i));
 }
}
