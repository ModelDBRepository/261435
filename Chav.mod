: HH components of chay-cook

TITLE ChayCook
 
UNITS {
        (molar) = (1/liter)
        (S) = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
         F = (faraday) (coulomb)
         R = (mole k) (mV-coulomb/degC)
        (mM) =  (millimolar)
}
 
NEURON {
        SUFFIX chav
        USEION k WRITE ik
        USEION na WRITE ina 
        NONSPECIFIC_CURRENT il
        RANGE ik,ina, il, gna, gk, gl, ena, ek, el, SC1, SC2,hinit
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
    v   (mV)
    dt  (ms)
	celsius = 35.0 (degC)
	ena = 30
	ek = -75
	el = -40
	SC1 = 1e0
	SC2 = 1e0
	gna = 4.0e-3
	gk = 0.3e-3
	gl = 0.003e-3
	hinit = -1
}
 
STATE {
  h n
}
 
ASSIGNED {
	minf
	hinf
	ninf
	htau
	ntau
	ina
	ik
	il
}
 
BREAKPOINT {
    SOLVE states METHOD cnexp
    ina = gna*pow(minf,3)*h*(v-ena)
    ik = gk*pow(n,4)*(v-ek)
    il = gl*(v-el)
	
}

INITIAL {
        rates(v)
        if(hinit >0){
           h = hinit
        } else {
		   h = hinf
		}
        n = ninf
}

DERIVATIVE states {  :Computes state variables m, h, and n 
        rates(v)
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
}
 
UNITSOFF

PROCEDURE rates(vf){
	LOCAL 	am,bm,ah,bh,an,bn,vbar
	vbar = 127.0/105.0*vf + 8265.0/105.0
	am = 0.1*alpha(vbar, 50.0, 10.0)
	bm = 4.0*beta(vbar,25.0,18.0)
	ah = 0.07*beta(vbar,25.0,10.0)
	bh = 1.0*boltz(vbar,55.0,10.0)
	bn = 0.125*beta(vbar,45.0,80.0)
	an = 0.01*alpha(vbar,55.0,10.0)
	minf = am/(am+bm)
	hinf = ah/(ah+bh)
	ninf = an/(an+bn)
	htau = SC2*12.5/(ah+bh)
	ntau = SC1*12.5/(an+bn)


}

 
FUNCTION boltz(x,y,z) {
                boltz = 1/(1 + exp(-(x - y)/z))
}

FUNCTION alpha(x,y,z){
    if (0){
		alpha = z+(y-x)/2.0 + pow(y-x,2)/(12.0*z)
	} else {
		alpha = (y-x)/(exp((y-x)/z)-1)
	}
}

FUNCTION beta(x,y,z){
	beta = exp((y-x)/z)
}
UNITSON
COMMENT
double alpha_m(v)
 double v;
 {double vbar;
 vbar = (127.0/105.0)*v + (8265.0/105.0);
 return( (TSF*0.1*(50.0 - vbar))/(exp((50.0 - vbar)/10.0)-1.0));}

double beta_m(v)
 double v;
{double vbar;
vbar = (127.0/105.0)*v + (8265.0/105.0);
return ( TSF*4*exp((25 - vbar)/18));}

double alpha_n(v)
 double v;
{double vbar;
vbar = (127.0/105.0)*v + (8265.0/105.0);
return ( TSF*0.01*(55.0 - vbar)/(exp((55.0 - vbar)/10.0)-1.0)); }

double beta_n(v)
 double v;
{double vbar;
vbar = (127.0/105.0)*v + (8265.0/105.0);
return(TSF*0.125*exp((45.0 - vbar)/80.0)); }

double alpha_h(v)
 double v;
{double vbar;
vbar = (127.0/105.0)*v + (8265.0/105.0);
return ( TSF*0.07*exp((25.0 - vbar)/10.0));} /*20*/

double beta_h(v)
 double v;
{double vbar;
vbar = (127.0/105.0)*v + (8265.0/105.0);
return ( TSF*1/(exp((55.0-vbar)/10.0)+1.0));}

ENDCOMMENT

