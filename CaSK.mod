
TITLE  calcium and SK
 
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
        SUFFIX cask
        USEION ca WRITE ica
        USEION k WRITE ik
        RANGE ica, ik, sinf, shalf, sslope, kc, RHO,gca,gsk,stau,cinit,sinit,km
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
    v   (mV)
    dt  (ms)
	celsius = 35.0 (degC)
	gca = 0.004e-3 (S/cm2)
	gsk = 0.03e-3 (S/cm2)
	alpha = 1
	eca = 140.0
	ek = -75.0
	RHO = 1.5e-4 : -4+3
	kc = 0.00425
	km = 0.5
	shalf = -40
	sslope = 3.333 :1/0.3
	stau = 9400.0
	cinit = 0.1
	sinit = 0
}
 
STATE {
  c
  s
}
 
ASSIGNED {
	ica
	ik
	:stau
	sinf
}
 
BREAKPOINT {
    SOLVE states METHOD cnexp
	ica = gca*s*(v-eca)
	ik = gsk*(v-ek)*c/(c+km)
	
}

INITIAL {
        rates(v)
        c = cinit
        if (sinit > 0 && sinit < 1){
          s = sinit
        } else {
		  s = sinf
		}
        
}

DERIVATIVE states {  :Computes state variables m, h, and n 
        rates(v)
        s' = 1.0*(sinf-s)/stau
        c' = RHO*(kc*s*(eca-v)-c)
}
 
UNITSOFF

PROCEDURE rates(vf){
	sinf = boltz(vf, shalf, sslope)
}

 
FUNCTION boltz(x,y,z) {
                boltz = 1/(1 + exp(-(x - y)/z))
}

UNITSON

