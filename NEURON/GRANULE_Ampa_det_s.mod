TITLE 

COMMENT
ENDCOMMENT

NEURON {
	POINT_PROCESS GRANULE_Ampa_det_s
	NONSPECIFIC_CURRENT i
	RANGE Q10_diff,Q10_channel
	RANGE R, g, ic
	RANGE Cdur, Erev 
	RANGE r1FIX, r2, r3,r4,r5,gmax,r1,r6,r6FIX,kB
	RANGE tau_1, tau_rec, tau_facil, U	 
	RANGE PRE,T,Tmax
		
	RANGE diffuse,Trelease,lamd
	RANGE M,Diff,Rd
	
	RANGE tspike
	RANGE nd, syntype, gmax_factor
}

UNITS {
	(nA) = (nanoamp)	
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(pS) = (picosiemens)
	(nS) = (nanosiemens)
	(um) = (micrometer)
	PI	= (pi)		(1)
    }
    
    PARAMETER {
	syntype
	gmax_factor = 1
	Q10_diff	= 1.1
	Q10_channel	= 2.4
	: Parametri Postsinaptici
	gmax		= 1200   (pS)	

	U 			= 0.42 (1) 	< 0, 1 >
	tau_rec 	= 35.1 (ms) 	< 1e-9, 1e9 > 	 
	tau_facil 	= 10.8 (ms) 	< 0, 1e9 > 	

	M			= 21.515	: numero di (kilo) molecole in una vescicola		
	Rd			= 1.03 (um)
	Diff		= 0.223 (um2/ms)
		 
	Cdur		= 0.3	(ms)		 
	r1FIX		= 5.4		(/ms/mM) 	 				
	r2			= 0.82	(/ms)			 
	r3			= 0		(/ms)		 
	r4			= 0		(/ms)		 
	r5			= 0.013	(/ms)			 
	r6FIX		= 1.12	(/ms/mM)		 
	Erev		= 0	(mV)
	kB			= 0.44	(mM)

	: Parametri Presinaptici
	tau_1 		= 3 (ms) 	< 1e-9, 1e9 >
	 

	
	u0 			= 0 (1) 	< 0, 1 >	: se u0=0 al primo colpo y=U
	Tmax		= 1  (mM)	
	
	: Diffusion			
	diffuse		= 1
	
	lamd		= 20 (nm)
	nd			= 1
	celsius (degC)
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	ic 		(nA)		: current = g*(v - Erev)
	g 		(pS)		: conductance
	r1		(/ms)
	r6		(/ms)
	T		(mM)

	Trelease	(mM)
	tspike[101]	(ms)
	x 
	tsyn		(ms)
	PRE[101]
	
	Mres		(mM)	
	numpulses
	t0
	gbar_Q10 (mho/cm2)
	Q10 (1)
}

STATE {	
	C
	O
	D
}

INITIAL {
	C=1
	O=0
	D=0
	T=0 (mM)
	numpulses=0
	Trelease=0 (mM)
	tspike[0]=1e12	(ms)
	 
	gbar_Q10 = Q10_diff^((celsius-30)/10)
	Q10 = Q10_channel^((celsius-30)/10)

	: fattore di conversione che comprende molecole -> mM
	: n molecole/(Na*V) -> M/(6.022e23*1dm^3)

	Mres = 1e3 * ( 1e3 * 1e15 / 6.022e23 * M ) : (M) to (mM) so 1e3, 1um^3=1dm^3*1e-15 so 1e15	
	numpulses=0

	FROM i=1 TO 101 { PRE[i-1]=0 tspike[i-1]=0 }  
	tspike[0]=1e12	(ms)
	if(tau_1>=tau_rec){ 
		printf("Warning: tau_1 (%g) should never be higher neither equal to tau_rec (%g)!\n",tau_1,tau_rec)
		tau_rec=tau_1+1e-5
	} 
}

FUNCTION diffusione(){	 
	LOCAL DifWave,i	
	DifWave=0
	FROM i=1 TO numpulses{
		t0=tspike[i-1]
		if(t>t0){
			DifWave=DifWave+PRE[i-1]*Mres*exp(-Rd*Rd/(4*Diff*(t-t0)))/((4*PI*Diff*(1e-3)*lamd)*(t-t0))^nd
		}
	}	
	diffusione=DifWave
}

BREAKPOINT {

	if ( diffuse && (t>tspike[0]) ) { Trelease= T + diffusione() } else { Trelease=T }

	SOLVE kstates METHOD sparse
	g = gmax * gbar_Q10 * O
	i = (1e-6) * g * (v - Erev) * gmax_factor
	ic = i
}

KINETIC kstates {
	r1 = r1FIX * Trelease^2 / (Trelease + kB)^2
	r6 = r6FIX * Trelease^2 / (Trelease + kB)^2
	~ C  <-> O	(r1*Q10,r2*Q10)
	~ O  <-> D	(r3*Q10,r4*Q10)
	~ D  <-> C	(r5*Q10,r6*Q10)
	CONSERVE C+O+D = 1
}


NET_RECEIVE(weight, on, nspike, t0 (ms),y, z, u, tsyn (ms)) {

: *********** ATTENZIONE! ***********
:
: Qualora si vogliano utilizzare impulsi di glutammato saturanti e'
: necessario che il pulse sia piu' corto dell'intera simulazione
: altrimenti la variabile on non torna al suo valore di default.

INITIAL {
	y = 0
	z = 0
	u = u0
	tsyn = t
	nspike = 1
}
   if (flag == 0) { 
		: Qui faccio rientrare la modulazione presinaptica
		nspike = nspike + 1
		if (!on) {
			t0 = t
			on = 1				
			z = z*exp( - (t - tsyn) / (tau_rec) )	 
			z = z + ( y*(exp(-(t - tsyn)/tau_1) - exp(-(t - tsyn)/(tau_rec)))/((tau_1/(tau_rec))-1) )  
			y = y*exp(-(t - tsyn)/tau_1)			
			x = 1-y-z
				
			if (tau_facil > 0) { 
				u = u*exp(-(t - tsyn)/tau_facil)
				u = u + U * ( 1 - u )							
			} else { u = U }
			
			y = y + x * u
			
			T=Tmax*y
			PRE[numpulses]=y	: PRE[numpulses]=y
			
			:PRE=1	: Istruzione non necessaria ma se ommesso allora le copie dell'oggetto successive alla prima non funzionano!
			:}
			: all'inizio numpulses=0 !			
			
			tspike[numpulses] = t
			numpulses=numpulses+1
			tsyn = t
			
		}
		net_send(Cdur, nspike)	 
    }
    if (flag == nspike) { 
		t0 = t
		T = 0
		on = 0
    }
}

