: $Id: pregen2.mod 1.2 Wed, 02 Feb 2005 16:10:51 +0100 sergio $
: comments at end

NEURON	{ 
  POINT_PROCESS Golgi_SpikeGenerator2
  RANGE y
  RANGE fast_invl, slow_invl, burst_len, start, end, delay,on
  RANGE noise,time
}

PARAMETER {
	:burst_len	= 40		: burst length (# spikes)
	start		= 100000 (ms)	: start of first interburst interval
	end		= 1e10 (ms)	: time to stop bursting
	:noise		= 1		: amount of randomeaness (0.0 - 1.0)
	:delay		= 4 (ms)	: to avoid 
}

ASSIGNED {
	y
	:burst
	event (ms)
	:burst_off (ms)
	:burst_on (ms)
	toff (ms)
	on

	time[5000] (ms)	
	indice
}

:PROCEDURE inits(){
:}

INITIAL {	
	:init()
	indice=0
	on = 1
	toff = 1e9
	y = -90
	:burst = 0
	event = start 
	event_time()
	while (on == 1 && event < 0) {
		event_time()
	}
	if (on == 1) {
		net_send(event, 1)
	}
}	

PROCEDURE event_time() {
	event=event+time[indice] 
	indice=indice+1
	if (event > end) {
		on = 0
	}
}

NET_RECEIVE (w) {
:printf("Pregen receive t=%g flag=%g\n", t, flag) 
	if (flag == 1 && on == 1) {
		y = 20
		net_event(t)
		event_time()
		net_send(event - t, 1)
		net_send(.1, 2)
	}
	if (flag == 2) {
		y = -90
	}
}

COMMENT
Presynaptic spike generator
---------------------------

This mechanism has been written to be able to use synapses in a single
neuron receiving various types of presynaptic trains.  This is a "fake"
presynaptic compartment containing a fast spike generator.  The trains
of spikes can be either periodic or noisy (Poisson-distributed), and 
either tonic or bursting.

Parameters;
   noise: 	between 0 (no noise-periodic) and 1 (fully noisy)
   fast_invl: 	fast interval, mean time between spikes (ms)
   slow_invl:	slow interval, mean burst silent period (ms), 0=tonic train
   burst_len: 	mean burst length (nb. spikes)

Written by Z. Mainen, modified by A. Destexhe, The Salk Institute

Modified by Michael Hines for use with CVode

Modified by Michael Hines to use logical event style with NET_RECEIVE
ENDCOMMENT

