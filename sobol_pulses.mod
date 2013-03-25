TITLE Current pulses drawn from a quasi-random sequence  
: Delivers a current pulse pulse using sobol quasi-random phase sampling.
: Uses a frequency estimation method to deliver the pulse at the proper phase.
: Parameters are:
:       delay            - delay (ms)
:       pw               - pulse width (ms)
:       amp              - amplitude of current pulse (nA)
:       spkCount         - number of spikes between pulses
:       F                - target firing frequency if using PI (to turn off use 0)
:       gp               - parameter for the PI controller
:       gi               - parameter for the PI controller
:
: 2012 Theoretical Neurobiology and Neuroengineering, University of Antwerp.
: Joao Couto and Daniele Linaro


VERBATIM

#define SOBOL_MAXBIT 30
#define SOBOL_MAXDIM 6

float sobseq(int init)
{
  int j,k,l;
  unsigned long i,im,ipp;
  static float fac;
  static unsigned long in = 0,ix = 0, *iu[SOBOL_MAXBIT+1];
  static unsigned long mdeg[SOBOL_MAXDIM+1] = {0,1,2,3,3,4,4};
  static unsigned long ip[SOBOL_MAXDIM+1] = {0,0,1,1,2,1,4};
  static unsigned long iv[SOBOL_MAXDIM*SOBOL_MAXBIT+1] = {0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};
  
  if (init) {
    if (iv[1] != 1)
      return -1.0;
    fac = 1.0/(1L << SOBOL_MAXBIT);
    for (j=1,k=0; j<=SOBOL_MAXBIT; j++,k+=SOBOL_MAXDIM)
      iu[j] = &iv[k];
    for (k=1; k<=SOBOL_MAXDIM; k++) {
      for (j=1; j<=mdeg[k]; j++)
	iu[j][k] <<= (SOBOL_MAXBIT-j);
      for (j=mdeg[k]+1; j<=SOBOL_MAXBIT; j++) {
	ipp = ip[k];
	i = iu[j-mdeg[k]][k];
	i ^= (i >> mdeg[k]);
	for (l=mdeg[k]-1; l>=1; l--) {
	  if (ipp & 1)
	    i ^= iu[j-l][k];
	  ipp >>= 1;
	}
	iu[j][k] = i;
      }
    }
    return 0;
  }
  im = in++;
  for (j=1; j<=SOBOL_MAXBIT; j++) {
    if (!(im & 1))
      break;
    im >>= 1;
  } 
  if (j > SOBOL_MAXBIT)
    fprintf(stderr, "SOBOL_MAXBIT (%d) too small in sobseq.\n", SOBOL_MAXBIT);
  im = (j-1)*SOBOL_MAXDIM;
  ix ^= iv[im+1];
  return ix*fac;
}

ENDVERBATIM

UNITS {
    (nA) = (nanoampere)
}

NEURON {
    POINT_PROCESS SobolPulses
    RANGE iPulse,iPI,F,delay,estimatedFreq, pw, amp, spkCount,gp,gi,maxPIcount
    ELECTRODE_CURRENT i
}

ASSIGNED {
    i                   (nA)
    tnext               (ms)
    tLastSpk            (ms)
    estimatedFreq       (1/s)
    count               (1)
    countPI             (1)
    erri                (1)
    errp                (1)
    cval                (nA)
    iPI                 (nA)
    iPulse              (nA)
}

PARAMETER {
    delay       = 1000  (ms)
    pw          = 0.5   (ms)
    amp         = 0.25  (nA)
    spkCount    = 10    (1)
    tau         = 0.01  (s)
    gp          = 0     (1)
    gi          = 0     (1)
    maxPIcount  = 3     (1)
    F           = 30    (1/s)
}   

INITIAL {
    tnext = -1
    VERBATIM
    sobseq(1);
    ENDVERBATIM
    erri = 0
    errp = 0
    cval = 0     (nA) 
    iPulse  = 0     (nA) 
    iPI     = 0     (nA) 
}

NET_RECEIVE(weight) { 
    LOCAL w,isi,frac
    INITIAL {
        count         = 0
        estimatedFreq = 0       : to be used in estimating the time delays
        tLastSpk      = 0
        countPI       = maxPIcount
    }
    if ( tLastSpk>0 ) {
        :Frequency estimate
        isi = (t - tLastSpk)/1000.0  : in s
        w   = exp(-isi/tau)
        estimatedFreq = (1-w)/isi + w*estimatedFreq
        :estimatedFreq = 1.0/isi
        :Count spikes
        count = count + 1
        countPI = countPI + 1
    }
    if (t<tnext){
        tnext = -1
    }
    if (count >= spkCount && t>delay) {
        countPI = 0
        count   = 0
	    frac    = sobseq(0)
        tnext   = t + frac*(1000.0/estimatedFreq)
        at_time(tnext)
	:    printf("%12.3f %12.6f %12.3f %.12.6f\n", t, estimatedFreq, tnext, frac)
    }
    if (countPI>=maxPIcount){
        :Frequency clamp
        errp = F - estimatedFreq
        erri = erri + errp
        cval = (gp*errp + gi*erri)/1000.0
        :printf("cval = %12.3f\n",cval)
    }
    :else {
	:printf("%12.3f %12.6f\n", t, estimatedFreq)
    :}
    tLastSpk = t
}

BREAKPOINT {
    iPI = cval
    iPulse = 0
    if (t>tnext && t<=(tnext+pw)){
        iPulse = amp
        at_time(tnext+pw)
    }
    i=iPI+iPulse
    :printf("%12.3f %12.4f %12.4f\n",t,iPI,iPulse)
}
