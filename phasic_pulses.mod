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

UNITS {
    (nA) = (nanoampere)
}

NEURON {
    POINT_PROCESS PhasicPulses
    RANGE iPI,iPulse,F,delay,estimatedFreq, pw, amp, spkCount,gp,gi,maxPIcount,phase
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
    tau         = 0.05  (s)
    gp          = 0     (1)
    gi          = 0     (1)
    maxPIcount  = 3     (1)
    F           = 30    (1/s)
    phase       = 0.5   (1)
}   

INITIAL {
    tnext   = -1
    erri    = 0
    errp    = 0
    cval    = 0     (nA) 
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
	    frac    = phase 
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
}
:BREAKPOINT {
:   i = cval
:   if (t>tnext && t<=(tnext+pw)){
:       i = cval+amp
:       at_time(tnext+pw)
:   }  
:}
