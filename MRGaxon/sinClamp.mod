COMMENT
Adapted from the Neuron Forum...
ENDCOMMENT

NEURON {
        POINT_PROCESS SinClamp
        RANGE del, dur, pkamp, freq, phase, bias
        ELECTRODE_CURRENT i
}

UNITS {
        (nA) = (nanoamp)
             }

PARAMETER {
        del=5   (ms)
        dur=200   (ms)
        pkamp=1 (nA)
        freq=1  (Hz)
        phase=0 (rad)
        bias=0  (nA)
        PI=3.14159265358979323846
}

ASSIGNED {
        i (nA)
}

BREAKPOINT {
       if (t < del) {
      i=0   
   }else{  
            if (t < del+dur) {
           i = pkamp*sin(2*PI*freq*(t-del)/1000+phase)+bias
      }else{  
           i = 0
}}}
 

