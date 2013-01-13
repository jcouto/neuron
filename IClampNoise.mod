TITLE Normal distribution

COMMENT
-----------------------------------------------------------------------------

    Point current with noise coming from a normal distribution
    ==================================================

IMPLEMENTATION

  This mechanism is implemented as a nonspecific current defined as a point process.

PARAMETERS

  The mechanism takes as input the following parameters:

     del (ms)
     dur (ms)
     qmp (nA)

-----------------------------------------------------------------------------
ENDCOMMENT


NEURON {
  POINT_PROCESS IClampNoise
  RANGE i,del,dur,amp
  ELECTRODE_CURRENT i
}

UNITS {
  (nA) = (nanoamp)
}

PARAMETER {
  del=500    (ms)
  dur=1000000   (ms)
  amp=0 (nA)
}

ASSIGNED {
  ival (nA)
  i (nA)
  noise (nA)
  on (1)
}

INITIAL {
  i = 0
  on = 0
  net_send(del, 1)
}

PROCEDURE seed(x) {
  set_seed(x)
}

BEFORE BREAKPOINT {
  if  (on) {
    noise = normrand(0,amp*1(/nA))*1(nA)
    ival = noise
  } else {
    ival = 0
  }
}

BREAKPOINT {
  i = ival
}

NET_RECEIVE (w) {
  if (flag == 1) {
    if (on == 0) {
      : turn it on
      on = 1
      : prepare to turn it off
      net_send(dur, 1)
    } else {
      : turn it off
      on = 0
    }
  }
}
