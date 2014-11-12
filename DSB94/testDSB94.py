#!/usr/bin/env python
from neuron import h
h.load_file('DSB94tmplt')
h.load_file('stdrun.hoc')
import numpy as np
import pylab as plt

def main():
    pkj = h.DSB94()
    h.dt = 0.02
    Vrest = -74
    h.celsius = 37
    rec={}
    rec['v'] = h.Vector()
    rec['t'] = h.Vector()
    rec['v'].record(pkj.soma(0.5)._ref_v)  
    rec['t'].record(h._ref_t)
    iclamp = h.IClamp(pkj.soma(0.5))
    iclamp.amp = 0.300
    iclamp.dur = 1000
    iclamp.delay = 0
    h.finitialize(Vrest)
    h.tstop = 700
    h.run()
    plt.plot(rec['t'],rec['v'])
    plt.show()
    from ipdb import set_trace
    set_trace()
if __name__ == '__main__':
    main()
