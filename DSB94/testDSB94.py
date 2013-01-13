#!/usr/bin/env python
from neuron import h
h.load_file('DSB94template')
import numpy as np
import pylab as plt

def main():
    pkj = h.DSB94()
    h.dt = 0.02
    Vrest = -50
    h.celsius = 37
    rec={}
    rec['v'] = h.Vector()
    rec['t'] = h.Vector()
    rec['v'].record(pkj.soma(0.5)._ref_v)  
    rec['t'].record(h._ref_t)
    h.load_file('stdrun.hoc')
    vclamp = h.IClamp(pkj.soma(0.5))
    vclamp.amp = 1
    vclamp.dur = 20+np.random.uniform(40)
    vclamp.delay = 0
    h.finitialize(Vrest)
    h.tstop = 400
    h.run()
    plt.plot(rec['t'],rec['v'])
    plt.show()

if __name__ == '__main__':
    main()
