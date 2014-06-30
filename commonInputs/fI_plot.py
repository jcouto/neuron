#!/usr/bin/env python

from CommonInput import *
#from DSB94 import *
import numpy as np
from neuron import h
import getopt
import sys
import time
import os
import ipdb
import h5py as h5


def main():
    gid  = 0
    model = 'kr'
    cell_size = 120
    N = 20
    neurons = [KhaliqRaman(gid,
                           {'diameter':cell_size,'length':cell_size,
                            'nSynapses':0,'channelNoise':False})
               for gid in range(N)]
    baseline = []
    pulse = []
    for n,current in zip(neurons,np.linspace(0.04,0.07,N)):
        baseline.append(h.IClamp(n.soma(0.5)))
        baseline[-1].amp = -0.17
        baseline[-1].dur = 2000
        baseline[-1].delay = 0
        pulse.append(h.IClamp(n.soma(0.5)))
        pulse[-1].amp = current
        pulse[-1].dur = 1000
        pulse[-1].delay = 500
 
    for i,n in enumerate(neurons):
        n.addSomaticVoltageRecorder()
    time = h.Vector()
    time.record(h._ref_t)

    h.load_file('stdrun.hoc')
    h.tstop = 2000
    h.run()
    
    import pylab as plt
    fig = plt.figure()
    for i,n in enumerate(neurons):
#        fig.add_subplot(len(neurons),1,i)
        plt.plot(time,n.somaticVoltage()+(i-1)*5)
    plt.show()
    ipdb.set_trace()

    h.quit()

if __name__ == '__main__':
    main()
