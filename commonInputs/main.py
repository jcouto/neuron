#!/usr/bin/env python

from CommonInput import *
from neuron import h

def main():
    fixedInput = {'probability': 0.5, 'spikeTimes': [100,215]}

    neurons = [KhaliqRaman(i) for i in range(2)]
    for n in neurons:
        n.addFixedInput(fixedInput['probability'], fixedInput['spikeTimes'])
        n.addPoissonInputs(clusterSize=5, stimulusProps={'frequency': 0.1, 'noise': 1.})
        n.addSomaticVoltageRecorder()

    time = h.Vector()
    time.record(h._ref_t)

    h.load_file('stdrun.hoc')
    h.tstop = 400
    h.run()
    for n in neurons:
        #print('Neuron [%d] emitted %d spikes.' % (n.ID, len(n.spikeTimes)))
        saveNeuron('neuron.h5', n)

if __name__ == '__main__':
    main()
