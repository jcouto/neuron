#!/usr/bin/env python

import sys
import os
import time
import tables as tbl
import h5utils as h5
import numpy as np
from ConfigParser import ConfigParser
from neuron import h
from CommonInput import *

### FUNCTIONS
def printHelp(scriptname):
    print('\nUsage: ' + scriptname + ' [options]')
    print('\twhere options are:');
    print('\t\t-f --conf-file: parse the specified configuration file.')
    print('\t\t-h --help: print this help message.')
    print('Authors: Daniele Linaro -- daniele@tnb.ua.ac.be')
    print('         Joao Couto -- joao@tnb.ua.ac.be\n')

def packTime():
    tm = time.localtime()
    return '%d-%02d-%02d %02d:%02d:%02d' % (tm.tm_year, tm.tm_mon, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec)
    
def logger(message, coreId):
    if coreId == 0:
        print('[' + packTime() + '] ' + message)

###

### PROPERTIES CLASSES
class SimulationProps (tbl.IsDescription):
    dt = tbl.Float64Col()
    tend = tbl.Float64Col()
    nNeurons = tbl.Int32Col()
    nCores = tbl.Int32Col()

class InputProps (tbl.IsDescription):
    poissonRate = tbl.Float64Col()
    poissonNoise = tbl.Float64Col()
    synapticClusterSize = tbl.Int32Col()
    fixedInputProbability = tbl.Float64Col()
    initialPulseDur = tbl.Float64Col()
    initialPulseDelay = tbl.Float64Col()
    initialPulseAmp = tbl.Float64Col()
###

# figure out which configuration file to use
configFile = 'common_input.cfg'
for k,arg in enumerate(sys.argv):
    if arg == '-f' or arg == '--config-file':
        if k < len(sys.argv)-1:
            configFile = sys.argv[k+1]
        else:
            print('Missing argument to %s option.' % arg)
            h.quit()
if not os.path.isfile(configFile):
    print('%s: no such file.' % configFile)
    h.quit()
    
# NEURON's parallel context used to manage the simulation
pc = h.ParallelContext()

# the id of this core
myId = int(pc.id())

# total number of cores used in the simulation
numberOfCores = int(pc.nhost())

# open the configuration file
fid = open(configFile,'r')
config = ConfigParser()
config.readfp(fid)

### read the configuration file

# properties of the connection between presynaptic neuron and synapse
connectionProps = {'weight': config.getfloat('Synapses','weight'),
                   'delay': config.getfloat('Synapses','delay')}

# properties of the input
inputProps = {'synapticClusterSize': config.getint('Input','synapticClusterSize'),
              'poissonRate': config.getfloat('Input','poissonRate'),
              'poissonNoise': config.getfloat('Input','poissonNoise'),
              'fixedInputProbability': config.getfloat('Input','fixedInputProbability'),
              'initialPulseDur': config.getfloat('Input','initialPulseDur'),
              'initialPulseDelay': config.getfloat('Input','initialPulseDelay'),
              'initialPulseAmp': config.getfloat('Input','initialPulseAmp')}

# properties of the simulation
simulationProps = {'tend': config.getfloat('Simulation','tend'), 'dt': h.dt,
                   'nNeurons': config.getint('Simulation','nNeurons'),
                   'nCores': numberOfCores}
try:
    simulationProps['dt'] = config.getfloat('Simulation','dt')
except:
    pass
try:
    outfile = config.get('Simulation','outfile')
except:
    outfile = 'common_input.h5'

# create the neurons
logger('Building neurons...', myId)
neuronType = config.get('Neuron','type')
if neuronType.lower() == 'khaliqraman':
    # properties of the synapses
    synapseProps = {'name': config.get('Synapses','type'),
                    'Erev': config.getfloat('Synapses','Erev')}
    try:
        synapseProps['tauRise'] = config.getfloat('Synapses','tauRise')
        synapseProps['tauDecay'] = config.getfloat('Synapses','tauDecay')
    except:
        pass
    # properties of the neurons
    neuronProps = {'length': config.getfloat('Neuron','length'),
                   'diameter': config.getfloat('Neuron','diameter'),
                   'nSynapses': config.getint('Neuron','nSynapses')}
    # the neurons
    neurons = [KhaliqRaman(i, neuronProps, synapseProps) 
               for i in range(myId, simulationProps['nNeurons'], numberOfCores)]
else:
    print('Unknown neuron type [%s].' % neuronType)
    sys.exit(1)

# the times at which presynaptic spikes will be delivered
fixedSpikeTimes = [100,215]

# we inject a brief random pulse of current to the soma in order to have
# the neurons firing out of phase.
initialPulses = []   

logger('Adding inputs to neurons...', myId)
for n in neurons:
    #n.addFixedInput(inputProps['fixedInputProbability'], fixedSpikeTimes)
    n.addPoissonInputs(inputProps['synapticClusterSize'],
                       {'frequency': inputProps['poissonRate'],
                        'noise': inputProps['poissonNoise']},
                       connectionProps)
    n.addSomaticVoltageRecorder()
    stim = h.IClamp(n.soma(0.5))
    stim.dur = np.random.uniform(high=inputProps['initialPulseDur'])
    stim.amp = np.random.uniform(low=-inputProps['initialPulseAmp'],high=inputProps['initialPulseAmp'])
    stim.delay = np.random.uniform(high=inputProps['initialPulseDelay'])
    initialPulses.append(stim)

logger('Started the simulation...', myId)
# run the simulation    
h.load_file('stdrun.hoc')
h.tstop = simulationProps['tend']
h.dt = simulationProps['dt']
pc.set_maxstep(1)
h.stdinit()
pc.psolve(simulationProps['tend'])
pc.barrier()
# the simulation is finished.

logger('Saving data...', myId)
# save everything.
if myId == 0:
    fid = h5.H5File(outfile, 'w', 'Common input simulation')
    fid.createGroup('/','Properties')
    fid.writeTable('/Properties', 'Simulation', SimulationProps, 'Simulation properties', simulationProps)
    fid.createGroup('/','Input')
    fid.writeTable('/Input', 'Properties', InputProps, 'Input properties', inputProps)
    fid.writeArray('/Input', 'PresynapticSpikes', tbl.Float64Atom(), fixedSpikeTimes)
    fid.close()
pc.barrier()

for i in range(numberOfCores):
    if myId == i:
        for n in neurons:
            saveNeuron(outfile, n, saveVoltage=True)
    pc.barrier()

logger('Terminating program...', myId)
pc.runworker()
pc.done()
h.quit()


