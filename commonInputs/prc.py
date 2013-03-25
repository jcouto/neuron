#!/usr/bin/env python

from CommonInput import *
from DSB94 import *
import numpy as np
from neuron import h
import pylab as plt
import getopt
import sys
import time
import os

def usage():
    """
    Prints a help message.
    """
    print('++++++++++++++++++++++++++++++++++++++++++++++')
    print('Script to calculate the phase response curve.')
    print('Usage: %s [option <value>]'% sys.argv[0])
    print('Available options:')
    print('\t -h : display help message')
    print('\t -i : global id, used for recording.')
    print('\t -g : filename')
    print('\t -m : model [kr (Khaliq-Raman - default)|not implemented]')
    print('\t -f : frequency at which to compute the prc [30]')
    print('\t -s : cell diameter in micrometers [20]')
    print('\t -d : duration of the recording. If not specified ~ 1500 trials are used.')
    print('\t -a : amplitude of the perturbation pulse')
    print('\t -w : pulse width of the perturbation pulse')
    print('\t -n : amplitude of the noise (0pA - default [Not implemented!])')
    print('\t -p : plot (on/off)')

def makeOutputFilename(prefix='', extension='.out'):
    filename = prefix
    if prefix != '' and prefix[-1] != '_':
        filename = filename + '_'
    now = time.localtime(time.time())
    filename = filename + '%d%02d%02d%02d%02d%02d' % \
        (now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec)
    if extension[0] != '.':
        extension = '.' + extension
    suffix = ''
    k = 0
    while os.path.exists(filename + suffix + extension):
        k = k+1
        suffix = '_%d' % k
    return filename + suffix + extension

def main():
    available_models = ['kr','dsb94']
    try:
        opts,args = getopt.getopt(sys.argv[1:],'hpi:m:g:f:d:s:a:w:n:',['help'])
    except getopt.GetoptError, err:
        print('Could not complete task!')
        print(err)
        usage()
        sys.exit(1)
    # Default parameters
    gid  = 0
    model = 'kr'
    frequency = 30.0
    filename = makeOutputFilename(prefix='', extension='.h5')
    duration = None
    cell_size = 20
    pulse_amp = 0.05
    pulse_width = 1
    noise_amp = 0
    spk_count = 6
    plot = False
    for o,a in opts:
        if o == '-h':
            usage()
            sys.exit(0)
        elif o=='-p':
            plot = True
        elif o == '-m':
            if not a.lower() in available_models:
                print('Model %s not available.' % a)
                sys.exit(1)
            else:
                model = a.lower()
        elif o == '-g':
            filename = a
        elif o == '-i':
            gid = float(a)
        elif o == '-f':
            frequency = float(a)
        elif o == '-d':
            duration = float(a)
        elif o == '-s':
            cell_size = float(a)
        elif o == '-a':
            pulse_amp = float(a)
        elif o == '-w':
            pulse_width = float(a)
        elif o == '-n':
            noise_amp = float(a)
    ntrials = 1500
    N = 1
    if model == 'kr':
        neurons = [KhaliqRaman(gid,
                               {'diameter':cell_size,'length':cell_size,'nSynapses':0},
                               prcProps={'frequency':frequency,
                                         'gp':0.01, 'gi':0.1,
                                         'pulseAmp':pulse_amp,
                                         'pulseWidth':pulse_width,
                                         'spkCount':spk_count,
                                         'delay':1000}) 
                   for gid in range(N)]
        for n in neurons:
            n.soma.gkbar_hpkj = n.soma.gkbar_hpkj * 1.40
    elif model == 'dsb94':
        #fixedInput = {'probability': 0.5, 'spikeTimes': np.cumsum(np.random.poisson(30,140))}
        #for n in neurons:
        #    n.addFixedInput(fixedInput['probability'], fixedInput['spikeTimes'])

        neurons = [DSB94(gid,
                               {'nSynapses':0},
                               prcProps={'frequency':frequency,
                                         'gp':0, 'gi':0,
                                         'pulseAmp':pulse_amp,
                                         'pulseWidth':pulse_width,
                                         'spkCount':spk_count}) 
                   for gid in range(N)]
        baseline = []
        for n in neurons:
            baseline.append(h.IClamp(n.soma(0.5)))
            baseline[-1].amp = 0.0
            baseline[-1].dur = duration
            baseline[-1].delay = 10
 
    for n in neurons:
        n._addPRCestimator()

    if plot:
        iPID = []
        iPulse = []
        estimatedFreq = []
        for n in neurons:
            n.addSomaticVoltageRecorder()
            iPID.append(h.Vector())
            iPulse.append(h.Vector())
            estimatedFreq.append(h.Vector())
            iPID[-1].record(n._prc_sobol._ref_iPI)
            iPulse[-1].record(n._prc_sobol._ref_iPulse)
            estimatedFreq[-1].record(n._prc_sobol._ref_estimatedFreq)
        time = h.Vector()
        time.record(h._ref_t)

    h.load_file('stdrun.hoc')
    if duration is None:
        duration = np.ceil(ntrials*spk_count/frequency*1000.0)
    print('Going to run for %dms'%duration)
    h.tstop = duration
    h.run()
    print neurons[0].perturbationTimes()
    print('Computed prc for frequency %3.0fHz.'%(frequency))
    for n in neurons:
        saveNeuron(filename, n)
    
    if plot:
        fig = plt.figure()
        for i,n in enumerate(neurons):
            fig.add_subplot(len(neurons),1,i)
            plt.plot(time,n.somaticVoltage(),'k')
            plt.plot(time,np.array(iPID[i])*1000,'r')
            plt.plot(time,np.array(iPulse[i])*1000,'b')
            plt.plot(time,np.array(estimatedFreq[i]),'g')
            plt.plot(n.perturbationTimes(),np.ones(np.shape(n.perturbationTimes()))*20.,'|r')
            plt.plot(n.spikeTimes(),np.ones(np.shape(n.spikeTimes()))*20.,'|b')
        plt.show()

    h.quit()

if __name__ == '__main__':
    main()
