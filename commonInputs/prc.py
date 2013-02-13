#!/usr/bin/env python

from CommonInput import *
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

def makeOutputFilename(prefix='', extension='.out'):
    filename = prefix
    if prefix != '' and prefix[-1] != '_':
        filename = filename + '_'
    now = time.localtime(time.time())
    filename = filename + '%d%02d%02d-%02d%02d%02d' % \
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
    available_models = ['kr']
    try:
        opts,args = getopt.getopt(sys.argv[1:],'hi:m:g:f:d:s:a:w:n:',['help'])
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
    for o,a in opts:
        if o == '-h':
            usage()
            sys.exit(0)
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
                
    N = 1500
    if model == 'kr':
        neurons = [KhaliqRaman(gid,
                               {'diameter':cell_size,'length':cell_size,'nSynapses':0},
                               prcProps={'frequency':frequency,
                                         'gp':0.01, 'gi':0.1,
                                         'pulseAmp':pulse_amp,
                                         'pulseWidth':pulse_width,
                                         'spkCount':spk_count},) for i in range(N)]
    for n in neurons:
        n._addPRCestimator()

    # Uncomment to record voltage
    #for n in neurons:
    #    n.addSomaticVoltageRecorder()
    #time = h.Vector()
    #time.record(h._ref_t)
    
    h.load_file('stdrun.hoc')
    if duration is None:
        duration = N*spk_count*1./frequency*1000.
    h.tstop = duration
    h.run()
    print('Computed prc for frequency %3.0f.'%(f))
    
    for n in neurons:
        saveNeuron(filename, n)

    # Uncomment to plot
    #fig = plt.figure()
    #for i,n in enumerate(neurons):
    #    fig.add_subplot(len(neurons),1,i)
    #    plt.plot(time,n.somaticVoltage(),'k')
    #    plt.plot(n.perturbationTimes(),np.ones(np.shape(n.perturbationTimes()))*20.,'|r')
    #    plt.plot(n.spikeTimes(),np.ones(np.shape(n.spikeTimes()))*20.,'|b')
    #plt.show()


if __name__ == '__main__':
    main()
