#!/usr/bin/env python

from CommonInput import *
#from DSB94 import *
import numpy as np
from neuron import h
import getopt
import sys
import time
import os
#import ipdb
import h5py as h5

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
    print('\t -m : model [kr (Khaliq-Raman - default)|kr_cn|not implemented]')
    print('\t -f : frequency at which to compute the prc [30]')
    print('\t -s : cell diameter in micrometers [20]')
    print('\t -d : duration of the recording. If not specified ~ 1500 trials are used.')
    print('\t -a : amplitude of the perturbation pulse')
    print('\t -w : pulse width of the perturbation pulse')
    print('\t -g : Recording filename')
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
    available_models = ['kr','dsb94', 'kr_cn']
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
    saveVoltage = True
    duration = None
    cell_size = 20
    pulse_amp = 0.05
    pulse_width = 1
    noiseAmplitude = 0
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
            noiseAmplitude = float(a)
    ntrials = 5000
    N = 1
    if duration is None:
        duration = np.ceil(ntrials*spk_count/frequency*1000.0)
    if model == 'kr':
        neurons = [KhaliqRaman(gid,
                               {'diameter':cell_size,'length':cell_size,'nSynapses':1000,'channelNoise':False},
                               prcProps={'frequency':frequency,
                                         'gp':0.01, 'gi':0.1,
                                         'pulseAmp':pulse_amp,
                                         'pulseWidth':pulse_width,
                                         'spkCount':spk_count,
                                         'delay':1000}) 
                   for gid in range(N)]
        for n in neurons:
            n.addGaussianCurrent(noiseAmplitude,0,duration)
                
    elif model == 'kr_cn':
        neurons = [KhaliqRaman(gid,
                               {'diameter':cell_size,
                                'length':cell_size,
                                'nSynapses':0, 'channelNoise': True},
                               prcProps={'frequency':frequency,
                                         'gp':0.01, 'gi':0.1,
                                         'pulseAmp':pulse_amp,
                                         'pulseWidth':pulse_width,
                                         'spkCount':spk_count,
                                         'delay':1000}) 
                   for gid in range(N)]
        #for n in neurons:
        #    n.soma.gkbar_hpkj = n.soma.gkbar_hpkj * 1.40
    elif model == 'dsb94':
        
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

    if saveVoltage:
       # iPID = []
        iPulse = []
        #estimatedFreq = []
#       gaussCurrent = []
        for i,n in enumerate(neurons):
            n.addSomaticVoltageRecorder()
 #           iPID.append(h.Vector())
            iPulse.append(h.Vector())
#            estimatedFreq.append(h.Vector())
#           gaussCurrent.append(h.Vector())
#            iPID[-1].record(n._prc_sobol._ref_iPI)
            iPulse[i].record(n._prc_sobol._ref_iPulse)
#           gaussCurrent[-1].record(n._IClampNoise._ref_i)
#            estimatedFreq[-1].record(n._prc_sobol._ref_estimatedFreq)
        time = h.Vector()
        time.record(h._ref_t)

    h.load_file('stdrun.hoc')
    print('Going to run for %dms'%duration)
    h.tstop = duration
    h.run()
#    print neurons[0].perturbationTimes()
    print('Computed prc for frequency %3.0fHz.'%(frequency))
    downsampleFactor=4
    for ii,n in enumerate(neurons):
        saveNeuron(filename, n, saveVoltage=False, simulationName='PRC %s'%(model))
        fid = h5.File(filename, 'a')
        groupName = '/Neurons/ID_' + str(n.ID)
        group = fid[groupName]
        ds = group.create_dataset('Perturbation',data= np.array(iPulse[ii])[::downsampleFactor])
        dt = np.unique(np.round(np.diff(np.array(time)[::downsampleFactor]),3))[0]
        tend = np.array(time)[::downsampleFactor][-1]
        ds.attrs['dt'] = dt
        ds.attrs['tend'] = tend

        ds = group.create_dataset('Voltage', data = n.somaticVoltage()[::downsampleFactor])
        ds.attrs['dt'] = dt
        ds.attrs['tend'] = tend
        fid.close()

    if plot:
        import pylab as plt
        fig = plt.figure()
        for i,n in enumerate(neurons):
            fig.add_subplot(len(neurons),1,i)
            plt.plot(time,n.somaticVoltage(),'k')
#           plt.plot(time,np.array(gaussCurrent[i])*1000.0,'k')
        #    plt.plot(time,np.array(iPID[i])*1000,'r')
            plt.plot(time,np.array(iPulse[i])*1000,'b')
         #   plt.plot(time,np.array(estimatedFreq[i]),'g')
         #   plt.plot(n.perturbationTimes(),np.ones(np.shape(n.perturbationTimes()))*20.,'|r')
         #   plt.plot(n.spikeTimes(),np.ones(np.shape(n.spikeTimes()))*20.,'|b')
            isi = np.diff(n.spikeTimes())
            print('The CV is %f.'%(np.std(isi)/np.mean(isi)))
        plt.show()
#    ipdb.set_trace()
    h.quit()

if __name__ == '__main__':
    main()
