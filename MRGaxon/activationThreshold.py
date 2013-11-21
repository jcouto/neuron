#!/usr/bin/env python
from MRGnodeHFS import *
import getopt
import os.path as path
import numpy as np

__all__=['calculateStimStep','searchForThreshold']

def calculateStimStep(A,Apast,stimStep, resolution = 0.9e-3, verbose = False):
    '''
    "A" and "Apast" are the current trial and the previous trial results, respectively.
    Implemetation of the binomial search.
        - If the last two trials are the different cut step in half.
        - If there was a block and the step is smaller than the resolution, return nan to stop algorithm.
        - Evaluate the signal of the step depending on the current trial.
    Returns the step to be added to the previous trial stimulation amplitude.
    '''

    stimStep = np.abs(stimStep)
    if not (A==Apast):
        stimStep = stimStep/2.0
        if verbose:
            print('Last two trials were different.Cutting in half.')
    if stimStep<resolution and A:
        if verbose:
            print "Condiction satisfied. Reached last step."
        return np.nan
    if A:
        return -stimStep
    else:
        return stimStep

def searchForThreshold(par,recpar,rec,k,startpoint,resolution=1e-3,
                       caseExpr='len(spks)>0',verbose=False):
    '''
    Performs the binomial search.
    par and recpar are the simulation parameters, k the key to be replaced, and 
    startpoint the value to start the search. resolution specifies when to stop 
    the search and caseExpr is the expression to validate (default is len(spks)>0).
    '''
    pastCase = False
    stimAmp = startpoint
    stimStep = startpoint
    lastNodeName = 'spk' + str(np.max(recpar['nodes']))
    while not np.isnan(stimAmp):
        # Runs the simulation and records data when the condition is satisfied.
        resetRecorder(rec,False)
        par[k] = stimAmp
        if verbose:
            print('Running '+str(stimAmp) + ' on '+ k +
                  ' for condition ' + caseExpr)
        updateMRGaxon(par,False)
        runMRGaxon()
        spks = np.array(rec['spiketimes'][lastNodeName])
        case = eval(caseExpr)
        stimStep = calculateStimStep(case,pastCase,stimStep,
                                         resolution,verbose)
        pastCase = case
        stimAmp = stimAmp + stimStep
    print('Condition satisfied for minimal amplitude at %fmA.'%par[k])
    if recpar['record']:
        append_fiber_to_file(rec,par,recpar)

def main():
    '''
    Runs a binomial search for activation threshold.
    '''
    print('Binomial search for activation threshold.')
    
    verbose = False
    verbose_level1 = True
    counter = 0
    for v in sys.argv:
        counter+=1
        if path.basename(__file__) in v:
            break
    filename = sys.argv[counter]
    par, recpar = readConfigurations(filename)
    createMRGaxon(par,verbose)
    rec = recordMRGaxon(recpar,verbose)
    frequencies = np.array([100,200,500,2000,4000,6000,9000,15000,25000])
    for freq in frequencies:
        par['HFSfrequency'] = freq
        par['HFSamp'] = 0.5
        par['HFSdur'] = 1000.0/freq
        print('Processing frequency : %4.1f'%(freq))
        resetRecorder(rec)
        stimStep = 1
        searchForThreshold(par,recpar,rec,k='HFSamp',startpoint=stimStep,
                           resolution=1e-3,verbose=verbose_level1)
if __name__=='__main__':
    main()
    
