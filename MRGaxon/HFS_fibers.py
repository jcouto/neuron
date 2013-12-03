#!/usr/bin/env python
from MRGnodeHFS import *
import getopt
import os.path as path
import numpy as np
from neuron import h

__all__=['randomPositions','simulatePopulation','calculateInterNodeLength']

def randomPositions(scale=[1,1,1],N=100,minimum=0,maximum=1,fixed_seed=True,verbose=False):
    ''' Creates a N random fiber locations.
    Positions are relative to maximum allowed. 
    Max allowed depends on the fiber diameter.
    '''
    fixed_seed=True
    if fixed_seed:
        np.random.seed(916792037)
    positions=[]
    for i in range(0,N):
        p,r=extract3coord()
        while r<minimum or r>maximum:
            p,r=extract3coord(scale)
        if verbose:
            print('Found position ' + str(p) + ' at distance ' + str(r))
        positions.append(p)
    positions=np.array(positions)
    positions=positions.reshape(N,3)
    r=np.sqrt(np.sum(positions**2,axis=1))
    idx=r.argsort()
    return positions[idx,:]

def extract3coord(scale=[1,1,1]):
    '''
    Returns a randomized 3 coordinate position and the distance.
    '''
    p=np.array([np.random.uniform(0,scale[0],1)[0],
                np.random.uniform(0,scale[1],1)[0],
                np.random.uniform(0,scale[2],1)[0]])
    r=np.sqrt(np.sum(p**2,axis=0))
    return p,r

def calculateInterNodeLength(fiberD=5.7):
    '''
    Returns the inter-node distance in microms given the fiber diameter.
    '''
    paralength1=3  
    nodelength=1.0
    if (fiberD==5.7):
        deltax=500 
        paralength2=35 
    elif (fiberD==7.3):
        deltax=750 
        paralength2=38 
    elif (fiberD==8.7):
        deltax=1000
        paralength2=40 
    elif (fiberD==10.0): 
        deltax=1150 
        paralength2=46 
    elif (fiberD==11.5):
        deltax=1250
        paralength2=50
    elif (fiberD==12.8):
        deltax=1350
        paralength2=54
    elif (fiberD==14.0):
        deltax=1400
        paralength2=56
    elif (fiberD==15.0):
        deltax=1450
        paralength2=58
    elif (fiberD==16.0):
        deltax=1500
        paralength2=60 
    else:
        print('Fiber size not valid.')
        return False
    interlength = (deltax-nodelength-(2*paralength1)-(2*paralength2))/6
    return 2*paralength1+2*paralength2+6*interlength+nodelength

def simulatePopulation(par,recpar,rec,XYZ=None, fixed_seed=True,verbose=False, singleRun = None):
    ''' Runs a simulation on the specified fiber positions XYZ, if None is specified, the default parameters are used (randomized with fixed random seed).
    par,recpar and rec are the simulation parameters, recording parameters and recording vectors respectivelly.
    '''
    scale = np.array([3000,0.5*calculateInterNodeLength(par['fiberD']),3000])
    if XYZ is None:
        XYZ = randomPositions(scale,N=100,minimum=635,maximum=3000,fixed_seed=fixed_seed)
    fiber=[]
    if singleRun is None:
        for ii in range(0,XYZ.shape[0]):
            if verbose:
                print('Running fiber ' + str(ii)  + ' - ' + str(XYZ[ii,:]))
            par['HFSx']=XYZ[ii,0]
            par['HFSy']=XYZ[ii,1]
            par['HFSz']=XYZ[ii,2]
            updateMRGaxon(par,False)
            runMRGaxon()
            if recpar['record']:
                #gname=str(par['HFSfrequency'])+'Hz'
                append_fiber_to_file(rec,par,recpar)
            resetRecorder(rec,False)
    else:
        if singleRun >= 100:
            print('Fiber (%d) exceeds limits (100).' % (singleRun))
            return
        if verbose:
            print('Running fiber ' + str(singleRun)  + ' - ' + str(XYZ[singleRun,:]))
        par['HFSx']=XYZ[singleRun,0]
        par['HFSy']=XYZ[singleRun,1]
        par['HFSz']=XYZ[singleRun,2]
        updateMRGaxon(par,False)
        runMRGaxon()
        if recpar['record']:
            append_fiber_to_file(rec,par,recpar)
        resetRecorder(rec,False)

def main():
    '''
    Runs 100 fibers in series at a specified condition.
    Or a single particular fiber if a separate parameter is specified.
    The first parameter is always the name of the cfg file.  
    '''
    verbose = False
    verbose_level1 = True
    plot = False
    if verbose: 
        print('Running HFS simulation.')
    if plot:
        import pylab as plt
    counter = 0
    for v in sys.argv:
        counter+=1
        if path.basename(__file__) in v:
            break
    try: 
        filename = sys.argv[counter]
    except:
        print('This function takes the filename of the configuration'
              ' file as input.')
        sys.exit(1)
    try:
        singleRun = int(sys.argv[counter + 1])
    except:
        print('Running an entire population.') 
        singleRun =  None
    if not path.exists(filename):
        print('File (%s) does not exist.' % (filename))
        sys.exit(1)
    par, recpar = readConfigurations(filename)
    createMRGaxon(par,verbose)
    rec = recordMRGaxon(recpar,verbose)
    simulatePopulation(par,recpar,rec,None,True,True, singleRun)
    h.quit()

if __name__=='__main__':
    main()
