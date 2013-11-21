#!/usr/bin/env python
from MRGnodeHFS import *
import getopt
import os.path as path
import numpy as np
from activationThreshold import *
from HFS_fibers import *

def main():
    '''
    Runs a binomial search for block threshold at HFS.
    To run this use the configuration example in 'template.cfg'
    '''
    print('Binomial search for block threshold.')
    
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
    frequencies = np.append(np.array([3e3]),np.arange(5e3,30e3,2.0e3))
    for freq in frequencies:
        print('Processing frequency : %4.1f'%(freq))
        par['HFSfrequency'] = freq
        par['HFSamp'] = 0.5
        #par['HFSy'] = 0.5*calculateInterNodeLength(par['fiberD'])
        stimStep = 1.0
        resetRecorder(rec)
        searchForThreshold(par,recpar,rec,k='HFSamp',
                           startpoint=stimStep,
                           resolution=1e-3,
                           caseExpr='len(spks[(spks>30.0) & (spks<50.0)])==0',
                           verbose=verbose_level1)
if __name__=='__main__':
    main()

