#!/usr/bin/env python
from neuron import h
import h5py as h5
import os.path as path
import numpy as np
import pylab as plt
import getopt
import ConfigParser as cp
import time
######################## Default parameters ##############################
#MODEL PARAMETERS
g_par                       = {}               # global parameters
g_par["channelDescription"] = 0
g_par["axonnodes"]          = 51
g_par["fiberD"]             = 10
#SIMULATION PARAMETERS
g_par["dt"]                 = 0.002
g_par["tstop"]              = 50
#HFS PARAMETERS
g_par["HFSreferenceNode"]   = 25
g_par["HFSdur"]             = 50
g_par["HFSfrequency"]       = 5e3	
g_par["HFSpolarity"]        = 1.0
g_par["HFSdelay"]           = 0
g_par["HFSpulsewidth"]      = 0.09
g_par["HFSamp"]	            = 1.154		#mA
g_par["HFSwaveform"]        = 0			#1 sinusoid
g_par["HFSx"]               = 0			# all positions in micrometers
g_par["HFSy"]               = 0
g_par["HFSz"]               = 1000.0
#INTRACELLULAR PARAMETERS
g_par["intrinsicNode"]      = 0
g_par["intrinsicDur"]       = 0.1 
g_par["intrinsicAmp"]       = 1
# for this simulation I need only that it is delivered at t=40ms, it is loaded after MRGnodeHFS because of the objrefs and new Vectors...
g_vecpar                    = {}                # default intracellular pattern
g_vecpar["pattern"]	    = np.array([40.0])
g_vecpar["patternLag"]	    = np.array([ii-0.005 for ii in g_vecpar["pattern"]])
#RECORD DATA
g_recpar		    = {}                # parameters to record data
g_recpar["record"]	    = False
g_recpar["filename"]        = "./data/simulation.h5"

##########################################################################
######################## Parse cfg files #################################
def getOptions(cfg,metadata,section,verb=False):
    '''
    Parse configuration file from default values.
       - "cfg" is a ConfigParser file object
       - "metadata" a dictionary with default parameters
       - "section" is the section to look for in the cfg file
    This function looks for the same type of the value in the metadata dict.
    If it is a float or a list, it will look for a value or evaluate an expression.
    '''
    for option in metadata.keys():
        if type(metadata[option]) == str:
            try: 
                metadata[option] = cfg.get(section,option)    
            except cp.NoSectionError :
                if verb:
                    print 'getOptions:'+ section +', Using default value for '+option
            except cp.NoOptionError:
                if verb:
                    print 'getOptions:'+ section +', Using default value for '+option
        elif type(metadata[option]) == int or  type(metadata[option]) == float or type(metadata[option]) == list or type(metadata[option])==type(np.array([1])):
            # Quick fix to use also np.array type. Change later.
            try: 
                metadata[option] = cfg.getfloat(section,option)    
            except cp.NoSectionError :
                if verb:
                    print 'getOptions:'+ section +', Using default value for '+option
            except cp.NoOptionError:
                if verb:
                    print 'getOptions:'+ section +', Using default value for '+option
            except ValueError:
                metadata[option] = eval(cfg.get(section,option))
                if verb:
                    print 'getOptions:'+ section +', Evaluating '+cfg.get(section,option)+' for '+option
        elif type(metadata[option]) == bool:
            try: 
                metadata[option] = cfg.getboolean(section,option)    
            except cp.NoSectionError :
                if verb:
                    print 'getOptions:'+ section +', Using default value for '+option
            except cp.NoOptionError:
                if verb:
                    print 'getOptions:'+ section +', Using default value for '+option
    return metadata

##########################################################################
############ Functions to interact with NEURON ###########################
def insertRecorders(segment, labels, rec=None):
    ''' 
    Inserts recorders for NEURON state variables. 
    Note: labels is a dictionary. Example {'v': '_ref_v'}.
          Specify 'rec' to append to previous recorders.
          Records also time if 'rec' is 'None'(default).
    (Adapted from Daniele Linaro)
    '''
    if rec is None:
        rec = {'t': h.Vector()}
        rec['t'].record(h._ref_t)
    for k,v in labels.items():
        rec[k] = h.Vector()
        rec[k].record(getattr(segment, v))
    return rec

def passValuesToNeuron(parameters,safelist=[None],verb=False):
    '''
    Passes parameters from a dictionary to NEURON.
    If the element is in the safelist, 
    it assumes that the vectors have been created as objref and new Vector() 
    '''
    for ii in parameters.keys():
        if type(parameters[ii]) != type(np.array([])):
            if verb:
                print(ii + " = "+str(parameters[ii]))
            h("{"+ii + " = "+str(parameters[ii])+"}")
        else:
            if ii in safelist:
                exec("h."+ii+".from_python("+str(parameters[ii])+")")
            if verb:
                print("h."+ii+".from_python("+str(parameters[ii])+")")
            else:
                if verb:
                    print "passValuesToNeuron: "+ii+" not in Vectors safelist!"
##########################################################################

def runModel(tstop,recorders = None,parameters = None,vectorparameters = None,verbose = False):
    '''
    Clears the recorders.
    Passes parameters to NEURON.
    Initializes the model and runs until tstop.
    '''
#   print "Resetting model."
    for ii in recorders.itervalues():
       	ii.clear()
    h.rec_electrode.clear()
    passValuesToNeuron(parameters,verbose)
    passValuesToNeuron(vectorparameters,vectorparameters.keys(),verbose)
    h('resetModel()')
    h.initialize()
    h.init()
    h.continuerun(tstop)

def wasBlocked(spks,window=[40,50]):
    '''
    Checks if there was a there is a value of the np.array "spks" in the between the values of the list "window"
    Returns a boolean. True for block and False if there was a spike in that window. 
    '''
    try:
       	tmp=spks[(spks>window[0]) & (spks<window[1])]
    except IndexError:
       	tmp=np.array([])
    if len(tmp)<1:
       	return True
    return False

def calculateStimStep(block,block_1,stimStep, resolution = 0.9e-3, verbose = False):
    '''
    "block" and "block_1" are the current trial and the previous trial results, respectively.
    Implemetation of the binomial search.
        - If the last two trials are the different cut step in half.
        - If there was a block and the step is smaller than the resolution, return nan to stop algorithm.
        - Evaluate the signal of the step depending on the current trial.
    Returns the step to be added to the previous trial stimulation amplitude.
    '''
    stimStep = np.abs(stimStep)
    if not (block==block_1):
        stimStep = stimStep/2.0
        if verbose:
            print "updateStimAmplitude: Last two trials were different.Cutting in half."
    if stimStep<resolution and block:
        if verbose:
            print "updateStimAmplitude: There was block. Last step."
        return np.nan
    if block:
        return -stimStep
    else:
        return stimStep

def saveDictToH5Attrs(ID,labels, verbose = False):
    '''
    Saves a dictionary to an H5 group attributes.
    '''
    for kk,oo in labels.items():
        if type(oo)==bool:
            ID.attrs[kk] = int(oo)
       	else:
            ID.attrs[kk] = oo
        if verbose:
            print("Saved %s as attribute of %s."%(kk, ID.name))

def saveDictToH5ds(ID,labels,nchunks = 10000, verbose = False):
    '''
    Saves a dictionary to an H5 group.
    '''
    for kk,oo in labels.items():
        if len(oo)<1:
            print("Warning: [saveDictToH5ds] %s is empty. Skipping!" %(kk))
        else:
            if type(oo)==type(np.array([1])):
                if kk in ID.keys():
                    if verbose:
                        print("Warning: [saveDictToH5ds] %s already exists in %s. Appending!" %(kk, ID.name))
                    newshape=(ID[kk].shape[0]+np.size(oo),)
                    ID[kk].resize(newshape)
                    ID[kk][newshape[0]-np.size(oo):,]=oo[:]
                else:
                    ds = ID.create_dataset(kk,shape=(np.size(oo),),dtype='f',compression='gzip',maxshape=(None,))
                    ds[:]=oo[:]
                    if verbose:
                        print("Creating dataset %s in %s."%(kk,ID.name))

def parseCfgFile(cfgName, par, recpar, vecpar, verbose = False):
    '''
    Parses the configuration file for this particular simulation.
    '''
    fid = open(cfgName)
    if verbose:
        print "Parsing "+cfgName+"."
    cfg = cp.ConfigParser()
    cfg.readfp(fid)
    fid.close()
    par = getOptions(cfg,par,'MRGnode',verb=verbose)
    recpar = getOptions(cfg,recpar,'Recording',verb=verbose)
    vecpar = getOptions(cfg,vecpar,'IntracellularPattern',verb=verbose)
    vecpar["patternLag"] = np.array([ii-0.005 for ii in vecpar["pattern"]]) # Fix the lag to be used by Vector.play()

def searchForBlockThreshold(h,keyToIterate,passpar,recorders,vectorpar, verbose = False, plot = False):
    '''
    Uses a binomial search algorithm to search for the threshold to block a fiber.
    '''

    print "Running "+keyToIterate+" - " + str(passpar[keyToIterate]) +"!"
    stimAmp             = 0.0           # starts at zero.
    stimStep		= 1
    wasBlock		= False
    I			= []
    ts			= []
    startFrequencyRun   = time.clock()
    while not np.isnan(stimAmp):
       	passpar['HFSamp']     = stimAmp     # update stimulation amplitude
       	startSingleRun        = time.clock()
       	runModel(passpar['tstop'], recorders, passpar, vectorpar, verbose)
        print("Trial took %.1f sec. Used %3.3f mA, %5.0f Hz" %( (time.clock()-startSingleRun), h.HFSamp, h.HFSfrequency))
       	I.append(h.HFSamp)
       	spks = np.array(recorders['apc'].to_python())
       	ts.append(spks)
       	tmpBlock = wasBlocked(spks,[30, passpar['tstop']])
#       print tmpBlock, wasBlock
       	stimStep = calculateStimStep(tmpBlock, wasBlock, stimStep, verbose)
#	print stimStep
       	stimAmp  = stimAmp + stimStep
       	wasBlock = tmpBlock
        if plot:
            fig = plt.figure()
            for ii in recorders.keys():#range(0,int(h.axonnodes)):
                if verbose:
                    print("Plotting %s."%(ii))
                    print recorders[ii].to_python()
                if ii in ['apc']:
                    plt.plot(spks,np.ones(spks.shape),'|r',lw=2,ms=30)
                elif not ii in ['t']:
                    plt.plot(recorders['t'],recorders[ii])
            plt.plot(recorders['t'],h.rec_electrode.to_python)
            plt.show()
    return startFrequencyRun

def initModel(h,par,vecpar,recpar, verbose):
    '''
    Initializes the model.
    Creates the axon and so on.
    '''
    passValuesToNeuron(par,verb=verbose)
    h('{load_file("MRGnodeHFS.hoc")}')
    passValuesToNeuron(vecpar,[a for a in vecpar.keys()],verb=verbose)
	#h('{load_file("nrngui.hoc")}')
	#h('{load_proc("nrnmainmenu")}')
    h('{buildModel()}')
    if verbose:
        print("Passed parameters and built model.")
    # insert recorders and record action potential timestamps!
    segmentsToRecordV = [];
    segmentsNames = [int(h.intrinsicNode), int(h.HFSreferenceNode), int(h.axonnodes-1)] #segments to record
    for ii in segmentsNames:
       	segmentsToRecordV.append(h.node[ii](0.5))
    rec = None
    for ii in range(0,len(segmentsToRecordV)):
       	rec = insertRecorders(segmentsToRecordV[ii],{'node'+str(segmentsNames[ii]):'_ref_v'},rec)
    h.node[int(h.axonnodes-1)].push()
    apc = h.APCount(h.node[int(h.axonnodes-1)](0.5))
    apc.thresh               = 0
    rec['apc']               = h.Vector()
    apc.record(rec['apc'])
    #rec['electrodeWaveform'] = h.rec_electrode
    if verbose:
        print("Inserted recorders and APCount. Returning recorders.")
    ##########################################################################
    return rec

def main():
    verbose = True
    plot    = True
    if plot:
        import pylab as plt
    ###############PARSE SCRIPT OPTIONS#######################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "hc:d:", ["help", "cfg","data"])
    except getopt.GetoptError, err:
        sys.exit("Help: Usage: \n        -cfg=cfgName \n        -data=dataFolder.\n")
    dataFolder = None
    cfgName    = None
    par		   = g_par.copy()
    vecpar	   = g_vecpar.copy()
    recpar	   = g_recpar.copy()
    for o, a in opts:
       	if o in ("-h","--help"):
            sys.exit("Help: Usage: \n        -cfg=cfgName \n        -data=dataFolder\n        Or nothing for model.cfg and CWD.\n")
        elif o in ("-c", "--cfg"):
            cfgName = a
            if verbose:
                print("Reading %s."%cfgName)
            #print par
            parseCfgFile(cfgName, par, recpar, vecpar)
        elif o in ("-d","--data"):
            dataFolder = a
        else:
            assert False, "unhandled option"
	##########################################################################
    rec = initModel(h, par, vecpar, recpar, verbose)
    ###############LOADS PARAMETERS AND BUILDS THE MODEL######################
    #################SINGLE TRIAL OR ITERATE##################################
    parkeys=np.array([a for a in par.keys()])
    keyToIterate = parkeys[np.nonzero([type(a)==type(np.array([1.0])) for a in par.itervalues()])]
    if not len(keyToIterate):
        print "No parameter to iterate."
        keyToIterate = 'HFSfrequency'
        par['HFSfrequency']=[par['HFSfrequency']]
    else:
        keyToIterate=keyToIterate[0]
        print "Iterating over "+keyToIterate
        print keyToIterate + "list:"
        print par[keyToIterate]
    for ii in par[keyToIterate]:
       	passpar = par.copy()
       	passpar[keyToIterate] = ii
       	startFrequencyRun=searchForBlockThreshold(h, keyToIterate, passpar, rec, vecpar, verbose, plot)
       	if recpar['record']:
            simulation = {}
            for kk in rec.keys():
                simulation[kk] = np.array(rec[kk].to_python())
            simulation['electrodeWaveform'] = np.array(h.rec_electrode.to_python())
            fname = recpar['filename']
            if not path.isabs(fname):
                fname = path.abspath(fname)
            fid = h5.File(fname)
            passpar['ComputingTime'] = time.clock()-startFrequencyRun
            groupname = 'simulation'+str(len(fid.keys()))
            group = fid.create_group(groupname)
            saveDictToH5ds(group, simulation, verbose)
            saveDictToH5Attrs(group, passpar, verbose)
            saveDictToH5Attrs(group, {'pattern':vecpar['pattern']})
            fid.close()
    h.quit()


if __name__ == "__main__":
    main()
