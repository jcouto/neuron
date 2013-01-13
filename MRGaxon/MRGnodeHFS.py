#!/usr/bin/env python
from neuron import h
import os.path as path
import numpy as np
import h5py as h5


__all__ = ['insert_nrn_recorders','pass_parameters_to_nrn','createMRGaxon',
           'recordMRGaxon','fix_patternLag_vector','record_node_voltage',
           'record_node_spikes', 'plotMRGaxon', 'readConfigurations',
           'resetRecorder','runMRGaxon','updateMRGaxon','append_fiber_to_file']

# Define global/default parameters that are used when probing the configuration file
g_par = { 
    'dt': 0.002,
    'tstop': 50,
    'channelDescription': 0,
    'axonnodes': 51,
    'fiberD': 10,
    'HFSreferenceNode': 25,
    'HFSdur': 50.0,
    'HFSfrequency': 200,
    'HFSpolarity': 1.0,
    'HFSdelay': 0,
    'HFSpulsewidth':0.09,
    'HFSamp': 1.0,
    'HFSwaveform': 0,
    'HFSx': 0.0,
    'HFSy': 0.0,
    'HFSz': 1000.0,
    'intrinsicNode': 0,
    'intrinsicDur': 0.1,
    'intrinsicAmp': 2.0,
    'pattern': np.array([40.0]),
    'patternLag': np.array([39.9])
    }
g_recpar = {
    'record': True,
    'nodes':np.array(range(0,g_par['axonnodes'])),
    'filename': 'data/simulation.h5'
    }

def insert_nrn_recorders(segment, labels, rec=None):
    ''' 
    Inserts recorders for NEURON state variables.
    Use one per segment.
    "labels" is a dictionary. 
    Example {'v': '_ref_v'}.
    Specify 'rec' to append to previous recorders.
    Records also time if 'rec' is 'None' (default).
    (Aknowledgements: Daniele Linaro)
    '''
    if rec is None:
        rec = {'t': h.Vector()}
        rec['t'].record(h._ref_t)
    for k,v in labels.items():
        rec[k] = h.Vector()
        rec[k].record(getattr(segment, v))
    return rec

def pass_parameters_to_nrn(parameters, exception = [], verb=False):
    '''
    Passes parameters from a dictionary to NEURON.
    If the the element is a vector it assumes that the a vector 
    has been created as objref and new Vector() in the hoc code.
    Items in 'exception' list are not submitted. 
    Set 'verb' to True for verbose.
    '''
    for k,v in parameters.iteritems():
        if k not in exception:
            if type(v) is not type(np.array([])):
                h("{"+ k + " = "+str(v)+"}")
                if verb:
                    print(k + " = "+str(v))
            else:
               # exec("h."+ k +".from_python("+str(v)+")")
                getattr(h,k).from_python(v)
                if verb:
                    print("h."+k+".from_python("+str(v)+")")

def fix_patternLag_vector(parameters):
    '''
    Changes the patternLag variable according to the pattern vector.
    This is used by the Play method of the vector. 
    Ideally this should be done in the NEURON script.
    '''
    if 'pattern' in parameters.keys():
        parameters["patternLag"] = np.array([ii-0.005 for ii in parameters["pattern"]])
    else:
        print('Key "pattern" not found in this dictionary. Did nothing.')

def record_node_spikes(nodenumber, rec=None, apc=None, threshold = -15):
    '''
    Records the action potentials of a particular set of nodes.
    Returns a "rec" dictionary.
    '''
    if rec is None:
        rec = {}
    if apc is None:
        apc = {}
    for n in nodenumber:
        apc['apc'+str(n)] = h.APCount(h.node[int(n)](0.5))
        apc['apc'+str(n)].thresh = threshold
        rec['spk'+str(n)] = h.Vector()
        apc['apc'+str(n)].record(rec['spk'+str(n)])
    return rec,apc

def record_node_voltage(nodenumber, rec=None):
    '''
    Records the membrane potential of a particular set of nodes.
    '''
    rec = None
    segments = []
    for n in nodenumber:
        segments.append(h.node[n](0.5))
    for seg,n in zip(segments,nodenumber):
       	rec = insert_nrn_recorders(seg,{'v_node'+str(n):'_ref_v'},rec)
    return rec

def createMRGaxon(par, verbose):
    '''
    Initializes the model.
    Creates the axon and stimulation according to the parameters.
    '''
    h('{load_file("stdgui.hoc")}')
    fix_patternLag_vector(par)
    pass_parameters_to_nrn(par,['pattern','patternLag'],verb=verbose)
    h('{load_file("MRGnodeHFS.hoc")}')
    pass_parameters_to_nrn(par, verb=verbose)    
    h('{buildModel()}')

def updateMRGaxon(par,verbose):
    '''
    Updates the parameters of the model.
    '''
    fix_patternLag_vector(par)
    pass_parameters_to_nrn(par,verb=verbose)
    h.resetModel()

def recordMRGaxon(recpar,verbose):
    '''
    Inserts the recorders as specified in recpar.
    '''
    k = recpar['nodes']
    rec = {}
    rec['voltage'] = record_node_voltage(k)
    rec['spiketimes'],rec['apcount'] = record_node_spikes(k)
    if verbose:
        print('Now recording from '+str(k))
    return rec

def resetRecorder(rec,verbose=False):
    '''
    Clears hoc vectors in spiketimes and voltage and resets apcounts.
    '''
    for k,o in rec['spiketimes'].iteritems():        
        if verbose:
            print('Reseting ' + k)
        o.clear()
    for k,o in rec['voltage'].iteritems():
        if verbose:
            print('Reseting ' + k)
        o.clear()
    for k,o in rec['apcount'].iteritems():
        if verbose:
            print('Setting apcount ' + k + ' to zero.' )
        o.n = 0
    

def plotRastergram(plt,spiketrains,offset,color=[0.1,0.1,0.1]):
    '''
    Plots a rastergram of the spike trains
    Spiketrains should be a list of spiketrains
    '''
    for i,sp in enumerate(spiketrains):
        X = np.tile(sp,(2,1))
        Y = np.transpose(np.tile(np.transpose([i,i+1]),(len(sp),1))+offset)
        plt.plot(X,Y,color=color,lw=1)
    
def plotMRGaxon(plt, rec, recpar,color=[0,0,0]):
    '''
    Plots the voltage traces and a rastergram of the spikes counting ordered by node.
    '''
    fig = plt.figure(figsize=(10,5))
    ax = []
    ax.append(fig.add_axes([0.1,0.1,0.8,0.2]))
    spiketimes = []
    n_sptrain = len(spiketimes)
    
    for offset,ii in enumerate(recpar['nodes']):
        spiketimes.append(rec['spiketimes']['spk'+str(ii)].to_python())
    plotRastergram(ax[-1], spiketimes, 0, color)

    ax.append(fig.add_axes([0.1,0.4,0.8,0.6]))
    voltages = rec['voltage']
    time = voltages['t']
    n_voltage = len(voltages) - 1
    for ii in recpar['nodes']:
        ax[-1].plot(time,voltages['v_node'+str(ii)])
    return fig

def process_configuration(cp, cfg, metadata, section):
    '''
    Parse configuration file from default values.
       - "cfg" is a ConfigParser file object
       - "metadata" a dictionary with default parameters
       - "section" is the section to look for in the cfg file
    This function looks for the same type of the value in the metadata dict.
    If it is a float or a list, it will look for a value or evaluate an expression.
    '''
    output = metadata.copy()
    if not cfg.has_section(section):
        return output
    for option in metadata.keys():
        if cfg.has_option(section,option):
            if type(metadata[option]) is str:
                output[option] = cfg.get(section,option)
            elif type(metadata[option]) is int or  \
                type(metadata[option]) is float or \
                type(metadata[option]) is list or \
                type(metadata[option])==type(np.array([1])):
                try: 
                    output[option] = cfg.getfloat(section,option)    
                except ValueError:
                    output[option] = eval(cfg.get(section,option))
            elif type(metadata[option]) == bool:
                output[option] = cfg.getboolean(section,option)    
    return output

def append_fiber_to_file(rec,par,recpar,group=None):
    '''
    Uses h5py to append a fiber to a file.
    '''
    fid = h5.File(recpar['filename'],'a')
    n_fiber = len(fid.keys())
    if group is None:
        gid = fid.create_group('fiber'+str(n_fiber))
    else:
        gid = fid.create_group(group)

def runMRGaxon():
    h.resetModel()
    h.run()

def readConfigurations(filename):
    '''
    Reads the parameters from the configuration file.
    Returns par and recpar or defaults if filename is not a valid file.
    '''
    if path.isfile(filename):
        fid = open(filename)
    else:
        print('File "'+filename+'" does not exist.')
        return g_par,g_recpar
    import ConfigParser as cp
    cfg = cp.ConfigParser()
    cfg.readfp(fid)
    fid.close()
    par = process_configuration(cp, cfg, g_par, 'MRGnode')
    recpar = process_configuration(cp, cfg, g_recpar, 'Recording')
    return par, recpar

def main():
    '''
    Executes a simulation taking the parameters from a configuration file.
    '''
    verbose = True
    plot = True
    counter = 0
    for v in sys.argv:
        counter+=1
        if path.basename(__file__) in v:
            break
    filename = sys.argv[counter]
    par, recpar = readConfigurations(filename)
    createMRGaxon(par,0)
    rec = recordMRGaxon(recpar,0)
    runMRGaxon()
    if recpar['record']:
        append_fiber_to_file(rec,par,recpar)
    if plot:
        import pylab as plt
        fig = plotMRGaxon(plt, rec, recpar)
        plt.show()
if __name__ == '__main__':
    main()
