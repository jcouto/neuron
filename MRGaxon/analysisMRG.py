#!/usr/bin/env python
import os
import os.path as path
import numpy as np
import h5py as h5
import pylab as plt

__all__ = ['MRGfile','zero_len','nonzero_len','fix_axes']

class MRGfile():
    def __init__(self, filename=None,verbose=False):
        """ Initializes the file object and reads the fiber attributes from the file.
        """
        
        self._file_ref = []
        if verbose:
            print('Reading data from file ' + bind+'.')
        
        self.names = []
        self.file = []
        self.HFSfrequency = []
        self.HFSamp = []
        self.HFSwaveform = []
        self.HFSx = []
        self.HFSy = []
        self.HFSz = []
        self.nChannels = []
        self.fiberD = []
        self.axonnodes = []
        self.recordedNodes = []
        self.output = []
        self.counter = -1
        if not filename is None:
            self.append(filename)

    def append(self,filename,verbose=False):
        ''' Appends a file and the fibers in that file to the class.
        '''
        if verbose:
            print('Reading data from file ' + filename+'.')
        if path.isfile(filename):
            try: 
                self._file_ref.append(h5.File(filename, 'r'))
            except IOError:
                    if not os.access(filename, os.R_OK):
                        print('File %s not readable'%(filename))
                    else: 
                        print('Could not append file (IOError): %s '%(filename))
            else:
                self.counter+=1
                for v,k in self._file_ref[-1].iteritems():
                    self.HFSfrequency.append(k.attrs['HFSfrequency'])
                    self.HFSamp.append(k.attrs['HFSamp'])
                    self.HFSwaveform.append(k.attrs['HFSwaveform'])
                    self.HFSx.append(k.attrs['HFSx'])
                    self.HFSy.append(k.attrs['HFSy'])
                    self.HFSz.append(k.attrs['HFSz'])
                    try: 
                        self.nChannels.append(k.attrs['n_na'])
                    except: # This was not always like this...  
                        self.nChannels.append(1.0)
                    self.fiberD.append(k.attrs['fiberD'])
                    self.axonnodes.append(k.attrs['axonnodes'])
                    self.recordedNodes.append(k.attrs['nodes'])
                    self.names.append(k)
                    self.file.append(self.counter)
                    output_name = 'spk'+str(np.max(k.attrs['nodes']))
                    if len(k['spiketimes/'+output_name]):
                        self.output.append(np.array(k['spiketimes'+'/'+output_name].value))
                    else:
                        self.output.append(np.array([]))
        else:
            print('File %s not found.'%(filename)) 
    
    def spiketrains(self,index=None):
        '''Returns the output spiketimes for the specified indexes
        '''
        tmp = []
        for idx in index:
            tmp.append(self.output[idx])
        return tmp
    def voltage(self,index=None, node=None):
        '''Returns the voltage recorded at a specific node for the given indexes.
        '''
        tmp = []
        for idx in index:
            if node is None:
                # Use the output node.
                output_name = 'v_node'+str(np.max(self.names[idx].attrs['nodes']))
            else:
                output_name = 'v_node'+str(node)
            v = self.names[idx]['voltage/'+output_name][()]
            tmp.append(v)
        return np.array(tmp).T

    def time(self,fiber = [0]):
        '''Returns the time vector at a given fiber.
        If fiber is not specified, it returns it for fiber 1.
        '''
        tmp = []
        for idx in fiber:
            output_name = 't'
            v = self.names[idx]['voltage/'+output_name][()]
            tmp.append(v)
        return np.array(tmp).T

    def close(self):
        for fid in self._file_ref:
            fid.close()


def nonzero_len(A,B=None):
    '''Given that A and B are 2 lists, min_nonzero_len(A,B) returns the subset of A that 
    '''
    if B is None:
        B=A
    A=np.array(A)
    return A[np.nonzero([len(i) for i in B])[0]]

def zero_len(A,B=None):
    '''Given that A and B are 2 lists, min_nonzero_len(A,B) returns the subset of A that 
    '''
    if B is None:
        B=A
    A=np.array(A)
    print np.nonzero([len(i)==0 for i in B])
    return A[np.nonzero([len(i)==0 for i in B])[0]]

def fix_axes(ax=None,xlabel='',ylabel='',xloc='bottom',yloc='left',
             xposition=('outward',0),yposition=('outward',0),
             lw=1,fontsize=10,xcolor='black',ycolor='black'):
    notxloc='bottom'
    notyloc='right'
    if xloc=='bottom':
        notxloc='top'
    if yloc=='right':
        notyloc='left'
    if ax==None:
        ax=plt.gca()
    ax.tick_params(axis='both', direction='out')
    ax.spines[notxloc].set_visible(False)
    ax.spines[notyloc].set_visible(False)
    ax.spines[xloc].set_visible(True)
    ax.spines[yloc].set_visible(True)
    ax.spines[xloc].set_lw(lw)
    ax.spines[xloc].set_color(xcolor)
    ax.spines[xloc].set_position(xposition)
    ax.spines[yloc].set_lw(lw)
    ax.spines[yloc].set_color(ycolor)
    ax.spines[yloc].set_position(yposition)
    ax.yaxis.set_ticks_position(yloc)
    ax.xaxis.set_ticks_position(xloc)
    xax = ax.xaxis
    xax.set_label_text(xlabel)
    xax.set_label_position(xloc)
    xax.set_label_position(xloc)
    
    xax.label.set_color(xcolor)
    yax = ax.yaxis
    yax.set_label_text(ylabel)
    yax.set_label_position(yloc)
    yax.set_label_position(yloc)
    yax.label.set_color(ycolor)
    for label in [xax.label,yax.label]:
        label.set_size(fontsize)
    for label in ax.get_xticklabels():
            label.set_color(xcolor)
            label.set_fontsize(fontsize)
    for label in ax.get_yticklabels():
            label.set_color(ycolor)
            label.set_fontsize(fontsize)
    for label in ax.get_xticklines():
            label.set_color(xcolor)
            label.set_lw(lw)
            label.set_markersize(4)
    for label in ax.get_yticklines():
            label.set_color(ycolor)
            label.set_lw(lw)
            label.set_markersize(4)
