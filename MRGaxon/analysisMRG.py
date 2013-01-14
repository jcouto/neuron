#!/usr/bin/env python
import os.path as path
import numpy as np
import h5py as h5

__all__ = ['MRGfile','nonzero_len','fix_axes']

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
        self.fiberD = []
        self.axonnodes = []
        self.output = []
        self.counter = -1
        if not filename is None:
            self.append(filename)

    def append(self,filename,verbose=False):
        ''' Appends a file and the fibers in that file to the class.
        '''
        if verbose:
            print('Reading data from file ' + filename+'.')
        self._file_ref.append(h5.File(filename))
        self.counter+=1
        for v,k in self._file_ref[-1].iteritems():
            self.HFSfrequency.append(k.attrs['HFSfrequency'])
            self.HFSamp.append(k.attrs['HFSamp'])
            self.HFSwaveform.append(k.attrs['HFSwaveform'])
            self.HFSx.append(k.attrs['HFSx'])
            self.HFSy.append(k.attrs['HFSy'])
            self.HFSz.append(k.attrs['HFSz'])
            self.fiberD.append(k.attrs['fiberD'])
            self.axonnodes.append(k.attrs['axonnodes'])
            self.names = k
            self.file.append(self.counter)
            output_name = 'spk'+str(np.max(k.attrs['nodes']))
            if len(k['spiketimes/'+output_name]):
                self.output.append(np.array(k['spiketimes'+'/'+output_name].value))
            else:
                self.output.append(np.array([]))
                
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

def fix_axes(ax=None,xlabel='',ylabel='',xloc='bottom',yloc='left',xposition=('outward',10),yposition=('outward',10),lw=1,fontsize=10,xcolor='black',ycolor='black'):
    notxloc='bottom'
    notyloc='right'
    if xloc=='bottom':
        notxloc='top'
    if yloc=='right':
        notyloc='left'
    if ax==None:
        ax=plt.gca()
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
    for label in ax.get_xticklabels():
            label.set_color(xcolor)
            label.set_fontsize(fontsize)
    for label in ax.get_yticklabels():
            label.set_color(ycolor)
            label.set_fontsize(fontsize)
    for label in ax.get_xticklines():
            label.set_color(xcolor)
            label.set_lw(lw)
            label.set_markersize(6)
    for label in ax.get_yticklines():
            label.set_color(ycolor)
            label.set_lw(lw)
            label.set_markersize(6)
