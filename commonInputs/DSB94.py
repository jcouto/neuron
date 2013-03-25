1#!/usr/bin/env python
from CommonInput import *
import os
template = 'DSB94tmplt'

if os.path.isfile(template):
    h.load_file(template)
else:
    print('Could not load template [%s]'%template)
#    h.quit()

class DSB94(Neuron):
    """
    de Schutter E, Bower JM.(1994)
    An active membrane model of the cerebellar Purkinje cell. 
    I. Simulation of current clamps in slice. 
    The Journal of neuroscience, (71), 375-400.
    """
    class NeuronProperties (tbl.IsDescription):
        name = tbl.StringCol(32)
        description = tbl.StringCol(128)
        nSynapses = tbl.Int32Col()

    class SynapseProperties (tbl.IsDescription):
        name = tbl.StringCol(32)
        description = tbl.StringCol(128)
        tauRise = tbl.Float64Col()
        tauDecay = tbl.Float64Col()
        Erev = tbl.Float64Col()
        gbar = tbl.Float64Col()

    class PhaseResponseCurveProperties (tbl.IsDescription):
        frequency = tbl.Float64Col()
        gp = tbl.Float64Col()
        gi = tbl.Float64Col()
        pulseAmp = tbl.Float64Col()
        pulseWidth = tbl.Float64Col()
        spkCount = tbl.Float64Col()

    def __init__(self, ID,neuronProps={'nSynapses': 1474},
                 synapseProps={'name': 'biexp', 'Erev': 0., 'tauRise':0.5, 'tauDecay':1.2},
                 prcProps={'frequency':None, 'gp':0.01, 'gi':0.1, 'pulseAmp':0.05,'pulseWidth':1,'spkCount':6},
                 verbose=False):
        """
        Sets the parameters of the neuron and the calls the constructor of the parent class.
        Parameters:
        ID - the unique identifier of the neuron
        nsyn - number of synapses
        neuronProps - a dictionary with the properties of the neuron, namely its length, diameter and
        the number of synapses
        synapseProps - a dictionary with the properties of the synapses attached to the neuron, namely
        the type and the reversal potential (in case of an alpha synapse, also the rise and fall times
        are required)
        """

        self._neuronProps = neuronProps
        self._neuronProps['name'] = 'DSB94'
        self._neuronProps['description'] = 'De Schutter and Bower Purkinje neuron model'
        
        self._synapseProps = synapseProps
        if self._synapseProps['name'].lower() == 'ampa':
            self._synapseProps['description'] = 'AMPA synapse'
            self._synapseProps['tauRise'] = -1.   # not used
            self._synapseProps['tauDecay'] = -1.  # not used
        elif self._synapseProps['name'].lower() == 'biexp':
            self._synapseProps['description'] = 'Bi-exponential synapse'
            self._synapseProps['tauRise'] = 0.5
            self._synapseProps['tauDecay'] = 1.2
            self._synapseProps['gbar'] = 0.7
        else:
            del self._synapseProps
            raise Exception('Unknown synaptic model [%s]' % synapseProps['name'])
        
        self._prcProps = prcProps 
        
        Neuron.__init__(self, ID, verbose)
       

    @property
    def dendrite(self):
        return self._dendrite
    @property
    def synapseProps(self):
        return self._synapseProps
    @property

    def nsyn(self):
        """
        TODO: replace with nSynapses
        Returns the number of synapses, this is the number of synapses to be added, not the total number of synapses if _addSynapses(location=['dendrite'|'soma']) is called. 
        """
        return self._neuronProps['nSynapses']

    def synapseCluster(self, size):
        for i in range(self.nsyn/size):
            yield i*size + np.arange(size)

    def _build(self):
        if self.verbose:
            print('>>> Building a De Schutter and Bower (1994) neuron model.')
        self._cell = h.DSB94(str(self._ID),'DSB94','neuroConstruct version')
        self._soma = self._cell.soma
        self._dendrite = self._cell.dendrite_group
        h.celsius = 37
        self._dendriticSynapseSections = [] 
        self._resetDendriticSynapseSections()        

    def _resetDendriticSynapseSections(self):
        """
        Sets the equal probability to draw a synapse from each location
        TODO: Use a probability (0-1) instead of a Boolean.
        """
        if self.verbose:
            print("Reseting _unusedSectionIndexes.")
        self._dendriticSynapseSections = [False for tmp in self._dendrite]

    def _draw(self,lst):
        """
        Draws an element randomly from a list.
        """
        randomIndex = int(np.random.uniform(len(lst)))
        if self.verbose:
            print("Selected element [%d]."%randomIndex)
        return lst[randomIndex]

    def _unusedSectionIndexes(self):
        """
        Returns the indexes of unused (dendritic) sections.
        """    
        tmp = [i for i, e in enumerate(self._dendriticSynapseSections) if e == 0] 
        self._resetDendriticSynapseSections()                    
        return tmp
        
    def _drawDendrite(self):
        """
        Draws a synaptic location that has not yet been used in the dendrite.
        In case the maximal number of synapses has been reached, it restarts the draw.
        """
        idx = self._draw(self._unusedSectionIndexes())
        self._dendriticSynapseSections[idx]=True
        return self._dendrite[idx]

    def _addSynapses(self,location='dendrite'):
        """
        Adds all synapses to the neuron. Distribute the synapses randomly thoughout the dendrite.
        """
        if self.verbose:
            print('>>> Adding synapses to the model.')
        # TODO: take synapse parameters from:
        # Jaeger, D., De Schutter, E., & Bower, J. M. (1997).
        # The role of synaptic and voltage-gated currents in the control of Purkinje cell
        # spiking: a modeling study.
        # The Journal of Neuroscience, 17(1), 91-106.
        for i in range(self.nsyn):
            if location=='dendrite':
                if self.verbose:
                    print("Adding dendritic synapse.")
                sec = self._drawDendrite()
            else:
                if self.verbose:
                    print("Adding somatic synapse.")
                sec = self.soma
            if self._synapseProps['name'].lower() == 'ampa':
                syn = h.AMPA_S(sec(0.5))
                h.Erev_AMPA_S = self._synapseProps['Erev']
            elif self._synapseProps['name'].lower() == 'biexp':
                syn = h.Exp2Syn(sec(0.5))
                syn.tau1 = self._synapseProps['tauRise']
                syn.tau2 = self._synapseProps['tauDecay']
                syn.e = self._synapseProps['Erev']
            else:
                print('Unknown synapse type [%s].' % synapseProps['type'])
                return
            self.synapses.append(syn)

    def _setNumberSynapses(self, nsyn):
        """
        Sets the number of synapses.
        """
        self._neuronProps['nSynapses'] = nsyn

    def addDendriticSynapses(self, nsyn):
        """
        Adds nsyn synapses to the dendrite and random places.
        """
        self._setNumberSynapses(nsyn)
        self._addSynapses(location='dendrite')
 
    def addSomaticSynapses(self, nsyn):
        """
        Adds nsyn synapses to the soma.
        """
        self._setNumberSynapses(nsyn)
        self._addSynapses(location='soma')
    
    def setPRCfrequency(self,frequency):
        """
        Changes the PRC estimator frequency property.
        """
        self._prcProps['frequency']=frequency
        if not self._prc_sobol is None:
            self._prc_sobol.F = frequency
        return

    def _addPRCestimator(self):
        """
        Adds a PRC estimator from intrinsic properties.
        """
        prcProps = {}
        prcProps['F'] = self._prcProps['frequency']
        prcProps['gp'] = self._prcProps['gp']
        prcProps['gi'] = self._prcProps['gi']
        prcProps['amp'] = self._prcProps['pulseAmp']
        prcProps['pw'] = self._prcProps['pulseWidth']
        prcProps['spkCount'] = self._prcProps['spkCount']
        self.addPRCestimator(prcProps,-20)


def main():
    import pylab as plt
    N = 1
    nsyn = 0
    fixedInput = {'probability': 0.5, 'spikeTimes': np.cumsum(np.random.poisson(30,140))}
    neurons = [DSB94(i,{'nSynapses':nsyn}) for i in range(1)]
    #for n in neurons:
    #    n.addFixedInput(fixedInput['probability'], fixedInput['spikeTimes'])
    #    n.addPoissonInputs(clusterSize=5, stimulusProps={'frequency': 0.1, 'noise': 1.})
    #    n.addSomaticVoltageRecorder()

    time = h.Vector()
    time.record(h._ref_t)

    h.load_file('stdrun.hoc')
    h.tstop = 100
    h.run()
    filename = 'DSBneuron.h5'
    for n in neurons:
        print('Neuron [%d] emitted %d spikes.' % (n.ID, len(n.spikeTimes())))
        #os.remove(filename)
        saveNeuron(filename, n)
        plt.plot(time,n.somaticVoltage())
    plt.show()
    h.quit()

if __name__ == '__main__':
    main()
