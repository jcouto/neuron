#!/usr/bin/env python

from neuron import h
import numpy as np
import h5utils as h5
import tables as tbl

class Neuron:
    """
    Base class for implementing model neurons that have two kinds of synaptic inputs:
    1) A fixed input delivered to a subset of the synapses. The times at which the
       synapses will be activated by the fixed input is hard-coded.
    2) Many Poisson spike trains delivered to small clusters of synapses. Each Poisson
       train is independent.
    """

    class NeuronProperties (tbl.IsDescription):
        pass

    class SynapseProperties (tbl.IsDescription):
        pass

    class PhaseResponseCurveProperties (tbl.IsDescription):
        pass

    def __init__(self, ID, verbose=False):
        """
        Constructs the morphology of the neuron, adds synapses to it and finally adds
        a spike counter that keeps track of the times at which the neuron fired.
        """
        self._ID = ID
        self._verbose = verbose
        self._soma = None
        self._netcons = {}
        self._synapses = []
        self._prc_netcon = None
        self._prc_sobol = None
        self._prc_sobol_times=None
        self._build()
        self._addSynapses()
        self._addSpikeCounter()
        self._neuronPropsClass = self.NeuronProperties
        self._synapsePropsClass = self.SynapseProperties
        self._prcPropsClass = self.PhaseResponseCurveProperties
    @property
    def ID(self):
        """
        Returns the unique identifier of this neurons
        """
        return self._ID

    @property
    def verbose(self):
        """
        Returns whether the neuron should print informational messages.
        """
        return self._verbose

    @property
    def soma(self):
        """
        Returns a section corresponding to the soma of the model.
        """
        return self._soma

    @property
    def synapses(self):
        """
        Returns a list that contains all the synapses of the model.
        """
        return self._synapses

    @property
    def connections(self):
        """
        Returns a dictionary whose keys are the id's of the synapses
        and whose values are lists containing all the NetCon objects
        connected to the corresponding synapse.
        """
        return self._netcons

    @property
    def numSynapses(self):
        """
        Returns the number of synapses in the model neuron.
        """
        return len(self.synapses)

    @property
    def neuronProperties(self):
        return self._neuronProps

    @property
    def neuronPropertiesClass(self):
        return self._neuronPropsClass

    @property
    def synapseProperties(self):
        return self._synapseProps

    @property
    def synapsePropertiesClass(self):
        return self._synapsePropsClass
    @property
    def phaseResponseCurveProperties(self):
        return self._prcProps

    @property
    def phaseResponseCurvePropertiesClass(self):
        return self._prcPropsClass

    def spikeTimes(self):
        """
        Returns a Numpy array containing the times at which the cell spiked.
        """
        return np.array(self._spikeTimes)

    def somaticVoltage(self):
        """
        Returns a Numpy array containing the trace of the membrane voltage.
        """
        return np.array(self._somaticVoltage)

    def fixedInputSpikeTimes(self):
        """
        Returns the (fixed) times at which presynaptic spikes where delivered to this neuron.
        """
        try:
            spks = np.array(self._fixedInput['spikeTimes'])
        except:
            spks = np.array([])
        return spks

    def synapseCluster(self, size):
        """
        Should return a list of length size that contains the indices of neighboring
        synapses. This is used to connect the same Poisson stimulus to a ``logical''
        group of synapses.
        """
        raise NotImplementedError()

    def addSomaticVoltageRecorder(self):
        """
        Adds a recorder of the somatic voltage.
        """
        self.removeSomaticVoltageRecorder()
        self._somaticVoltage = h.Vector()
        self._somaticVoltage.record(self.soma(0.5)._ref_v)

    def removeSomaticVoltageRecorder(self):
        """
        Removes the recorder of the somatic voltage.
        """
        try:
            del self._somaticVoltage
        except:
            pass

    def connectStimulus(self, stimulus, synapseIndexes,
                        connectionProps={'weight': 0.0001, 'delay': 1.}):
        """
        Connects a stimulus to one or more synapses.
        Parameters:
        stimulus - an object that is able to generate spikes.
        synapseIndexes - a list containing the index of the synapses to which
        the stimulus should be connected.
        connectionProps - a dictionary containing the parameters of the connections,
        namely the weight and the delay.
        """
        for idx in synapseIndexes:
            self._connectSingleStimulus(stimulus, idx, connectionProps)

    def addFixedInput(self, probability, spikeTimes, connectionProps={'weight': 0.0001, 'delay': 1.}):
        """
        Adds a fixed input to a random subset of synapses.
        Parameters:
        probability - the probability that the fixed input will contact a synapse.
        spikeTimes - an array containing the instants at which presynaptic spikes are emitted.
        connectionProps - a dictionary containing the parameters of the connections,
        namely the weight and the delay.
        """
        self._fixedInput = {
            'spikeTimes': h.Vector(spikeTimes),
            'stimulus': h.VecStim()
            }
        self._fixedInput['stimulus'].play(self._fixedInput['spikeTimes'])
        idx = np.nonzero(np.random.uniform(size=self.numSynapses) < probability)[0]
        self.connectStimulus(self._fixedInput['stimulus'], idx, connectionProps)

    def addPoissonInputs(self, clusterSize=1, stimulusProps={'frequency': 1., 'noise': 1.},
                         connectionProps={'weight': 0.0001, 'delay': 1.}):
        """
        Adds a population of presynaptic Poisson generators to groups of synapses.
        Parameters:
        clusterSize - the size of the groups of synapses that will receive the same input.
        stimulusProps - a dictionary with the properties of the poissonian stimuli, namely the
        frequency (in Hertz) and the noise (see the documentation of NetStim for the
        meaning of this parameter.
        connectionProps - a dictionary containing the parameters of the connections,
        namely the weight and the delay.
        """
        self._stimuli = []
        for idx in self.synapseCluster(clusterSize):
            ns = h.NetStim()
            ns.interval = 1000./stimulusProps['frequency']
            ns.number = 1e9
            ns.start = 0
            ns.noise = stimulusProps['noise']
            self.connectStimulus(ns, idx, connectionProps)
            self._stimuli.append(ns)

    def _build(self):
        """
        Constructs the morphology of the neuron and add all mechanisms.
        """
        raise NotImplementedError()

    def _addSynapses(self):
        """
        Add all synapses to the model neuron.
        """
        raise NotImplementedError()

    def _removeSynapses(self):
        """
        Removes all synapses to the model neuron.
        """
        try:
            del self._synapses
            del self._netcons
            self._synapses = []
            self._netcons = {}
        except:
            pass

    def _addSpikeCounter(self, threshold=-20):
        """
        Adds a spike counter to the model neuron.
        """
        if self.verbose:
            print('>>> Adding a spike counter to the model.')
        self._spikeTimes = h.Vector()
        self._apc = h.APCount(self.soma(0.5))
        self._apc.thresh = threshold
        self._apc.record(self._spikeTimes)

    def _connectSingleStimulus(self, stimulus, synapseIndex, connectionProps):
        """
        Connects a stimulus to a given synapse.
        """
        syn = self.synapses[synapseIndex]
        nc = h.NetCon(stimulus, syn)
        nc.weight[0] = connectionProps['weight']
        nc.delay = connectionProps['delay']
        if not synapseIndex in self.connections:
            self.connections[synapseIndex] = []
        self.connections[synapseIndex].append(nc)
    def addPRCestimator(self, prcProps, threshold=-20):
        """
        Connects a PRCestimator to the neuron soma.
        """
        self._prc_sobol = h.SobolPulses(0.5,sec=self.soma)
        self._prc_netcon = h.NetCon(self.soma(0.5)._ref_v,self._prc_sobol,sec=self.soma)
        self._prc_netcon.delay = 0
        self._prc_netcon.threshold = threshold
        self._prc_sobol.gp = prcProps['gp']
        self._prc_sobol.gi = prcProps['gi']
        self._prc_sobol.F = prcProps['F']
        self._prc_sobol.amp = prcProps['amp']
        self._prc_sobol.pw = prcProps['pw']
        self._prc_sobol.delay=0
        self._prc_sobol.spkCount=prcProps['spkCount']
        self._prc_sobol_times = h.Vector()
        nc = h.NetCon(self._prc_sobol._ref_iPulse, None,sec=self._soma)
        nc.delay = 0
        nc.threshold = self._prc_sobol.amp/2
        nc.record(self._prc_sobol_times)
    def _addPRCestimator(self):
        """
        Adds a PRC estimator from intrinsic properties.
        """
        raise NotImplementedError()
        
    def perturbationTimes(self):
        """
        Returns the times of the perturbation.
        """
        if self._prc_sobol_times is None:
            return np.array([])
        else:
            return np.array(self._prc_sobol_times)

class KhaliqRaman(Neuron):
    """
    Single compartment model of a Purkinje cell based on the paper:
    Khaliq, Z., Gouwens, N., & Raman, I. (2003).
    The Contribution of Resurgent Sodium Current to High-Frequency Firing in
    Purkinje Neurons: An Experimental and Modeling Study.
    The Journal of neuroscience, (12), 4899-4912.
    """
    class NeuronProperties (tbl.IsDescription):
        name = tbl.StringCol(32)
        description = tbl.StringCol(128)
        length = tbl.Float64Col()
        diameter = tbl.Float64Col()
        nSynapses = tbl.Int32Col()

    class SynapseProperties (tbl.IsDescription):
        name = tbl.StringCol(32)
        description = tbl.StringCol(128)
        tauRise = tbl.Float64Col()
        tauDecay = tbl.Float64Col()
        Erev = tbl.Float64Col()
    class PhaseResponseCurveProperties (tbl.IsDescription):
        frequency = tbl.Float64Col()
        gp = tbl.Float64Col()
        gi = tbl.Float64Col()
        pulseAmp = tbl.Float64Col()
        pulseWidth = tbl.Float64Col()
        spkCount = tbl.Float64Col()
        delay = tbl.Float64Col()

    def __init__(self, ID, neuronProps={'length': 20, 'diameter': 20, 'nSynapses': 100},
                 synapseProps={'name': 'ampa', 'Erev': 0.},
                 prcProps={'frequency':None, 'gp':0.01, 'gi':0.1, 'pulseAmp':0.05,'pulseWidth':1,'spkCount':6},
                 verbose=False):
        """
        Sets the parameters of the neuron and the calls the constructor of the parent class.
        Parameters:
        ID - the unique identifier of the neuron
        neuronProps - a dictionary with the properties of the neuron, namely its length, diameter and
        the number of synapses
        synapseProps - a dictionary with the properties of the synapses attached to the neuron, namely
        the type and the reversal potential (in case of an alpha synapse, also the rise and fall times
        are required)
        """
        self._neuronProps = neuronProps
        self._neuronProps['name'] = 'KhaliqRaman'
        self._neuronProps['description'] = 'Khaliq-Raman Purkinje neuron model'
        
        self._synapseProps = synapseProps
        if self._synapseProps['name'].lower() == 'ampa':
            self._synapseProps['description'] = 'AMPA synapse'
            self._synapseProps['tauRise'] = -1.   # not used
            self._synapseProps['tauDecay'] = -1.  # not used
        elif self._synapseProps['name'].lower() == 'biexp':
            self._synapseProps['description'] = 'Bi-exponential synapse'
        else:
            del self._synapseProps
            raise Exception('Unknown synaptic model [%s]' % synapseProps['name'])
        
        self._prcProps = prcProps
        
        Neuron.__init__(self, ID, verbose)

    @property
    def length(self):
        return self._neuronProps['length']

    @property
    def diameter(self):
        return self._neuronProps['diameter']

    @property
    def nSynapses(self):
        return self._neuronProps['nSynapses']

    def synapseCluster(self, size):
        for i in range(self.nSynapses/size):
            yield i*size + np.arange(size)

    def _build(self):
        if self.verbose:
            print('>>> Building a Khaliq-Raman neuron model.')
        self._soma = h.Section()
        self._soma.insert('naRsg')
        self._soma.insert('kpkj')
        self._soma.insert('kpkj2')
        self._soma.insert('kpkjslow')
        self._soma.insert('bkpkj')
        self._soma.insert('cadiff')
        self._soma.insert('cap')
        self._soma.insert('lkpkj')
        self._soma.insert('hpkj')
        self._soma.L = self.length
        self._soma.diam = self.diameter
        self._soma.ena = 60
        self._soma.ek = -88

    def _addSynapses(self):
        """
        Adds all synapses to the neuron. Being this model composed by just a single
        compartment, all synapses are added to the soma.
        """
        if self.verbose:
            print('>>> Adding synapses to the model.')
        # TODO: take synapse parameters from:
        # Jaeger, D., De Schutter, E., & Bower, J. M. (1997).
        # The role of synaptic and voltage-gated currents in the control of Purkinje cell
        # spiking: a modeling study.
        # The Journal of Neuroscience, 17(1), 91-106.
        self.soma.push()
        for i in range(self.nSynapses):
            if self.synapseProperties['name'].lower() == 'ampa':
                syn = h.AMPA_S(self.soma(0.5))
                h.Erev_AMPA_S = self.synapseProperties['Erev']
            elif self.synapseProperties['name'].lower() == 'biexp':
                syn = h.Exp2Syn(self.soma(0.5))
                syn.tau1 = self.synapseProperties['tauRise']
                syn.tau2 = self.synapseProperties['tauDecay']
                syn.e = self.synapseProperties['Erev']
            else:
                raise Exception('Unknown synaptic model [%s]' % synapseProps['name'])
                return
            self.synapses.append(syn)
        h.pop_section()
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


def saveNeuron(filename, neuron, saveVoltage=False, simulationName = 'Common input simulation'):
    fid = h5.H5File(filename, 'a', simulationName)
    groupName = '/Neurons/ID_' + str(neuron.ID)
    fid.createGroup('/','Neurons')
    fid.createGroup('/Neurons','ID_' + str(neuron.ID))
    fid.writeTable(groupName, 'Cell', neuron.neuronPropertiesClass, 'Neuron properties', neuron.neuronProperties)
    fid.writeTable(groupName, 'Synapses', neuron.synapsePropertiesClass, 'Synapses properties', neuron.synapseProperties)

    try:
        fid.writeTable(groupName,'PRC',neuron.phaseResponseCurvePropertiesClass,
                       'Phase Response Curve Properties',neuron.phaseResponseCurveProperties)
    except:
        print('Could not record PRC properties.')
    try:
        fid.writeArray(groupName, 'Perturbations', tbl.Float64Atom(), neuron.perturbationTimes())
    except:
        print('Did not record perturbations.')

    fid.writeArray(groupName, 'Spikes', tbl.Float64Atom(), neuron.spikeTimes())

    if saveVoltage:
        try:
            fid.writeArray(groupName, 'Voltage', tbl.Float64Atom(), neuron.somaticVoltage())
        except:
            pass
    fid.close()


    
