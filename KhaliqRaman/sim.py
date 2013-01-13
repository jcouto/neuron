from neuron import h
import array
import model,electrodes
import numpy as np
import tables as tb
import os

# constants
constants={'simtime':10000,'baseline':0,'pulsedelay':2500,'pulseamp':0,'pulsedur':0.5,'gaussianNoise':0.1,'ModelName':'Khaliq_Raman deterministic','nTrials':1}

#create all shapes and add electrodes
model.create_shape()
electrodes.apply(model,constants)
electrodes.detect_spk(model)
#spks=np.zeros((constants['nTrials'],),dtype=np.object)
electrodes.baseline.amp=constants["baseline"]
# create data file
pathToData='/Users/joao/Documents/projectPurkinjePrc/modelling/newKhaliqRaman/influenceOfNoiseOnTheISI/data/'
folderlist=os.listdir(pathToData)
filename=pathToData+'data'+str(len(folderlist)+1)+'.h5' 
h5file = tb.openFile(filename, mode = "w",title = constants['ModelName'])
PRCt=h5file.createGroup("/",'PRCt','Phase response curve with the traditional method') #root group
for n in range(0,1):
	data=h5file.createGroup(PRCt,'ds'+str(n),'Dataset '+str(n)) #children groups create dataset
	# saving the parameters of the dataset
	par=h5file.createGroup(data,'parameters','Parameters of the simulation')
	h5file.createArray(par,'tstop',constants["simtime"],'Trial length (ms)')
	h5file.createArray(par,'holding',constants["baseline"],'Holding current (nA)')
	h5file.createArray(par,'pulseAmp',constants["pulseamp"],'Amplitude of the perturbation pulse (nA)')
	h5file.createArray(par,'pulseDel',constants["pulsedelay"],'Delay before perturbation (ms)')
	h5file.createArray(par,'pulseDur',constants["pulsedur"],'Duration of the perturbation pulse (ms)')
	h5file.createArray(par,'GnoiseStd',constants["gaussianNoise"],'Standard deviation of the injected gaussian noise (nA)')
	h5file.createArray(par,'nTrials',constants["nTrials"],'Number of trials')
	h5file.createArray(par,'modelName',constants["ModelName"],'Name/specification of the model')
	# run trials
	for ii in range(0,constants['nTrials']):
		# run the simulation (1 trial)
		h.load_file("stdrun.hoc")
		print "---> Trial " + str(ii)
		f=open('/dev/random')
		a=array.array('L')
		a.fromfile(f,1)
		f.close()
		useThisSeed=int(a[0]/10000)
		electrodes.set_random(useThisSeed)
		electrodes.gaussianNoise.std=constants["gaussianNoise"];
		electrodes.apc.n=0;
		h.dt=0.025
		h.tstop=constants['simtime']
		h.init()
		h.run(h.tstop)
		spks=(electrodes.rec_spk.to_python())
		h5file.createArray(data,'t'+str(ii),spks,'Trial '+str(ii))
		h5file.createArray(data,'t'+str(ii)+'_seed',useThisSeed,'Seed used in trial '+str(ii))
h5file.close()	

