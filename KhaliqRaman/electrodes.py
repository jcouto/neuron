from neuron import h
import array

constants={'simtime':10000,'baseline':0,'pulsedelay':300,'pulseamp':0,'pulsedur':0.5,'gaussianNoise':0.14}

def addBaseline(model,amplitude=0,duration=1e50,delay=0):
'''
    Adds a constant current with delay 0 to the soma in model.
    Amplitude in nA.
    Delay and duration in ms.
'''
	model.soma.push()
	global baseline
	baseline=h.IClamp(0.5)
	baseline.dur=duration
	baseline.amp=amplitude      # in nA!!
	baseline.delay=delay
    print "Baseline current is " + str(amplitude/1000) + "pA."

def addGaussianNoise(model,amplitude=0,duration=1e50,delay=0):	
 '''
    Adds gaussian distributed noise to introduce spike jitter.
    Amplitude in nA.
    Delay and duration in ms.
'''
    global gaussianNoise
	
	gaussianNoise=h.IClampNoise(0.5)
	gaussianNoise.delay=delay
	gaussianNoise.dur=duration
	gaussianNoise.std=amplitude
	
    print "Baseline current is " + str(amplitude/1000) + "pA."

def recordAll(model):
	global rec_t,rec_c,rec_v
	rec_t=h.Vector()
	rec_v=h.Vector()
	#rec_c=h.Vector()
	rec_t.record(h._ref_t)
	rec_v.record(model.soma(0.5)._ref_v)
	#rec_c.record(gaussianNoise._ref_i)
	
def detect_spk(model):
	''' inserts an APCount mechanism on the soma of the cell spescified in model'''
	global rec_spk
	rec_spk=h.Vector()
	apc= h.APCount(model.soma(0.5))
	apc.thresh = -10 #set threshold to a value of your choice
	apc.record(rec_spk)

def set_random(value): 
	''' Sets the specified value on the random seed... Doesn't do much.... '''
	global gaussianNoise
	gaussianNoise.seed(value)
	print "Random seed of the Gaussian Noise set to " + str(value) + "."


