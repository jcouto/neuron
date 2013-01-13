from neuron import h
import matplotlib.pyplot as plt 
import numpy as np

h('{load_file("pyr3_template")}') # template from Traub et al. 2003
h('{load_file("stdrun.hoc")}')		# to make run() accessible...

def pyr3(x0=0, y0=0, z0=0):
	try:
		return h.pyr3(44, x0, y0, z0)
	except:
		print "Could not make pyramidal neuron..."
		return None

def insertRecorders(segment, labels, rec=None):
    ''' Inserts recorders for NEURON state variables. 
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

def insertIClamp(segment,amplitude=0,duration=1e60,delay=0):
    ''' Inserts an IClamp to segment. '''
    iclamp          = h.IClamp(segment)
    iclamp.amp      = amplitude
    iclamp.dur      = duration
    iclamp.delay    = delay
    print "Neuron: Inserted Iclamp, amplitude: "+str(amplitude)+"nA, duration: "+str(duration)+"ms, delay: "+str(delay)+"ms."
    return iclamp

#############################
# General parameters
plotFlag	  = True
N			  = 10			# number of cells
cell_distance = 150			# distance between cells (on x axis)
h.v_init	  = -70
tstop		  = 500.0

# Build model
cells		  = []
for ii in range(0,N):
	cells.append(pyr3(ii*cell_distance))
	
# Insert somatic Iclampstruct in each section
iClmp = []
jj=0
for ii in cells:
	iClmp.append(insertIClamp(ii.comp[1](0.5),amplitude = 0.03*jj, duration=tstop,delay = 3))
	jj+=1

# Insert recorders
jj = 0
rec = None
for ii in cells:
	rec = insertRecorders(ii.comp[1](0.5),{'cell'+str(jj):'_ref_v'},rec)
	jj+=1	

# Run
h.init()
h.tstop = tstop
h.run()

# Plot in stack
if plotFlag:
	for ii in range(0,len(cells)):
		plt.plot(rec['t'],np.array(rec['cell'+str(ii)])+80.0*ii, label='cell'+str(ii))
	ax = plt.gca()
	ax.legend(loc=1, borderaxespad=0.)
	plt.show()



