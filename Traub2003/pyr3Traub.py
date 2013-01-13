import numpy as np
import h5py as h5
import sys
import os.path
import time
import getopt
from neuron import h
import ConfigParser as cp
h('{load_file("pyr3_template")}') # template from Traub et al. 2003
h('{load_file("stdrun.hoc")}')		# to make run() accessible...


def pyr3(x0=0, y0=0, z0=0,fig=0):
	try:
		return h.pyr3(fig, x0, y0, z0)
	except:
		print "Could not make pyramidal neuron..."
		return None

def getOptions(cfg,metadata,section,verb=True):
    '''Parse configuration file from default values.
	   "cfg" is a ConfigParser file object
	    "metadata" a dictionary with default parameters
		"section" is the section to look for in the cfg file
		This function looks for the same type of the value in the metadata dict.
			 If it is a float or a list, it will look for a value or evaluate an expression.'''
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
		
def insertRecorders(segment, labels, rec=None):
    ''' Inserts recorders for NEURON state variables. 
    Note: labels is a dictionary. Example {'v': '_ref_v'}.
          Specify 'rec' to append to previous recorders.
          Records also time if 'rec' is 'None'(default).
	rec = insertRecorders(cell.comp[1](0.5),{'cell'+str(cell_id):'_ref_v'}) 
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
    #print "Neuron: Inserted Iclamp, amplitude: "+str(amplitude)+"nA, duration: "+str(duration)+"ms, delay: "+str(delay)+"ms."
    return iclamp

def createNeuron(x0=0,y0=0,z0=0,fig=0, cell_id=0,additional=[]):
	cell=pyr3(x0,y0,z0,fig)	
	for ii in additional:
		if ii['key'] in ['IClamp']:
			iClmp = insertIClamp(cell.comp[1](0.5),amplitude = ii['amp'], duration=ii['dur'],delay = ii['delay'])
			ii['object'] = iClmp
	# record spikes
	nc =h.NetCon(cell.comp[1](0.5)._ref_v,None) 
	nc.delay = 0
	nc.threshold = 0
	rec = {'spks'+str(cell_id):h.Vector()}
	nc.record(rec['spks'+str(cell_id)])
	return [cell, rec, additional]

def saveDictToH5Attrs(ID,labels):
	for kk,oo in labels.items():
		if type(oo)==bool:
			ID.attrs[kk] = int(oo)
		else:
			try:
				ID.attrs[kk] = oo
			except:
				try:
					ds = ID.create_dataset(kk,shape=(len(oo),),dtype='f') 
					ds[:] = oo.flatten() 
				except:
					print("Could not record dataset: %s in object %s"%(kk,ID.name))

def saveDataToH5(pc,par,pariclamp,gids,recs):
	timeRun			= time.clock()-par['timeRun']
	timeRec			= time.clock()
	pc.barrier()	# wait until all workers are done.
	timeAbsRun		= time.clock()
	filename		= os.path.abspath(par['datadir'])+'/'+par['filename']
	if pc.id()==0:
		try:
			fid = h5.File(filename,'a')
			fid.create_group('data')
		except IOError:
			print("File (%s) could not be created!"%(filename))
			sys.exit("Exit without saving!!!")
		except ValueError:
			print("File (%s) already contained data!"%(filename))
			sys.exit("Exit without saving!!!")
		fid.close()
	pc.barrier()	# make sure the workers don't try to access the file when it hasn't been created 
	for ii in range(0,int(pc.nhost())):
		pc.barrier()
		if (pc.id()==ii):
			fid = h5.File(filename,'a')
			processor = fid.create_group("processor"+str(ii))
			processor.attrs['timeRun'] = timeRun
			processor.attrs['cellIDs']	   = gids
			data = fid['data']
			for jj in recs:
				for kk in jj.keys():
					d=np.array(jj[kk].to_python())
					try:
						ds = data.create_dataset(kk,shape=(np.size(d),),dtype='f',compression='gzip',maxshape=(None,)) 
						ds[:] = d[:]
						ds.attrs['name'] = kk
					except:
						if not 't'==kk:
							print "Could not save, the following datasets where found when recording "+kk+':'
							print data.keys()
			processor.attrs['timeRec'] = time.clock()-timeRec
			fid.close()
		pc.barrier()	# one worker at a time...
	if pc.id()==0:
		fid			= h5.File(filename,'a')
		data=fid['data']
		config		= fid.create_group('parameters')
		saveDictToH5Attrs(config,par)
		iclmp		= config.create_group('IClamp')
		saveDictToH5Attrs(iclmp,pariclamp)
		data.attrs['timeAbs'] = time.clock()-par['timeRun']
		data.attrs['timeRun'] = timeAbsRun - par['timeRun']
		fid.close()
		print("Saved data to file (%s)."%(filename))

def spikeout(pc,gids,recs):
	pc.barrier() #wait for hosts to get here
	if (pc.id()==0):
		print "cell (id)\t processor\t spiketimes (ms)"
	for rank in range(0,int(pc.nhost())):
		if rank == pc.id():
			for ii in range(len(gids)):
				#print "Looking for cell "+str(gids[ii])+" in processor "+str(rank)
				for jj in recs:
					try:
						for spk in np.array(jj["spks"+str(gids[ii])]):
							print str(gids[ii])+"\t"+str(rank)+'\t'+str(spk)
					except:
						pass
	pc.barrier()

# Config file under "Simulation".
g_par					= {}
g_par['N']				= 6
g_par['cellSpacing']	= 150.0
g_par['v_init']			= -70.0
g_par['tstop']			= 500.0
g_par['datadir']		= './'
g_par['filename']		= 'data.h5'
g_par['record']         = True
g_par['traubFigure']	= 44
 
# Config file under "IClamp".
g_iclamppar				= {}
g_iclamppar['id']		= np.arange(1,g_par['N'])
g_iclamppar['delay']	= np.ones((g_par['N'],1))*0.0
g_iclamppar['amp']		= np.linspace(0,0.5,g_par['N'])
g_iclamppar['dur']		= np.ones((g_par['N'],1))*g_par['tstop']



def main():
	###############PARSE SCRIPT OPTIONS#######################################
	print sys.argv[1:]
	counter = 0
	for ii in sys.argv:
		counter+=1
		if ii in [os.path.basename(__file__)]:
			break
	try:
		opts, args = getopt.gnu_getopt(sys.argv[counter:], "hc:d:", ["help", "cfg","datadir"])
	except getopt.GetoptError, err:
		sys.exit("Help: Usage: \n        -cfg cfgName \n        -datadir dataFolder.\n")
	dataFolder = None
	cfgName    = None
	par		   = g_par
	iclamppar  = g_iclamppar
	for o, a in opts:
		if o in ("-h","--help"):
			sys.exit("Help: Usage: \n        -cfg=cfgName \n        -datadir=dataFolder\n        Or nothing for model.cfg and CWD.\n")
		elif o in ("-c", "--cfg"):
			cfgName = a
			fid = open(cfgName)
			print "Parsing "+cfgName+"."
			cfg = cp.ConfigParser()
			cfg.readfp(fid)
			fid.close()
			par = getOptions(cfg,par,'Simulation',verb=False)
			iclamppar = getOptions(cfg,iclamppar,'IClamp',verb=False)
		elif o in ("-d","--datadir"):
			dataFolder = a
		else:
			assert False, "unhandled option"
	##########################################################################
    ###############LOADS PARAMETERS AND BUILDS THE MODEL######################
	# Start ParallelContext()
	par['timeRun'] = time.clock()
	pc = h.ParallelContext()
	st = pc.time()
	# Build model
	gids = []
	cells = []
	additionals = []
	recs = []
	# Splits cells accross processors.
	for ii in range(int(pc.id()),int(par['N']),int(pc.nhost())):
		additional = []
		if float(ii) in iclamppar['id']:
			pp = np.where(iclamppar['id']==float(ii))
			additional.append({'key':'IClamp','amp':iclamppar['amp'][pp],'dur':iclamppar['dur'][pp],'delay':iclamppar['delay'][pp]})
		[cell, rec, additional] = createNeuron(x0=ii*par['cellSpacing'],fig=par['traubFigure'],cell_id=ii, additional=additional)
		print "Cell id "+str(ii)+" on processor ",str(pc.id())
		cells.append(cell)
		additionals.append(additional)
		recs.append(rec)
		gids.append(ii)
		pc.set_gid2node(ii, pc.id()) # assign neuron to gid.
	# Run
	h.tstop = par['tstop']
	pc.set_maxstep(1)		# Sets every machine's maximum step size to the minimum delay of those netcons. Usefull if the netcons are connected between processors.
	h.stdinit()
	pc.psolve(par['tstop'])		# Run sim on all workers. Similar to cvode.solve(tstop)
	#spikeout(pc,gids,recs)
	if g_par['record']:
		saveDataToH5(pc, par, iclamppar,gids, recs)
	pc.runworker()			# Master node is now the only worker.	
	pc.done()
	h.quit()

if __name__ == "__main__":
    main()

