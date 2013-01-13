from neuron import h
#load session files
# this needs to be changed when using another computer than my mac!
h.load_file("/Users/joao/NeuronLibrary/KhaliqRamanPkj/naRsg_det.ses")
#create some variables 

global somaLength,somaDiam,soma
somaLength=20
somaDiam=20*13
def create_shape():
	global soma,somaLength,somaDiam
	print "Creating single neuron..."
	soma=h.Section()
	soma.L=somaLength
	soma.diam=somaDiam
	soma.nseg=1
	print "... and appending mechanisms..."
	soma.insert('naRsg_det')
	soma.insert('kpkj')
	soma.insert('kpkj2')
	soma.insert('kpkjslow')
	soma.insert('bkpkj')
	soma.insert('cadiff')
	soma.insert('cap')
	soma.insert('lkpkj')
	soma.insert('hpkj')
	for seg in soma:
		print "fixing ena and ek!"
		seg.ena=60
		seg.ek=-88


