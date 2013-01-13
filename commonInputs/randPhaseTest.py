from common_input import *
h.load_file('stdrun.hoc')
import matplotlib.pylab as plt

def main(N=4):
    nrn=[]
    [nrn.append(KhaliqRaman()) for i in range(0,N)]
    h.finitialize()
    vclamps=[]
    vrest=-65
    for nn in nrn:
        nn.addSomaticVoltageRecorder()
        vclamps.append(h.VClamp(nn.soma(0.5)))
    for vclamp in vclamps:
        vclamp.amp[0] = -65
        vclamp.dur[0] = np.random.uniform(40)
    time = h.Vector()
    time.record(h._ref_t)
    h.tstop = 100
    h.run(100)
    for nn in nrn:
        plt.plot(time,nn.somaticVoltage)
    plt.show()

if __name__ == '__main__':
    main()
