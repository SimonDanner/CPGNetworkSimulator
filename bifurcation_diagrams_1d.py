import CPGNetworkSimulator as nsim
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib as npml
from scipy.stats import circstd, circmean
import h5py
from scoop import futures

neurons = ["RGF_NaP_L_hind", "RGF_NaP_R_hind", "RGF_NaP_L_front", "RGF_NaP_R_front"]
simulator = nsim.CPGNetworkSimulator("./models/MLR_45.txt",["a","b"],(neurons,))

sim='cnf_glu_ppn_i'
run=True

if sim=='cnf_glu':
    variable_name = 'd0_CnF_Glu'
    range_ = [3.95,5.3675]
    ofilename='cnf_glu_1d.h5'

if sim=='cnf_glu_ppn_i':
    variable_name = 'd0_CnF_Glu'
    range_ = [4.2,6.8775]
    ofilename='cnf_glu_ppn_i_1d.h5'
    simulator.updateVariable('PPN_Glu_to_LGPi_Glu1',  0.0)
    simulator.updateVariable('PPN_Glu_to_cLGPi_Glu1', 0.0)
    simulator.updateVariable('PPN_GAT_to_PPN_Glu',    0.0)

steps = 2000

dt=0.001
duration = 10.0
phase_diffs = [(0,1),(2,3),(0,2),(0,3),(1,3),(1,2)]

simulator.setAlpha(0.0)
time_vec = np.arange(0.0,duration,dt)

v=range_[0]+(np.arange(0,steps,1))/(steps-1.0)*(range_[1]-range_[0])

simulator.updateVariable(variable_name,v[0])


#import IPython; IPython.embed()

def run_sim():
    out = np.zeros((len(time_vec),len(neurons)))
    for ind_t,t in enumerate(time_vec):
        simulator.step(dt)
        act = simulator.getAct()
        out[ind_t,:]=act[0]
    return out

run_sim()
IC=simulator.getState()

def calc_phase(time_vec,out):
    os_=((np.diff((out>0.1).astype(np.int),axis=0)==1).T)
    onsets=npml.repmat(time_vec[:-1],len(neurons),1)[os_]
    leg=(npml.repmat(np.arange(len(neurons)),len(time_vec)-1,1).T)[os_]
    times=np.stack((onsets,leg),1)
    times=times[times[:,0].argsort()]

    pdur=np.diff(times[(times[:,1]==0),0])
    phases=np.zeros((len(phase_diffs)))
    std_phases=np.zeros((len(phase_diffs)))
    for i,(x,y) in enumerate(phase_diffs):
        times_=times[(times[:,1]==x) | (times[:,1]==y)]
        indices = np.where((times_[:-2,1]==x) & (times_[1:-1,1]==y) & (times_[2:,1]==x))
        pdur_=times_[[ind+2 for ind in indices],0]-times_[indices,0]
        phases_=((times_[[ind+1 for ind in indices],0]-times_[indices,0])/pdur_).T
        phases[i]=circmean(phases_[-6:-1],1,0)
        std_phases[i]=circstd(phases_[-6:-1],1,0)

    fq = 1/np.mean(pdur[-6:-1])
    return (fq,phases,np.max(std_phases))

def do_iteration(j):
    simulator.updateVariable(variable_name,v[j])
    fq = 0
    phases_ = [] 
    stdp = 1.0
    its=0
    while stdp>0.005 and its <20:
        out = run_sim()
        fq, phases_, stdp = calc_phase(time_vec,out)
        its+=1
    return (fq,phases_)


def do_one_bifurcation():
    frequency = np.zeros((steps,2))*np.nan
    phases = np.zeros((steps,len(phase_diffs),2))*np.nan
    IChist=list()
    j_start_back=0
    go_up_on_nan=True
    simulator.setState(IC)

    simulator.updateVariable(variable_name,v[0])
    for r in range(10):
        run_sim()
    for j in range(0,steps):
        IChist.append(simulator.getState())
        fq, phases_ = do_iteration(j)
        frequency[j,0]=fq
        phases[j,:,0]=phases_
        j_start_back=j-1
        if np.isnan(fq) and not go_up_on_nan:
            j_start_back=j-2
            break
        if not np.isnan(fq):
            go_up_on_nan=False
    
    simulator.setState(IChist[-2])
    for j in np.arange(j_start_back,-1,-1):
        fq, phases_ = do_iteration(j)
        frequency[j,1]=fq
        phases[j,:,1]=phases_
        if np.isnan(fq):
            break
    return (frequency,phases)
    
if __name__ == "__main__":
    
    if run==True:
        fq, ph = do_one_bifurcation()
        h5f = h5py.File(ofilename, 'w')
        h5f.create_dataset('frequency', data=fq)
        h5f.create_dataset('phases', data=ph)
        h5f.create_dataset('v', data=v)
        h5f.close()
    else:
        h5f = h5py.File(ofilename, 'r')
        fq = h5f['frequency'][()]
        v = h5f['v'][()]
        ph = h5f['phases'][()]
        h5f.close()

    
    
    _, (ax1, ax2, ax3, ax4,ax5) = plt.subplots(5, 1, sharex='all')
    
    ax1.plot(v,fq[:,0],'b',linewidth=1,)
    ax1.plot(v,fq[:,1],'r',linewidth=1)
    ax1.set_ylim([0.0, 12.0])
    ax2.plot(v,ph[:,0,0],'b.',markersize=1.5)
    ax2.plot(v,ph[:,0,1],'r.',markersize=1.5)
    ax2.plot(v,1-ph[:,0,0],'b.',markersize=1.5)
    ax2.plot(v,1-ph[:,0,1],'r.',markersize=1.5)
    ax2.set_ylim([-0.05, 1.05])

    ax3.plot(v,ph[:,1,0],'b.',markersize=1.5)
    ax3.plot(v,ph[:,1,1],'r.',markersize=1.5)
    ax3.plot(v,1-ph[:,1,0],'b.',markersize=1.5)
    ax3.plot(v,1-ph[:,1,1],'r.',markersize=1.5)
    ax3.set_ylim([-0.05, 1.05])
    
    ax4.plot(v,ph[:,2,0],'b.',markersize=1.5)
    ax4.plot(v,ph[:,2,1],'r.',markersize=1.5)
    ax4.plot(v,ph[:,4,0],'b.',markersize=1.5)
    ax4.plot(v,ph[:,4,1],'r.',markersize=1.5)
    ax4.set_ylim([-0.05, 1.05])

    ax5.plot(v,ph[:,3,0],'b.',markersize=1.5)
    ax5.plot(v,ph[:,3,1],'r.',markersize=1.5)
    ax5.plot(v,ph[:,5,0],'b.',markersize=1.5)
    ax5.plot(v,ph[:,5,1],'r.',markersize=1.5)
    ax5.set_ylim([-0.05, 1.05])
    
    
    plt.show()   


#plt.plot(time_vec,out[:,:])
#plt.show()