import CPGNetworkSimulator as nsim
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib as npml
from scipy.stats import circstd, circmean
import h5py
from scoop import futures

 
simulator = nsim.CPGNetworkSimulator("./models/MLR_60.txt",["a","b"],(["RGF_NaP_L_hind", "RGF_NaP_R_hind", "RGF_NaP_L_front", "RGF_NaP_R_front"],))

variable_names = ('d0_CnF_GAT','d0_CnF_Glu')
ranges = ([1.15,1.85],[2.34,2.64])
ofilename='cnf_glu.h5'
steps = (100,100)

dt=0.001
duration = 10.0

alpha=0.0
simulator.setAlpha(alpha)
time_vec = np.arange(0.0,duration,dt)

v0=ranges[0][0]+(np.arange(0,steps[0],1))/(steps[0]-1.0)*(ranges[0][1]-ranges[0][0])
v1=ranges[1][0]+(np.arange(0,steps[1],1))/(steps[1]-1.0)*(ranges[1][1]-ranges[1][0])
#import IPython; IPython.embed()

def run_sim():
    out = np.zeros((len(time_vec),4))
    for ind_t,t in enumerate(time_vec):
        simulator.step(dt)
        act = simulator.getAct()
        out[ind_t,:]=act[0]
    return out

def calc_phase(time_vec,out):
    os_=((np.diff((out>0.1).astype(np.int),axis=0)==1).T)
    onsets=npml.repmat(time_vec[:-1],4,1)[os_]
    leg=(npml.repmat(np.arange(4.0),len(time_vec)-1,1).T)[os_]
    times=np.stack((onsets,leg),1)
    times=times[times[:,0].argsort()]

    pdur=np.diff(times[(times[:,1]==0),0])
    phases=np.zeros((3))
    std_phases=np.zeros((3))
    for i,(x,y) in enumerate([(0,1),(0,2),(0,3)]):
        times_=times[(times[:,1]==x) | (times[:,1]==y)]
        indices = np.where((times_[:-2,1]==x) & (times_[1:-1,1]==y) & (times_[2:,1]==x))
        pdur_=times_[[ind+2 for ind in indices],0]-times_[indices,0]
        phases_=((times_[[ind+1 for ind in indices],0]-times_[indices,0])/pdur_).T
        phases[i]=circmean(phases_[-6:-1],1,0)
        std_phases[i]=circstd(phases_[-6:-1],1,0)

    fq = 1/np.mean(pdur[-6:-1])
    return (fq,phases,np.max(std_phases))

def do_iteration(i,j):
    simulator.updateVariable(variable_names[0],v0[i])
    simulator.updateVariable(variable_names[1],v1[j])
    fq = 0
    phases_ = [] 
    stdp = 1.0
    while stdp>0.05:
        out = run_sim()
        fq, phases_, stdp = calc_phase(time_vec,out)
    return (fq,phases_)

simulator.updateVariable(variable_names[0],v0[0])
simulator.updateVariable(variable_names[1],v1[0])
run_sim()
IC=simulator.getState()


def do_one_bifurcation(i):
    frequency = np.zeros((steps[1],2))*np.nan
    phases = np.zeros((steps[1],3,2))*np.nan
    IChist=list()
    j_start_back=0
    go_up_on_nan=True
    simulator.setState(IC)
    simulator.updateVariable(variable_names[0],v0[i])
    simulator.updateVariable(variable_names[0],v1[0])
    run_sim()
    for j in range(0,steps[1]):
        IChist.append(simulator.getState())
        fq, phases_ = do_iteration(i,j)
        frequency[j,0]=fq
        phases[j,:,0]=phases_
        j_start_back=j
        if np.isnan(fq) and not go_up_on_nan:
            j_start_back=j-2
            break
        if not np.isnan(fq):
            go_up_on_nan=False
    
    simulator.setState(IChist[-2])
    for j in np.arange(j_start_back,-1,-1):
        fq, phases_ = do_iteration(i,j)
        frequency[j,1]=fq
        phases[j,:,1]=phases_
        if np.isnan(fq):
            break

    return (frequency,phases)
    
if __name__ == "__main__":
    frequency = np.zeros(steps+(2,))*np.nan
    phases = np.zeros(steps+(3,2))*np.nan
    paras=futures.map(do_one_bifurcation,[x for x in range(steps[0])])
    for i,(fq,ph) in enumerate(paras):
        frequency[i,:,:]=fq
        phases[i,:,:,:]=ph
    #import IPython; IPython.embed()
        
    h5f = h5py.File(ofilename, 'w')
    h5f.create_dataset('frequency', data=frequency)
    h5f.create_dataset('phases', data=phases)
    h5f.create_dataset('v0', data=v0)
    h5f.create_dataset('v1', data=v1)
    h5f.close()

    fig, ax = plt.subplots()
    plt.subplot(2,1,1)
    plt.pcolor(v0,v1,frequency[:,:,0])
    plt.subplot(2,1,2)
    plt.pcolor(frequency[:,:,1])

    fig, ax = plt.subplots()
    plt.subplot(3,2,1)
    plt.pcolor(v0,v1,phases[:,:,0,0])
    plt.subplot(3,2,3)
    plt.pcolor(v0,v1,phases[:,:,1,0])
    plt.subplot(3,2,5)
    plt.pcolor(v0,v1,phases[:,:,2,0])
    plt.subplot(3,2,2)
    plt.pcolor(v0,v1,phases[:,:,0,1])
    plt.subplot(3,2,4)
    plt.pcolor(v0,v1,phases[:,:,1,1])
    plt.subplot(3,2,6)
    plt.pcolor(v0,v1,phases[:,:,2,1])
    plt.show()   


#plt.plot(time_vec,out[:,:])
#plt.show()