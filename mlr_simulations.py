import tools.py_simulator as nsim
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib as npml
from scipy.stats import circstd, circmean
import h5py
import seaborn as sns
import pandas as pd

def plot_1d_bifurcation_diagram(v,fq,ph,gait):
    _, (ax1, ax2, ax3, ax4,ax5,ax6) = plt.subplots(6, 1, sharex='all')

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

    ax6.plot(v,gait[:,0],'b.',markersize=1.5)
    ax6.plot(v,gait[:,1],'r.',markersize=1.5)
    ax6.set_yticks([1,2,3,4])
    ax6.set_yticklabels(['walk','trot','gallop','bound'])
    ax6.set_ylim([.25, 4.75])
    
    plt.show()  

def plot_2d_bifurcation_diagram(v0,v1,frequency,phases):
    _, (ax1,ax2,ax3,ax4) = plt.subplots(4, 1, sharex='all')
    ax1.pcolor(v1,v0,frequency[:,:,0])
    ax2.pcolor(v1,v0,phases[:,:,0,0])
    ax3.pcolor(v1,v0,phases[:,:,2,0])
    ax4.pcolor(v1,v0,phases[:,:,3,0])
    plt.show()

def plot_noise_phase_fq(fqs,phases):
    fqedges = np.linspace(0.0,14,14*4+1)
    phedges = np.linspace(0.0,1.0,65)
    H = np.zeros((fqedges.size-1,phedges.size-1,4))
    _, ax = plt.subplots(4, 1, sharex='all')
    
    for i in range(4):
        for (fq,ph) in zip(fqs,phases):
            H_, _, _ = np.histogram2d(fq,ph[:,i],bins=(fqedges,phedges))
            H[:,:,i] += H_
            #print(_)     
        ax[i].pcolor(fqedges,phedges, H[:,:,i].T)
        


def plot_noise_gait_fq(fqs,gaits):
    fqedges = np.linspace(0.0,14,14*4+1)

    H = np.zeros((fqedges.size-1,4))
    _, ax1 = plt.subplots(1, 1, sharex='all')
    
    for i in range(4):
        for (fq,gait) in zip(fqs,gaits):
            H_, _ = np.histogram(fq[gait==i+1],fqedges)
            H[:,i] += H_
    
    H[:,2]+=H[:,3]
    H = H[:,0:3]
    ax1.pcolor(fqedges[:-1],[0.5,1.5,2.5,3.5], H[:,:].T)
    ax1.set_yticks([1,2,3])
    ax1.set_yticklabels(['walk','trot','gallop/bound'])

def plot_noise_phase_fq_kde(fqs,phases):
    fqs_=np.concatenate(fqs)
    ph_=np.concatenate(phases)
    f = pd.Series(fqs_,name='fq')
    lfh = pd.Series(ph_[:,0],name='lf-hind')
    lff = pd.Series(ph_[:,1],name='lf-fore')
    hl = pd.Series(ph_[:,2],name='homolateral')
    di = pd.Series(ph_[:,3],name='diagonal')
    
    sns.jointplot(f,lfh,kind="kde",xlim=(0,14.0),ylim=(-0.05,1.05),space=0)
    sns.jointplot(f,lff,kind="kde",xlim=(0,14.0),ylim=(-0.05,1.05),space=0)
    sns.jointplot(f,hl,kind="kde",xlim=(0,14.0),ylim=(-0.05,1.05),space=0)
    sns.jointplot(f,di,kind="kde",xlim=(0,14.0),ylim=(-0.05,1.05),space=0)
    plt.show()

def plot_noise_gait_fq_violin(ax,fqs,gaits):
    fqs_=np.concatenate(fqs)
    gaits_=np.concatenate(gaits)
    sns.violinplot(x=fqs_[gaits_!=0],y=pd.Categorical(gaits_[gaits_!=0]),ax=ax)
    ax.set_xlim([0.0, 14.0])

if __name__ == "__main__":
    neurons = ["RGF_NaP_L_hind", "RGF_NaP_R_hind", "RGF_NaP_L_front", "RGF_NaP_R_front"]
    filename = "./models/MLR_45.txt"

    cpg_sim = nsim.simulator(neurons=neurons,filename=filename)
    do_sim=22
    if do_sim==0:
        v,fq, ph, gait = cpg_sim.do_1d_bifurcation('d0_CnF_Glu_L',[3.95,5.3675],200)
        plot_1d_bifurcation_diagram(v,fq,ph,gait)   
        
    elif do_sim==1:
        cpg_sim.initialize_simulator()
        cpg_sim.sim.updateVariable('PPN_Glu_to_LGPi_Glu1',  0.0)
        cpg_sim.sim.updateVariable('PPN_Glu_to_cLGPi_Glu1', 0.0)
        cpg_sim.sim.updateVariable('PPN_GAT_to_PPN_Glu',    0.0)
        
        v,fq, ph, gait = cpg_sim.do_1d_bifurcation('d0_CnF_Glu_L',[4.2,6.8775],50)
        plot_1d_bifurcation_diagram(v,fq,ph,gait) 
    elif do_sim==2:
        v,fq, ph, gait = cpg_sim.do_1d_bifurcation('d0_LPGi_Glu',[2.45,2.45+1.1],50)
        plot_1d_bifurcation_diagram(v,fq,ph,gait)
    
    elif do_sim==11:
        v0,v1,fq,ph = cpg_sim.do_2d_bifurcation(('d0_CnF_GAT','d0_CnF_Glu_bl'),([1.1,3.1],[2.78,3.1]),(20,20))
        plot_2d_bifurcation_diagram(v0,v1,fq,ph)
    elif do_sim==12:
        v0,v1,fq,ph = cpg_sim.do_2d_bifurcation(('d0_PPN_GAT','d0_CnF_Glu_bl'),([1.1,3.0],[2.78,3.1]),(20,20))
        plot_2d_bifurcation_diagram(v0,v1,fq,ph)
    elif do_sim==13:
        v0,v1,fq,ph = cpg_sim.do_2d_bifurcation(('d0_LPGi_GAT','d0_CnF_Glu_bl'),([1.1,3.0],[2.78,3.1]),(20,20))
        plot_2d_bifurcation_diagram(v0,v1,fq,ph)

    elif do_sim==20:
        v,fqs,fl_phase_durs,ex_phase_durs,phases,gaits = cpg_sim.do_noise_simulation('d0_CnF_Glu_L',[3.95,5.3675],50,1.0)
        #plot_noise_phase_fq_kde(fqs,phases)
        #import IPython; IPython.embed()
        plot_noise_gait_fq(fqs,gaits)
        plot_noise_phase_fq(fqs,phases)
        plt.show()
    elif do_sim==21:
        v,fqs,fl_phase_durs,ex_phase_durs,phases,gaits = cpg_sim.do_noise_simulation('d0_PPN_Glu',[4.0,5.5],50,1.0)
        plot_noise_gait_fq(fqs,gaits)
        plot_noise_phase_fq(fqs,phases)
        
        plt.show()
    elif do_sim==22:
        sv=[('PPN_Glu_to_LGPi_Glu1',  0.0),('PPN_Glu_to_cLGPi_Glu1',  0.0),('PPN_GAT_to_PPN_Glu',  0.0)]
        v,fqs,fl_phase_durs,ex_phase_durs,phases,gaits = cpg_sim.do_noise_simulation('d0_CnF_Glu_L',[4.2,6.8775],50,1.0,sv)
        plot_noise_gait_fq(fqs,gaits)
        plot_noise_phase_fq(fqs,phases)
        plt.show()
    elif do_sim==23:
        v,fqs,fl_phase_durs,ex_phase_durs,phases,gaits = cpg_sim.do_noise_simulation('d0_LPGi_Glu',[2.45,2.45+1.1],20,1.0)
        plot_noise_gait_fq(fqs,gaits)
        plot_noise_phase_fq(fqs,phases)
        plt.show()
    elif do_sim==30:
        _, (ax1,ax2,ax3,ax4) = plt.subplots(4, 1, sharex='all')
        _,fqs1,_,_,_,gaits1 = cpg_sim.do_noise_simulation('d0_CnF_Glu_L',[3.95,5.3675],20,1.0)
        cpg_sim = nsim.simulator(neurons=neurons,filename=filename)
        _,fqs2,_,_,_,gaits2 = cpg_sim.do_noise_simulation('d0_PPN_Glu',[4.0,5.5],20,1.0)
        cpg_sim = nsim.simulator(neurons=neurons,filename=filename)
        sv=[('PPN_Glu_to_LGPi_Glu1',  0.0),('PPN_Glu_to_cLGPi_Glu1',  0.0),('PPN_GAT_to_PPN_Glu',  0.0)]
        _,fqs3,_,_,_,gaits3 = cpg_sim.do_noise_simulation('d0_CnF_Glu_L',[4.2,6.8775],50,1.0,sv)
        cpg_sim = nsim.simulator(neurons=neurons,filename=filename)
        _,fqs4,_,_,_,gaits4 = cpg_sim.do_noise_simulation('d0_LPGi_Glu',[2.45,2.45+1.1],20,1.0)
        cpg_sim = nsim.simulator(neurons=neurons,filename=filename)
        
        plot_noise_gait_fq_violin(ax1,fqs1,gaits1)
        plot_noise_gait_fq_violin(ax2,fqs2,gaits2)
        plot_noise_gait_fq_violin(ax3,fqs3,gaits3)
        plot_noise_gait_fq_violin(ax4,fqs4,gaits4)
        plt.show()
        
        

    #cpg_sim.sim.updateParameter('sigmaNoise',2.0)
