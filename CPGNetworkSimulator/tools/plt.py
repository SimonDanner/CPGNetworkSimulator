import matplotlib.pyplot as plt
import numpy as np

def plot_1d_bifurcation_diagram(v,fq,ph,gait,fqmax=12.0):
    _, (ax1, ax2, ax3, ax4,ax5,ax6) = plt.subplots(6, 1, sharex='all')

    ax1.plot(v,fq[:,0],'b',linewidth=1,)
    ax1.plot(v,fq[:,1],'r',linewidth=1)
    ax1.set_ylim([0.0, fqmax])
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
    
    #plt.show()  

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