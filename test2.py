import CPGNetworkSimulator as nsim
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib as npml
import os

neurons = ["RGF_NaP_L", "RGF_NaP_R", "RGE_NaP_L", "RGE_NaP_R"]
filename = "./models/V3_Model.txt"
dt = 0.001
duration = 100.0
init_dur = 20.0

fn = os.path.join(os.path.dirname(__file__),".",filename)
sim = nsim.CPGNetworkSimulator(fn,["a","b"],(neurons,))
#sim.setAlpha(alphainit)

time_vec = np.arange(0.0,duration+init_dur,dt)

eleak0 = np.array(sim.getEleak())
alphas = (0.001,0.13)

out = np.zeros((len(time_vec),len(neurons)*3))
for ind_t,t in enumerate(time_vec):
    if t < init_dur:
        alpha = alphas[0]
    else:
        alpha = (t-init_dur)*(alphas[1]-alphas[0])/duration + alphas[0]

    sim.setEleak(eleak0*(1.0-alpha))
    sim.dense_step(dt)
    act = sim.getAct()
    Iepsp = sim.getIepsp()
    Iipsp = sim.getIipsp()
    out[ind_t,:len(neurons)]=act[0]
    out[ind_t,len(neurons):len(neurons)*2]=Iepsp[0]
    out[ind_t,len(neurons)*2:]=Iipsp[0]

_, ax = plt.subplots(3, 2, sharex='all',sharey='row')
ax[0,0].plot(time_vec,out[:,0:2])
ax[0,0].set_title('RG-F f(V)')
ax[1,0].plot(time_vec,out[:,4:6])
ax[1,0].set_title('Iepsp')
ax[2,0].plot(time_vec,out[:,8:10])
ax[2,0].set_title('Iipsp')
ax[0,1].plot(time_vec,out[:,2:4])
ax[0,1].set_title('RG-E f(V)')
ax[1,1].plot(time_vec,out[:,6:8])
ax[1,1].set_title('Iepsp')
ax[2,1].plot(time_vec,out[:,10:12])
ax[2,1].set_title('Iipsp')
plt.show()
