import CPGNetworkSimulator as nsim
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib as npml
import os

neurons = ["RGF_NaP_L", "RGF_NaP_R", "RGE_NaP_L", "RGE_NaP_R"]
filename = "./models/Model_2CPG.txt"
dt = 0.001
duration = 10.0
alphainit = 0.1

fn = os.path.join(os.path.dirname(__file__),".",filename)
sim = nsim.CPGNetworkSimulator(fn,["a","b"],(neurons,))
sim.setAlpha(alphainit)

time_vec = np.arange(0.0,duration,dt)

out = np.zeros((len(time_vec),len(neurons)*3))
for ind_t,t in enumerate(time_vec):
    sim.dense_step(dt)
    act = sim.getAct()
    Iepsp = sim.getIepsp()
    Iipsp = sim.getIipsp()
    out[ind_t,:len(neurons)]=act[0]
    out[ind_t,len(neurons):len(neurons)*2]=Iepsp[0]
    out[ind_t,len(neurons)*2:]=Iipsp[0]

_, ax = plt.subplots(3, 1, sharex='all')
ax[0].plot(time_vec,out[:,0:2])
ax[1].plot(time_vec,out[:,4:6])
ax[2].plot(time_vec,out[:,8:10])
plt.show()
