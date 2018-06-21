import sys
sys.path.insert(0, './build')
import CPGNetworkSimulator as nsim
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

left = ["MnIP_L_hind", "MnGM_L_hind", "MnVL_L_hind", "MnTA_L_hind", "MnSO_L_hind", "MnBF_L_hind","MnGA_L_hind"]
right = ["MnIP_R_hind", "MnGM_R_hind", "MnVL_R_hind", "MnTA_R_hind", "MnSO_R_hind", "MnBF_R_hind","MnGA_R_hind"]
mnnames = list()
mnnames.append(left)
mnnames.append(right)
simulator = nsim.CPGNetworkSimulator("./models/4CPG9MN_20180614_FB_LPN_FB.txt",mnnames)
lscondL = nsim.LimbSensorCondition(7)
lscondR = nsim.LimbSensorCondition(7)

simulator.setLscond([lscondL,lscondR])

t = np.arange(0.0,10,0.001)
acts00=np.zeros((len(t),7))
simulator.updateVariable("fbIa_IP",1.0)
alpha = 0.0

for i in range(len(t)):
    alpha+=0.0001
    simulator.setAlpha(alpha)
    simulator.step(0.001)
    if (t[i]==5.0):
        simulator.updateVariable("FtoIP",10.0)
        print("updated variable")

    if (t[i]==3.0):
        lscondL.Ia=[1.0,0,0,0,0,0,0]
        simulator.setLscond([lscondL,lscondR])
        print("updated lscond")


    act = simulator.getAct()
    acts00[i]=act[0][:]

fig, ax = plt.subplots()
plt.plot(t,acts00[:,0])
plt.show()