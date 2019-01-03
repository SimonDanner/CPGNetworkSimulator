import tools.py_simulator as nsim
from tools.plt import plot_1d_bifurcation_diagram, plot_2d_bifurcation_diagram, plot_noise_phase_fq, plot_noise_gait_fq
from optparse import OptionParser
import matplotlib as mpl_tmp
mpl_tmp.use('WXCairo')

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import numpy as np

if __name__ == "__main__":
    neurons = ["RGF_NaP_L_hind", "RGF_NaP_R_hind",
               "RGF_NaP_L_front", "RGF_NaP_R_front"]
    neurons = ["RGF_NaP_L", "RGF_NaP_R",
               "RGE_NaP_L", "RGE_NaP_R","CINe2_R"]
    filename = "./models/V3test.txt"

    filename = "./models/Model_2CPG.txt"


    #Note: the resolution of all calculations has been reduced to speed up simulation time.
    #Change following three valuse to increase number of steps
    steps_1D = 30
    steps_2D = (20,20)
    steps_noise = 106

    parser = OptionParser()
    parser.add_option("-s", "--simulation", dest="sim",default=0,type=int)
    options, args = parser.parse_args()

    cpg_sim = nsim.simulator(neurons=neurons, filename=filename)
    
    if options.sim == 1:
        cpg_sim.initialize_simulator()
        cpg_sim.sim.setAlpha(0.3)
        cpg_sim.updateVariable('LR_i',0.25)
        v, fq, ph, gait = cpg_sim.do_1d_bifurcation(
                'LR_e2', [0.1, 5.0], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)

    elif options.sim == 11:
        cpg_sim.alphainit=0.5
        v0, v1, fq, ph, gaits = cpg_sim.do_2d_bifurcation( ('LR_e2', 'alpha'), ([0.1,5.0], [0.01, 1.0]), steps_2D)
        plot_2d_bifurcation_diagram(v0, v1, fq, ph)

    elif options.sim == 21:
        col_out = np.zeros((0,len(neurons)))
        cpg_sim.initialize_simulator()
        cpg_sim.updateVariable('LR_e2',3.0)
        for alpha in np.arange(0.02,1.0,0.05):
            cpg_sim.sim.setAlpha(alpha)
            out = cpg_sim.run_sim()
            col_out = np.concatenate((col_out,out[-1000:,:]))
        plt.plot(col_out[:,-2:])
        plt.show()
        print(out)

    elif options.sim == 31:
        cpg_sim = nsim.simulator(neurons=neurons, filename=filename,duration=3.0)
        col_out = np.zeros((0,len(neurons)))
        cpg_sim.initialize_simulator()
        cpg_sim.updateVariable('LR_i',0.7)
        cpg_sim.updateVariable('LR_e2',3.0)
        cpg_sim.sim.setAlpha(0.15)
        for i in range(20):
            cpg_sim.run_sim()
        out = cpg_sim.run_sim()
        col_out = np.concatenate((col_out,out))
        cpg_sim.updateVariable('dCINe2L',20.0)
        cpg_sim.updateVariable('dCINe2R',20.0)
        out = cpg_sim.run_sim()
        col_out = np.concatenate((col_out,out))
        cpg_sim.updateVariable('dCINe2L',0.0)
        out = cpg_sim.run_sim()
        col_out = np.concatenate((col_out,out))
        

        cpg_sim.updateVariable('dCINe2R',00.0)
        out = cpg_sim.run_sim()
        col_out = np.concatenate((col_out,out))

        _, ax = plt.subplots(1, 1, sharex='all')
        
        
        pc1 = PatchCollection([Rectangle((3000, -3.1), 3000, 4)], facecolor='k', alpha=0.25,
                         edgecolor='w')
        pc2 = PatchCollection([Rectangle((6000, -3.1), 3000, 4)], facecolor='k', alpha=0.125,
                         edgecolor='w')
        # Add collection to axes
        ax.add_collection(pc1)
        ax.add_collection(pc2)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        # place a text box in upper left in axes coords
        ax.text(4000, -3.0, "stimulation: all V3",  fontsize=10,verticalalignment='top', bbox=props)
        ax.text(6500, -3.0, "stimulation: only right V3",  fontsize=10,verticalalignment='top', bbox=props)
        for i in range(4):
            plt.plot(col_out[:,i]-1.0*i)
        plt.yticks(np.arange(-3,1,1), ('RGE_R','RGE_L','RGF_R','RGF_L'))
        plt.show()
        print(out)