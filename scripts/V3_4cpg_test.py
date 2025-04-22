"""
mlr_simulations.py

Script to calculate bifurcation diagrams of Danner SM, Shevtsova NA, Frigon A, Rybak IA. 
Computational modeling of spinal circuits controlling limb coordination and gaits in 
quadrupeds. eLife. 2017 Nov 22;6. pii: e31050. doi: 10.7554/eLife.31050. 

Use option -s to indicate simulation to be run.
"""
import CPGNetworkSimulator.tools.py_simulator as nsim 
from CPGNetworkSimulator.tools.plt import (plot_1d_bifurcation_diagram,
                       plot_2d_bifurcation_diagram, plot_noise_gait_fq,
                       plot_noise_phase_fq)
#import tools.py_simulator as nsim
#from tools.plt import plot_1d_bifurcation_diagram, plot_2d_bifurcation_diagram, plot_noise_phase_fq, plot_noise_gait_fq
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    neurons = ["RGF_NaP_L_hind", "RGF_NaP_R_hind",      # neurons to be read every time step 
               "RGF_NaP_L_front", "RGF_NaP_R_front"]
    filename = "./models/4CPG_2019_V3.txt" #  network model configuration file 
    #filename = "./models/eLife-V3.txt"
    steps_1D = 100       #number of steps for 1D-bifurcation diagrams
    steps_2D = (50,20) 

    parser = OptionParser()
    parser.add_option("-s", "--simulation", dest="sim",default=0,type=int)
    options, args = parser.parse_args()
    vs = [('V3toRGE', 2.0), ('V3toInE', 1.0), ('InEtoRGF', 0.8)]
    cpg_sim = nsim.simulator(neurons=neurons, filename=filename) # instantiate simulator

    if options.sim == 0:
        cpg_sim.initialize_simulator()
        cpg_sim.its_limit=1
        #cpg_sim.sim.updateVariable(vs[0][0],  vs[0][1])
        #cpg_sim.sim.updateVariable(vs[1][0],  vs[1][1])
        v, fq, ph, gait = cpg_sim.do_1d_bifurcation('alpha', [0.00, 1.05], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)
        plt.show()

    elif options.sim == 1:
        v0, v1, fq, ph, gaits = cpg_sim.do_2d_bifurcation(
            ('alpha', 'V3toInE'), ([0.02, 1.05], [0.5, 5.0]), steps_2D)
        plot_2d_bifurcation_diagram(v0, v1, fq, ph)

    elif options.sim == 2:
        cpg_sim.initialize_simulator()
        cpg_sim.sim.updateVariable(vs[0][0],  vs[0][1])
        for value in np.arange(1.0,2,0.11):
            cpg_sim.sim.updateVariable(vs[1][0],  value)
            v, fq, ph, gait = cpg_sim.do_1d_bifurcation('alpha', [0.02, 1.05], 50)
            plot_1d_bifurcation_diagram(v, fq, ph, gait)
            plt.savefig(vs[1][0] + ' ' + str(value)+'.png')

    elif options.sim == 3:
        cpg_sim.initialize_simulator()
        cpg_sim.sim.updateVariable(vs[0][0],  vs[0][1])
        cpg_sim.sim.updateVariable(vs[1][0],  vs[1][1])
        for value in np.arange(0.01,0.8,0.1):
            cpg_sim.sim.updateVariable(vs[2][0],  value)
            v, fq, ph, gait = cpg_sim.do_1d_bifurcation('alpha', [0.02, 1.05], 50)
            plot_1d_bifurcation_diagram(v, fq, ph, gait)
            plt.savefig(vs[2][0] + ' ' + str(value)+'.png')
