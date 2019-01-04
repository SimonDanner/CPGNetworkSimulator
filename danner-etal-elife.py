"""
mlr_simulations.py

Script to calculate bifurcation diagrams of Danner SM, Shevtsova NA, Frigon A, Rybak IA. 
Computational modeling of spinal circuits controlling limb coordination and gaits in 
quadrupeds. eLife. 2017 Nov 22;6. pii: e31050. doi: 10.7554/eLife.31050. 

Use option -s to indicate simulation to be run.
"""

import tools.py_simulator as nsim
from tools.plt import plot_1d_bifurcation_diagram, plot_2d_bifurcation_diagram, plot_noise_phase_fq, plot_noise_gait_fq
from optparse import OptionParser
import matplotlib.pyplot as plt

if __name__ == "__main__":
    neurons = ["RGF_NaP_L_hind", "RGF_NaP_R_hind",      # neurons to be read every time step 
               "RGF_NaP_L_front", "RGF_NaP_R_front"]
    filename = "./models/Danner-etal-eLife.txt" #  network model configuration file 

    steps_1D = 50       #number of steps for 1D-bifurcation diagrams

    parser = OptionParser()
    parser.add_option("-s", "--simulation", dest="sim",default=0,type=int)
    options, args = parser.parse_args()

    cpg_sim = nsim.simulator(neurons=neurons, filename=filename) # instantiate simulator
    if options.sim == 0:
        '''Intact model '''
        v, fq, ph, gait = cpg_sim.do_1d_bifurcation('alpha', [0.02, 1.05], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)
    if options.sim == 1:
        '''Simulate ablation of all V0V '''
        cpg_sim.initialize_simulator()
        cpg_sim.sim.updateVariable('V0VtoRGFdiagfh',  0.0)
        cpg_sim.sim.updateVariable('V0VtoRGFdiaghf',  0.0)
        cpg_sim.sim.updateVariable('inV0VtoRGF',  0.0)
        v, fq, ph, gait = cpg_sim.do_1d_bifurcation('alpha', [0.02, 1.05], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)
    if options.sim == 2:
        '''Simulate ablation of diagonal V0V '''
        cpg_sim.initialize_simulator()
        cpg_sim.sim.updateVariable('V0VtoRGFdiagfh',  0.0)
        cpg_sim.sim.updateVariable('V0VtoRGFdiaghf',  0.0)
        v, fq, ph, gait = cpg_sim.do_1d_bifurcation('alpha', [0.02, 1.05], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)
    if options.sim == 3:
        '''Simulate ablation of all V0 (V0D and V0V) '''
        cpg_sim.initialize_simulator()
        cpg_sim.sim.updateVariable('V0VtoRGFdiagfh',  0.0)
        cpg_sim.sim.updateVariable('V0VtoRGFdiaghf',  0.0)
        cpg_sim.sim.updateVariable('inV0VtoRGF',  0.0)
        cpg_sim.sim.updateVariable('V0DtoRGFdiagfh',  0.0)
        cpg_sim.sim.updateVariable('V0DtoRGF',  0.0)
        v, fq, ph, gait = cpg_sim.do_1d_bifurcation('alpha', [0.02, 1.05], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)
    if options.sim == 4:
        '''Simulate ablation of descending long propriospinal neurons '''
        cpg_sim.initialize_simulator()
        cpg_sim.sim.updateVariable('V0VtoRGFdiagfh',  0.0)
        cpg_sim.sim.updateVariable('V0DtoRGFdiagfh',  0.0)
        cpg_sim.sim.updateVariable('inFH',  0.0)
        cpg_sim.sim.updateVariable('V2aHomfh',  0.0)
        v, fq, ph, gait = cpg_sim.do_1d_bifurcation('alpha', [0.02, 1.05], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)