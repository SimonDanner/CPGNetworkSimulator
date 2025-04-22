"""
mlr_simulations.py

Script to calculate bifurcation diagrams of Danner SM, Shevtsova NA, Frigon A, Rybak IA. 
Computational modeling of spinal circuits controlling limb coordination and gaits in 
quadrupeds. eLife. 2017 Nov 22;6. pii: e31050. doi: 10.7554/eLife.31050. 

Use option -s to indicate simulation to be run.
"""

import CPGNetworkSimulator.tools.py_simulator as nsim
from CPGNetworkSimulator.tools.plt import plot_1d_bifurcation_diagram, plot_2d_bifurcation_diagram, plot_noise_phase_fq, plot_noise_gait_fq
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    neurons = ["RGF_NaP_L_hind", "RGF_NaP_R_hind",      # neurons to be read every time step 
               "RGF_NaP_L_front", "RGF_NaP_R_front"]
    filename = "./models/Danner-etal-eLife.txt" #  network model configuration file 

    steps_1D = 100       #number of steps for 1D-bifurcation diagrams

    parser = OptionParser()
    parser.add_option("-s", "--simulation", dest="sim",default=0,type=int)
    options, args = parser.parse_args()

    cpg_sim = nsim.simulator(neurons=neurons, filename=filename) # instantiate simulator
    if options.sim == 0:
        '''Intact model '''
        v, fq, ph,fl_dur,ex_dur, gait = cpg_sim.do_1d_bifurcation('alpha', [0.02, 1.05], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)
        plt.show()
    if options.sim == 1:
        '''Simulate ablation of all V0V '''
        cpg_sim.initialize_simulator()
        cpg_sim.sim.updateVariable('V0VtoRGFdiagfh',  0.0)
        cpg_sim.sim.updateVariable('V0VtoRGFdiaghf',  0.0)
        cpg_sim.sim.updateVariable('inV0VtoRGF',  0.0)
        v, fq, ph,fl_dur,ex_dur, gait = cpg_sim.do_1d_bifurcation('alpha', [0.02, 1.05], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)
    if options.sim == 2:
        '''Simulate ablation of diagonal V0V '''
        cpg_sim.initialize_simulator()
        cpg_sim.sim.updateVariable('V0VtoRGFdiagfh',  0.0)
        cpg_sim.sim.updateVariable('V0VtoRGFdiaghf',  0.0)
        v, fq, ph,fl_dur,ex_dur, gait = cpg_sim.do_1d_bifurcation('alpha', [0.02, 1.05], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)
    if options.sim == 3:
        '''Simulate ablation of all V0 (V0D and V0V) '''
        cpg_sim.initialize_simulator()
        cpg_sim.sim.updateVariable('V0VtoRGFdiagfh',  0.0)
        cpg_sim.sim.updateVariable('V0VtoRGFdiaghf',  0.0)
        cpg_sim.sim.updateVariable('inV0VtoRGF',  0.0)
        cpg_sim.sim.updateVariable('V0DtoRGFdiagfh',  0.0)
        cpg_sim.sim.updateVariable('V0DtoRGF',  0.0)
        v, fq, ph,fl_dur,ex_dur, gait = cpg_sim.do_1d_bifurcation('alpha', [0.02, 1.05], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)
    if options.sim == 4:
        '''Simulate ablation of descending long propriospinal neurons '''
        cpg_sim.initialize_simulator()
        cpg_sim.sim.updateVariable('V0VtoRGFdiagfh',  0.0)
        cpg_sim.sim.updateVariable('V0DtoRGFdiagfh',  0.0)
        cpg_sim.sim.updateVariable('inFH',  0.0)
        cpg_sim.sim.updateVariable('V2aHomfh',  0.0)
        v, fq, ph,fl_dur,ex_dur, gait = cpg_sim.do_1d_bifurcation('alpha', [0.02, 1.05], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)

    if options.sim == 5:
        '''Intact model '''
        cpg_sim.initialize_simulator()
        cpg_sim.sim.setAlpha(0.2)
        v, fq, ph,fl_dur,ex_dur, gait = cpg_sim.do_1d_bifurcation('alpha', [0.02, 1.05], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)

    if options.sim == 6:
        cpg_sim.initialize_simulator()
        cpg_sim.sim.setAlpha(0.05)
        dur = 200.0
    

        time_vec = np.arange(0.0,dur,cpg_sim.dt)
        alphas = np.linspace(0.0,-0.2,len(time_vec))
        out = np.zeros((len(time_vec),len(cpg_sim.neurons)))
        for ind_t,alpha in enumerate(alphas):
            cpg_sim.sim.updateVariable('dRGFH',alpha)
            act = cpg_sim.run_step_dense()
            out[ind_t,:]=act

        burst_durs = nsim.simulator.calc_burst_durations(time_vec,out)
        phase_dur,fl_phase_dur,ex_phase_dur,phases,ro = nsim.simulator.calc_phase(time_vec,out,[(0,1)])
        #import IPython;IPython.embed()
        _, ax = plt.subplots(len(neurons)+3, 1, sharex='all')
        for i in range(len(neurons)):
            ax[0].plot(burst_durs[i][:,1],burst_durs[i][:,0])
        ax[0].set_ylim(bottom=0.0)
        ax[1].plot(ro,1.0/phase_dur)
        ax[1].set_ylim(bottom=0.0)
        
        ax[2].plot(ro,phases,'b.',markersize=1.5)
        ax[2].set_ylim([-0.05, 1.05])
        for i in range(len(neurons)):
            ax[3+i].plot(time_vec,out[:,i])
            ax[3+i].set_ylim((0., 1.00))
        plt.show()