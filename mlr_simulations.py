"""
mlr_simulations.py

Script to perform all simulations illustrated in Ausborn J, Shevtsova NA,
Caggiano V, Danner SM, Rybak IA. "Computational modeling of brainstem 
circuits controlling locomotor frequency and gait" submitted to eLife.

Usage:
python3 [-m scoop] mlr_simulations.py -s [index of simulation]

e.g. "python3 mlr_simulations.py -s 0" calculates 1D bifurcation diagrams for CNF stimulation

For 1D bifurcation diagrams use index 0-3 (Figure 3, 4A)
For 2D bifurcation diagrams use index 10-13 (Figure 6)
For gait-frequency / phase-frequency plots under application of noise use index 20-23
(For more details see comments below)

2D bifurcation diagrams and noise-simulations support parallelization through scoop.
"""
import CPGNetworkSimulator.tools.py_simulator as nsim
from CPGNetworkSimulator.tools.plt import plot_1d_bifurcation_diagram, plot_2d_bifurcation_diagram, plot_noise_phase_fq, plot_noise_gait_fq
from optparse import OptionParser
import matplotlib.pyplot as plt

if __name__ == "__main__":
    neurons = ["RGF_NaP_L_hind", "RGF_NaP_R_hind",      # neurons to be read every time step 
               "RGF_NaP_L_front", "RGF_NaP_R_front"]
    filename = "./models/Ausborn-etal-eLife.txt" #  network model configuration file 

    #Note: the resolution of all calculations has been reduced to speed up simulation time.
    #Change following three values to increase number of steps
    steps_1D = 50       #number of steps for 1D-bifurcation diagrams
    steps_2D = (20,20)  #number of steps for 2D-bifurcation diagrams
    steps_noise = 106   #number of steps for noise simulations

    parser = OptionParser()
    parser.add_option("-s", "--simulation", dest="sim",default=0,type=int)
    options, args = parser.parse_args()

    cpg_sim = nsim.simulator(neurons=neurons, filename=filename) # instantiate simulator

    if options.sim == 0:
        """Bifurcation diagram unilateral stimulation of glutaminergic neurons in CNF"""
        v, fq, ph, gait = cpg_sim.do_1d_bifurcation(
            'd0_CnF_Glu_L', [3.95, 5.3675], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)
        plt.show()

    elif options.sim == 1:
        """Bifurcation diagram unilateral stimulation of glutaminergic neurons in PPN"""
        v, fq, ph, gait = cpg_sim.do_1d_bifurcation(
            'd0_PPN_Glu', [4.0, 5.5], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)
        plt.show()

    elif options.sim == 2:
        """Bifurcation diagram unilateral stimulation of glutaminergic neurons in CNF
            while PPN is inactivated"""
        cpg_sim.initialize_simulator()
        cpg_sim.sim.updateVariable('PPN_Glu_to_LGPi_Glu1',  0.0)
        cpg_sim.sim.updateVariable('PPN_Glu_to_cLGPi_Glu1', 0.0)
        cpg_sim.sim.updateVariable('PPN_GAT_to_PPN_Glu',    0.0)

        v, fq, ph, gait = cpg_sim.do_1d_bifurcation(
            'd0_CnF_Glu_L', [4.2, 6.8775], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)
        plt.show()

    elif options.sim == 3:
        """Bifurcation diagram unlilateral stimulation of LPGi"""
        v, fq, ph, gait = cpg_sim.do_1d_bifurcation(
            'd0_LPGi_Glu', [2.45, 2.45+1.1], steps_1D)
        plot_1d_bifurcation_diagram(v, fq, ph, gait)
        plt.show()

    elif options.sim == 10:
        """2D bifurcation diagram stimulation of unilateral CnF GAT vs bilateral CnF Glu"""
        v0, v1, fq, ph, gaits = cpg_sim.do_2d_bifurcation(
            ('d0_CnF_GAT', 'd0_CnF_Glu_bl'), ([1.1, 3.0], [2.78, 3.1]), steps_2D)
        plot_2d_bifurcation_diagram(v0, v1, fq, ph)
        plt.show()
        

    elif options.sim == 11:
        """2D bifurcation diagram stimulation of unilateral PPN GAT vs bilateral CnF Glu"""
        v0, v1, fq, ph, gaits = cpg_sim.do_2d_bifurcation(
            ('d0_PPN_GAT', 'd0_CnF_Glu_bl'), ([1.1, 3.0], [2.78, 3.1]), steps_2D)
        plot_2d_bifurcation_diagram(v0, v1, fq, ph)
        plt.show()
        
    elif options.sim == 12:
        """2D bifurcation diagram stimulation of unilateral LPGi GAT vs bilateral CnF Glu"""
        v0, v1, fq, ph, gaits = cpg_sim.do_2d_bifurcation(
            ('d0_LPGi_GAT', 'd0_CnF_Glu_bl'), ([1.1, 3.0], [2.78, 3.1]), steps_2D)
        plot_2d_bifurcation_diagram(v0, v1, fq, ph)
        plt.show()

    elif options.sim == 20:
        """Calculate frequency vs gait/phase difference plots """
        """Unilateral stimulation of CnF_Glu neurons """
        v, fqs, fl_phase_durs, ex_phase_durs, phases, gaits = cpg_sim.do_noise_simulation(
            'd0_CnF_Glu_L', [3.95, 5.3675], steps_noise, 1.0)
        plot_noise_gait_fq(fqs, gaits)
        plot_noise_phase_fq(fqs, phases)
        plt.show()
    elif options.sim == 21:
        """Unilateral stimulation of PPN_Glu neurons """
        v, fqs, fl_phase_durs, ex_phase_durs, phases, gaits = cpg_sim.do_noise_simulation(
            'd0_PPN_Glu', [4.0, 5.5], steps_noise, 1.0)
        plot_noise_gait_fq(fqs, gaits)
        plot_noise_phase_fq(fqs, phases)
        plt.show()
    elif options.sim == 22:
        """Unilateral stimulation of CnF_Glu neurons while PPN is inactivated"""
        sv = [('PPN_Glu_to_LGPi_Glu1',  0.0),
              ('PPN_Glu_to_cLGPi_Glu1',  0.0), ('PPN_GAT_to_PPN_Glu',  0.0)]
        v, fqs, fl_phase_durs, ex_phase_durs, phases, gaits = cpg_sim.do_noise_simulation(
            'd0_CnF_Glu_L', [4.2, 6.8775], steps_noise, 1.0, sv)
        plot_noise_gait_fq(fqs, gaits)
        plot_noise_phase_fq(fqs, phases)
        plt.show()
    elif options.sim == 23:
        """Unilateral stimulation of LPGi_Glu neurons """
        v, fqs, fl_phase_durs, ex_phase_durs, phases, gaits = cpg_sim.do_noise_simulation(
            'd0_LPGi_Glu', [2.45, 2.45+1.1], steps_noise, 1.0)
        plot_noise_gait_fq(fqs, gaits)
        plot_noise_phase_fq(fqs, phases)
        plt.show()