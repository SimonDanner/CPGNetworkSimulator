# CPGNetworkSimulator
CPGNetworkSimulator is a simulator for neural network models using activity-dependent neuron population models. The simulator is written in C++-11 and uses boost's odeint for efficient numeric integration of the system of differential equations. A python interface is provided to control the simulation, calculate bifurcation diagrams, integrate the simulation within other simulators, and analsye the output. 

This simulator has been developed for the simulations performed the following publications:
- Danner SM, Shevtsova NA, Frigon A, Rybak IA. Computational modeling of spinal circuits controlling limb coordination and gaits in quadrupeds. eLife. 2017 Nov 22;6. pii: e31050. doi: 10.7554/eLife.31050. 
- Danner SM, Wilshin SD, Shevtsova NA, Rybak IA. Central control of interlimb coordination and speed-dependent gait expression in quadrupeds. J Physiol. 2016 Dec 1;594(23):6947-6967. doi: 10.1113/JP272787. 
- Ausborn J, Shevtsova NA, Caggiano V, Danner SM, Rybak IA. Computational modeling of brainstem circuits controlling locomotor frequency and gait. eLife. in-revision.

## Installation
To install CPGNetworkSimulator locally run
```
root_dir> python setup.py build
root_dir> python -m pip install .
```
### Requirements: 
- python 3.5 or newer
- C++ toolchaing
- boost 1.65 or newer (only headers for odeint are required)
- yaml-cpp
- pybind11
- numpy
- scipy
- matplotlib

## Usage
mlr_simulation.py contains scripts to run all simulations of Ausborn et al.
```
python3 [-m scoop] mlr_simulations.py -s [index of simulation]
```
For more detailed instructions refer to the comments in mlr_simulations.py.
