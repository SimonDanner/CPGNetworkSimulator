import CPGNetworkSimulator as nsim
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib as npml
from scipy.stats import circstd, circmean
import os
from scoop import futures

import functools

class simulator:
    stdp_limit = 0.005
    its_limit = 20
    initialized = False
    def __init__(self,**kwargs):
        self.neurons = kwargs.get('neurons', ["RGF_NaP_L_hind", "RGF_NaP_R_hind", "RGF_NaP_L_front", "RGF_NaP_R_front"])
        self.filename = kwargs.get('filename', "./models/MLR_45.txt")
        self.dt = kwargs.get('dt',0.001)
        self.duration = kwargs.get('duration',10.0)
        self.phase_diffs = kwargs.get('phase_diffs',[(0,1),(2,3),(0,2),(0,3),(1,3),(1,2)])
        self.phase_diff_names = kwargs.get('phase_diff_names',['L-R hind','L-R fore','homolateral','diagonal','homolateral','diagonal'])
        self.alphainit = 0.0
        self.doPreRun = True

    def initialize_simulator(self):
        if not self.initialized:
            fn = os.path.join(os.path.dirname(__file__),"..",self.filename)
            self.sim = nsim.CPGNetworkSimulator(fn,["a","b"],(self.neurons,))
            self.sim.setAlpha(self.alphainit)
            self.setDuration(self.duration)
            self.initialized=True
            if self.doPreRun:
                for i in range(10):
                    self.run_sim()
                self.IC = self.sim.getState()
            
    def updateVariable(self,name,value):
        if(name=="alpha"):
            self.sim.setAlpha(value)
        else:
            self.sim.updateVariable(name,value)

    def setDuration(self,dur):
        self.duration = dur
        self.time_vec = np.arange(0.0,self.duration,self.dt)

    def run_sim(self):
        out = np.zeros((len(self.time_vec),len(self.neurons)))
        for ind_t,t in enumerate(self.time_vec):
            self.sim.step(self.dt)
            act = self.sim.getAct()
            out[ind_t,:]=act[0]
        return out

    @staticmethod
    def calc_phase(time_vec,out,phase_diffs):
        os_=((np.diff((out>0.1).astype(np.int),axis=0)==1).T)
        of_=((np.diff((out>0.1).astype(np.int),axis=0)==-1).T)
        onsets=npml.repmat(time_vec[:-1],out.shape[1],1)[os_]
        offsets=npml.repmat(time_vec[:-1],out.shape[1],1)[of_]
        leg_os=(npml.repmat(np.arange(out.shape[1]),len(time_vec)-1,1).T)[os_]
        leg_of=(npml.repmat(np.arange(out.shape[1]),len(time_vec)-1,1).T)[of_]
        
        times_os=np.stack((onsets,leg_os,np.arange(len(leg_os))),1)
        times_os=times_os[times_os[:,0].argsort()]
        times_of=np.stack((offsets,leg_of,np.arange(len(leg_of))),1)
        times_of=times_of[times_of[:,0].argsort()]

        times = np.concatenate((
                    np.concatenate((times_os,np.ones((len(times_os),1))*0.0),1),
                    np.concatenate((times_of,np.ones((len(times_of),1))*1.0),1)))
        times=times[times[:,0].argsort()]
        
        ref_onsets = times[np.logical_and(times[:,1]==0,times[:,3]==0)][:,0]
        phase_dur=np.append(ref_onsets[1:]-ref_onsets[:-1],np.nan)

        p = times[times[:,1]==0]
        indices = np.where(np.diff(p[:,3])==1)
        fl_phase_dur = np.zeros((len(ref_onsets)))
        fl_phase_dur[:] = np.nan
        fl_phase_dur[p[indices,2].astype(int)] = p[[ind+1 for ind in indices],0] - p[indices,0]
        ex_phase_dur = phase_dur-fl_phase_dur
        
        M = np.zeros((len(ref_onsets),out.shape[1]))
        M[:] = np.nan
        M[:,0]=ref_onsets
        

        for i in range(1,out.shape[1]):
            p = times[np.logical_and((times[:,1]==0) | (times[:,1]==i),times[:,3]==0)]
            indices = np.where(np.diff(p[:,1])==i)
            M[p[indices,2].astype(int),i] = p[[ind+1 for ind in indices],0]

        phases=np.zeros((len(ref_onsets),len(phase_diffs)))
        for i,(x,y) in enumerate(phase_diffs):
            phases[:,i] = ((M[:,y]-M[:,x])/phase_dur)  % 1.0

        if phases.shape[0]!=0:
            no_nan = ~np.isnan(np.concatenate(
                        (np.stack((phase_dur,fl_phase_dur,ex_phase_dur),1),phases),1
                        )).any(axis=1)
            return (phase_dur[no_nan],fl_phase_dur[no_nan],ex_phase_dur[no_nan],phases[no_nan])
        else:
            return (phase_dur,fl_phase_dur,ex_phase_dur,phases)  

    @staticmethod
    def classify_gait_simple(duty_factor,phases): 
        
        lr_hind = phases[:,0]
        homolateral = phases[:,2]
        diagonal = phases[:,3]

        lr_alt = np.logical_and(lr_hind>=0.25,lr_hind<=0.75)
        hl_alt = np.logical_and(homolateral >= 0.25,homolateral <= 0.75)
        di_alt = np.logical_and(diagonal >= 0.25,diagonal <= 0.75)
        diag_quarter = np.logical_or(
                        np.logical_and(diagonal>0.1,diagonal<=0.4),
                        np.logical_and(diagonal>=0.6,diagonal<0.9))
        walk = np.logical_and(np.logical_and(duty_factor>0.5,lr_alt),
                    diag_quarter)
        trot = np.logical_and(lr_alt, np.logical_and(
                    np.logical_or(diagonal <= 0.1,diagonal >= 0.9), 
                    hl_alt))
        gallop = np.logical_and(hl_alt,
                    np.logical_or(
                        np.logical_and(lr_hind>0.025,lr_hind<=0.25),
                        np.logical_and(lr_hind>=0.75,lr_hind<0.975)))   
        bound = np.logical_and(hl_alt,
                    np.logical_or(lr_hind<=0.025,lr_hind>=0.975))
        if len(duty_factor) > 0:
            gaits = np.zeros((len(duty_factor),),dtype=int)
            gaits[walk]=1
            gaits[trot]=2
            gaits[gallop]=3
            gaits[bound]=4
        else:
            gaits = [np.nan]
        return gaits      

    def do_iteration(self):
        mfq = np.nan
        mphases_ = [] 
        max_std_phases = 1.0
        its=0
        while max_std_phases > self.stdp_limit and its < self.its_limit:
            out = self.run_sim()
            phase_dur,fl_phase_dur,ex_phase_dur, phases = self.calc_phase(self.time_vec,out,self.phase_diffs)
            if len(phase_dur)>10:
                mphases_ = circmean(phases[-5:,:],1.0,0.0,0)
                max_std_phases = np.max(circstd(phases[-5:,:],1.0,0.0,0))
                mfq = 1.0/np.nanmean(phase_dur[-5:])
            its+=1
        gaits = self.classify_gait_simple((ex_phase_dur/phase_dur)[-5:],phases[-5:,:])
        return (mfq, mphases_, gaits[-1:])

    def do_1d_bifurcation(self,variable_name,range_,steps,updown=True):
        self.initialize_simulator()
        v=range_[0]+(np.arange(0,steps,1))/(steps-1.0)*(range_[1]-range_[0])
        
        frequency = np.zeros((steps,2))*np.nan
        gait = np.zeros((steps,2),dtype=int)*np.nan
        phases = np.zeros((steps,len(self.phase_diffs),2))*np.nan
        IChist=list()
        j_start_back=0
        go_up_on_nan=True
        #self.sim.setState(IC)

        self.updateVariable(variable_name,v[0])
        for j in range(0,steps):
            IChist.append(self.sim.getState())
            self.updateVariable(variable_name,v[j])

            fq, phases_, gait_ = self.do_iteration()
            if not np.isnan(fq):
                frequency[j,0]=fq
                gait[j,0]=gait_
                phases[j,:,0]=phases_
            j_start_back=j-1
            if np.isnan(fq) and not go_up_on_nan:
                j_start_back=j-2
                break
            if not np.isnan(fq):
                go_up_on_nan=False
                
        
        if updown:
            self.sim.setState(IChist[j_start_back])
            for j in np.arange(j_start_back,-1,-1):
                self.updateVariable(variable_name,v[j])
                fq, phases_, gait_ = self.do_iteration()
                frequency[j,1]=fq
                gait[j,1]=gait_
                phases[j,:,1]=phases_
                if np.isnan(fq):
                    break
        return (v,frequency,phases,gait)

    def do_1d_bifurcation_helper(self,at_value,bi_variable_name,bi_range,bi_steps,at_variable_name,updown=True):
        self.initialize_simulator()
        self.sim.setState(self.IC)
        self.its_limit=5
        self.updateVariable(at_variable_name,at_value)
        return self.do_1d_bifurcation(bi_variable_name,bi_range,bi_steps,updown)

    def do_2d_bifurcation(self,variable_names,ranges,steps,updown=False):
        v0=ranges[0][0]+(np.arange(0,steps[0],1))/(steps[0]-1.0)*(ranges[0][1]-ranges[0][0])
        frequency = np.zeros(steps+(2,))*np.nan
        gaits = np.zeros(steps+(2,))*np.nan
        phases = np.zeros(steps+(len(self.phase_diffs),2))*np.nan
        paras = futures.map(
                        functools.partial(
                                self.do_1d_bifurcation_helper,
                                bi_variable_name=variable_names[1],
                                bi_range=ranges[1],
                                bi_steps=steps[1],
                                at_variable_name=variable_names[0],
                                updown=updown),
                        v0)
        v1 = []
        for i,(v_,fq,ph,g) in enumerate(paras):
            frequency[i,:,:]=fq
            phases[i,:,:,:]=ph
            gaits[i,:,:]=g
            v1=v_
        return (v0,v1,frequency,phases,gaits)

    def do_noise_iteration(self,value, variable_name,sigma,set_variables=list()):
        self.initialize_simulator()
        self.sim.updateParameter('sigmaNoise',sigma)
        self.updateVariable(variable_name,value)
        for (name,value) in set_variables:
            self.updateVariable(name,value)
        
        out = self.run_sim()
        phase_dur,fl_phase_dur,ex_phase_dur, phases = self.calc_phase(self.time_vec,out,self.phase_diffs)
        gaits = self.classify_gait_simple((ex_phase_dur/phase_dur),phases)
        return (1/phase_dur,fl_phase_dur,ex_phase_dur, phases, gaits)

    def do_noise_simulation(self,variable_name,range_,steps,sigma,set_variables=list()):
        v=range_[0]+(np.arange(0,steps,1))/(steps-1.0)*(range_[1]-range_[0])
        self.duration = 100.0
        self.doPreRun = False
        fqs = list()
        fl_phase_durs = list()
        ex_phase_durs = list()
        phases = list()
        gaits = list()
        paras = futures.map(
                        functools.partial(self.do_noise_iteration,
                            variable_name=variable_name,
                            sigma=sigma,
                            set_variables=set_variables),
                        v)
        for i,(fq,fl,ex,ph,g) in enumerate(paras):
            fqs.append(fq)
            fl_phase_durs.append(fl)
            ex_phase_durs.append(ex)
            phases.append(ph)
            gaits.append(g)
        
        return (v,fqs,fl_phase_durs,ex_phase_durs,phases,gaits)
        