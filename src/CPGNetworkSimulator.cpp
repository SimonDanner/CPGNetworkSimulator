//
//  Solver.cpp
//  aoi2limb
//
//  Created by Simon Danner on 29/06/2017.
//  Copyright © 2017 Simon Danner. All rights reserved.
//

#include "CPGNetworkSimulator.hpp"

using namespace boost::numeric::ublas;


CPGNetworkSimulator::CPGNetworkSimulator(const std::string filename,const std::vector<std::vector<std::string>> mnnames_){
    mnnames=mnnames_;
    net=new Network(filename,mnnames);
    sys.net = net;
    initialize();
}

void CPGNetworkSimulator::initialize(){
    state=net->genInitialCond();
    integrate_const( stepper, sys, state , 0.0 , 10., 0.002);
    for(int i = 0;i<net->in_Act.size();++i){
        std::vector<double> v_;
        for(int j = 0;j<net->in_Act[0].size();++j){
            v_.push_back(net->transV[net->in_Act[i][j]]);
        }
        act.push_back(v_);
    }
}


void CPGNetworkSimulator::step(double dt_){
    dt=dt_;
    stepper.do_step(sys,state,t0,dt);
    for(int i = 0;i<net->in_Act.size();++i){
        for(int j = 0;j<net->in_Act[0].size();++j){
            act[i][j] = net->transV[net->in_Act[i][j]];
        }
    }
    
    t0+=dt;
}

bool CPGNetworkSimulator::updateVariable(const std::string var, double value){
    return net->updateVariable(var,value);
}