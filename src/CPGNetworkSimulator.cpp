/*  CPGNetworkSimulator.cpp: C++ interface for simulator
    Copyright (C) 2019  Simon Danner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "CPGNetworkSimulator.hpp"

using namespace boost::numeric::ublas;


CPGNetworkSimulator::CPGNetworkSimulator(const std::string filename,const std::vector<std::string> musclenames, const std::vector<std::vector<std::string>> mnnames_){
    mnnames=mnnames_;
    net=new Network(filename,musclenames,mnnames);
    sys.net = net;
    initialize();
}

void CPGNetworkSimulator::initialize(){
    state=net->genInitialCond();
    integrate_const( controlled_stepper, sys, state , 0.0 , 10., 0.002);
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
void CPGNetworkSimulator::controlled_step(double dt_){
    dt=dt_;
    integrate_const( controlled_stepper , sys , state , t0 , t0+dt , dt );
    //controlled_stepper.try_step(sys,state,t0,dt);
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

std::vector<double> CPGNetworkSimulator::setupVariableVector(const std::vector<std::string> variablenames){
    variableVectorNames = variablenames;
    variableVectorPointers = std::vector<double*>(variableVectorNames.size());
    std::vector<double> values(variableVectorNames.size());
    for(int i = 0;i < variableVectorNames.size();++i){
        variableVectorPointers[i]=net->getVariablePointers(variableVectorNames[i]);
        values[i] = *variableVectorPointers[i];
    }
    return values;
}

void CPGNetworkSimulator::updateVariableVector(const std::vector<double> values){
    if(variableVectorPointers.size()!=values.size()){
        std::cout << "CPGNetworkSimulator::updateVariableVector: values has wrong size" << std::endl;
        return;
    }
    for(int i = 0;i < variableVectorPointers.size();++i){
        (*variableVectorPointers[i])=values[i];
    }
};

void CPGNetworkSimulator::updateParameter(std::string name, double value){
    if(name=="sigmaNoise"){
        std::fill(net->sigmaNoise.begin(),net->sigmaNoise.end(),value);
    }
};