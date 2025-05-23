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


CPGNetworkSimulator::CPGNetworkSimulator(const std::string filename,const std::vector<std::string> musclenames, const std::vector<std::vector<std::string>> mnnames_,bool debug){
    mnnames=mnnames_;
    net=new Network(filename,musclenames,mnnames,debug);
    sys.net = net;
    initialize();
}

void CPGNetworkSimulator::initialize(){
    state=net->genInitialCond();
    dense_stepper.initialize( state , 0.0 , 0.001 );
    integrate_const( stepper, sys, state , 0.0 , 10., 0.002);
    for(int i = 0;i<net->in_Act.size();++i){
        std::vector<double> v_;
        std::vector<double> ie_;
        std::vector<double> ii_;
        for(int j = 0;j<net->in_Act[0].size();++j){
            v_.push_back(net->transV[net->in_Act[i][j]]);
            ie_.push_back(net->Iepsp[net->in_Act[i][j]]);
            ii_.push_back(net->Iipsp[net->in_Act[i][j]]);
        }
        act.push_back(v_);
        Iepsp.push_back(ie_);
        Iipsp.push_back(ii_);
    }
}


void CPGNetworkSimulator::step(double dt_){
    dt=dt_;
    stepper.do_step(sys,state,t0,dt);
    for(int i = 0;i<net->in_Act.size();++i){
        for(int j = 0;j<net->in_Act[0].size();++j){
            act[i][j] = net->transV[net->in_Act[i][j]];
            Iepsp[i][j] = net->Iepsp[net->in_Act[i][j]];
            Iipsp[i][j] = net->Iipsp[net->in_Act[i][j]];
        }
    }
    
    t0+=dt;
}

void CPGNetworkSimulator::step(double dt_,double error){
    dt=dt_;
    myvec xerr = myvec(state.size());
    double maxerr = 0.0;
    for (int i = 0;i<N_substeps;++i){
        double dt_sub = dt/double(N_substeps);
        stepper.do_step(sys,state,t0,dt_sub,xerr);
        for(int i = 0;i<xerr.size();++i){
            maxerr = std::max<double>(maxerr,std::abs(xerr[i]));
        }
        t0+=dt_sub;
    }
    if(maxerr>error){
        N_substeps*=2;
    }else if(maxerr<error/2.0){
        if(N_substeps>1){
            N_substeps/=2;
        }
    }
    if (N_substeps>32){
        N_substeps=32;
    }

    for(int i = 0;i<net->in_Act.size();++i){
        for(int j = 0;j<net->in_Act[0].size();++j){
            act[i][j] = net->transV[net->in_Act[i][j]];
            Iepsp[i][j] = net->Iepsp[net->in_Act[i][j]];
            Iipsp[i][j] = net->Iipsp[net->in_Act[i][j]];
        }
    }
    
}

void CPGNetworkSimulator::controlled_step(double dt_){
    dt=dt_;
    integrate_const( controlled_stepper , sys , state , t0 , t0+dt , dt );
    //controlled_stepper.try_step(sys,state,t0,dt);
    for(int i = 0;i<net->in_Act.size();++i){
        for(int j = 0;j<net->in_Act[0].size();++j){
            act[i][j] = net->transV[net->in_Act[i][j]];
            Iepsp[i][j] = net->Iepsp[net->in_Act[i][j]];
            Iipsp[i][j] = net->Iipsp[net->in_Act[i][j]];
        }
    }
    
    t0+=dt;
}
void CPGNetworkSimulator::dense_step(double dt_){
    dt=dt_;
    while ( t_last_dense < t0+dt){
        auto t = dense_stepper.do_step(sys);
        t_last_dense = t.second;
    }
    dense_stepper.calc_state( t0+dt , state );

    for(int i = 0;i<net->in_Act.size();++i){
        for(int j = 0;j<net->in_Act[0].size();++j){
            int index = net->in_Act[i][j];
            act[i][j] = net->transV[net->in_Act[i][j]];
            Iepsp[i][j] = net->Iepsp[net->in_Act[i][j]];
            Iipsp[i][j] = net->Iipsp[net->in_Act[i][j]];
        }
    }
    
    t0+=dt;
}

bool CPGNetworkSimulator::updateVariable(const std::string var, double value){
    return net->updateVariable(var,value);
}

double CPGNetworkSimulator::getVariableValue(const std::string var){
    return net->getVariableValue(var);
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