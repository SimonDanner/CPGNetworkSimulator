/*  CPGNetworkSimulator.hpp
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

#ifndef CPGNetworkSimulator_hpp
#define CPGNetworkSimulator_hpp

#include <stdio.h>
#include <boost/numeric/odeint/iterator/times_time_iterator.hpp>
#include "Network.hpp"
#include "typedefs.h"


struct OdeSystemNetwork
{
    Network *net;
    void operator()( const myvec &x , myvec &dxdt , double t){
        net->step(x,dxdt,t*1000.0);
        dxdt*=1000.0;
    }
};


typedef runge_kutta_cash_karp54< myvec > stepper_type;

struct UpdateList{
    std::string name;
    std::vector<std::string> vars;
    std::vector<double> initialValues;
    std::vector<double*> pointers;
    double scale = 1.0;
};

class CPGNetworkSimulator{
private:
    stepper_type stepper;
    controlled_runge_kutta< stepper_type> controlled_stepper = make_controlled( 1.0e-6 , 1.0e-6 , stepper_type() );;
    int N_last_update=0;
    double t0=0.0;
    OdeSystemNetwork sys;
    std::vector< std::vector<double> > act;
    std::vector< std::vector<std::string> > mnnames;
    std::vector<double*> variableVectorPointers;
    std::vector<std::string> variableVectorNames;
    double dt;
    myvec state;
    Network* net;
    void initialize();
public:
    CPGNetworkSimulator(const std::string filename,const std::vector<std::string> musclenames,const std::vector<std::vector<std::string>> mnnames_);
    void setAlpha(double alpha){net->alpha = alpha;};
    void step(double dt);
    void controlled_step(double dt);
    const std::vector<std::vector<double>>& getAct(){return act;};
    bool updateVariable(const std::string var, double value);
    void setLscond(std::vector<LimbSensorCondition>& ls_){
        net->setLscond(ls_);
    }
    void setBodyTilt(double apa, double apv,double lra, double lrv){
        net->body_tilt.anterior_posterior_angle=apa;
        net->body_tilt.anterior_posterior_velocity=apv;
        net->body_tilt.left_right_angle=lra;
        net->body_tilt.left_right_velocity=lrv;
    }
    std::vector<double> setupVariableVector(const std::vector<std::string> variablenames);
    void updateVariableVector(const std::vector<double> values);
    std::vector<double> getState(){return state.data();};
    void setState(std::vector<double> s){state = myvec(s);};
    void updateParameter(std::string name, double value);
};

#endif /* Solver_hpp */
