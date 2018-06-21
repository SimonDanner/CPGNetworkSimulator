//
//  Solver.hpp
//  aoi2limb
//
//  Created by Simon Danner on 29/06/2017.
//  Copyright © 2017 Simon Danner. All rights reserved.
//

#ifndef Solver_hpp
#define Solver_hpp

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
    int N_last_update=0;
    double t0=0.0;
    OdeSystemNetwork sys;
    std::vector<std::vector<double>> act;
    std::vector<std::vector<std::string>> mnnames;
    double dt;
    myvec state;
    Network* net;
    void initialize();

public:
    CPGNetworkSimulator(const std::string filename,const std::vector<std::vector<std::string>> mnnames_);
    void setAlpha(double alpha){net->alpha = alpha;};
    void step(double dt);
    std::vector<std::vector<double>>& getAct(){return act;};
    bool updateVariable(const std::string var, double value);
    void setLscond(std::vector<LimbSensorCondition>& ls_){
        net->setLscond(ls_);
        std::cout << ls_[0].Ia[0] << std::endl;
    }
};

#endif /* Solver_hpp */
