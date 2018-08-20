//
//  CPGNetworkSimulator.hpp
//
//  Created by Simon Danner on 29/06/2017.
//

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
};

#endif /* Solver_hpp */
