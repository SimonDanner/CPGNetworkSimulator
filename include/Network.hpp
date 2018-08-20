//
//  Network.hpp
//
//  Created by Simon Danner on 22/05/15.
//

#ifndef Network_hpp
#define Network_hpp

#include <stdio.h>
#include <list>
#include <fstream>
#include <string>
#include <map>
#include <chrono>
#include <random>


#include <boost/numeric/odeint.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/circular_buffer.hpp>
#include <ctime>
#include <algorithm>
#include <math.h>

#include "typedefs.h"


using namespace boost::numeric::odeint;

const double c_gBarLeak = 2.8;
const double c_gBarNaP = 5.0;//5.0;

const double c_ELeak = -60.0;
const double c_ENa = 50.0;
const double c_ESynE = -10.0;//-10.0;
const double c_ESynI = -75.0;


const double c_Cmem = 10.0;//20.00;//120.0;//3.0;//20.0

const double c_mV12 = -40;//-40
const double c_mk = -6.5;//-6.5;  //-6.0;

const double c_hV12 = -47; //-55;//-55.0; // -48
const double c_hk = 5.4;  //10.0; // 6
const double c_htau = 110.0; //200;//4000.0;
const double c_hTauK = 50.0; // 50 //-12.0;
const double c_hTauV12 = c_hV12;
const double c_hTau0 = 50.0;

const double c_Vmin = -50.0;//-50.0;
const double c_Vmax = 0.0;

const double c_tauOutput = 5.;

const double c_sigmaNoise = 0.005;

const double c_tauNoise = 10.0;

//const double settingPeriod = 10000.0; //ms

struct connection{
    int from;
    double *weight; //offset
    double *slope; //ms^-1
};
typedef boost::numeric::ublas::vector< std::list<connection> > connection_matrix;


class OrnsteinUhlenbeck{
private:
    int Nhist=128;
    unsigned seed;
    double mu=0.0;
    double sigma=1.0;
    double tau=1.0;
    int N=1;
    boost::circular_buffer<myvec> data;
    double time0=0.0;
    std::mt19937 generator;
    std::normal_distribution<double> distribution;
public:
    OrnsteinUhlenbeck(int n, double m, double s, double t);
    OrnsteinUhlenbeck(){};
    void calculateUntil(double t);
    myvec get(double t);
};

// sturcture describing a drive-connection
struct drive{
    int to;         // index of neuron
    double *weight; // pointer to connection weight
    double *offset; // pointer to offset
    drive(int toc, double *weightc, double *offsetc){
        to=toc;
        weight=weightc;
        offset=offsetc;
    }
};

struct feedback{
    int fromleg;
    int frommg;
    int to;
    double* weight;
    feedback(int to_,int fromleg_,int frommg_,double* weight_):to(to_),frommg(frommg_),fromleg(fromleg_),weight(weight_){}
};

struct feedback_body_tilt{
    int direction;
    double* cutoff;
    double* weight;
    int to;
    int type;
};

class Network{
private:
public:
    BodyTilt body_tilt;
    std::vector<LimbSensorCondition> lscond;
    
    std::map<int,std::string> names;
    myvec initial;
    myvec gBarLeak;
    myvec gBarNaP;
    
    myvec ELeak;
    myvec ENa;
    myvec ESynE;
    myvec ESynI;
    
    myvec Cmem;
    
    myvec mk;
    myvec mV12;
    
    myvec hk;
    myvec hV12;
    myvec htau;
    myvec hTauK;
    myvec hTauV12;
    myvec hTau0;
    
    myvec Vmax;
    myvec Vmin;
    
    myvec tauOutput;
    
    myvec sigmaNoise;
    myvec tauNoise;
    
    std::vector<int> outputFunction;
    
    int N_NaP;
    int N_norm;
    
    int iNaP=0;
    int iNorm=0;
    
    int simDirection=1;
    
    connection_matrix connE;
    connection_matrix connI;
    connection_matrix connC;
    
    std::list<drive> driveE;
    std::list<drive> driveI;
    
    std::list<feedback> feedbackIa;
    std::list<feedback> feedbackIb;
    std::list<feedback> feedbackII;
    std::list<feedback> feedbackCutaneous;
    std::list<feedback_body_tilt> feedbackBodyTilt;
    
    myvec transV;
    
    OrnsteinUhlenbeck randomGen;
    
    int N_Neurons;
    
    double simDuration = 100000;
    double settingPeriod = 10000;
    double scalingFactor = 1.0;
    
    bool stepwise = false;
    int nSteps = 50;
    
    double alpha = 0.0;
    
    double sf = 1.0;

    std::map<std::string,double*> variableMap;
    
    std::vector<std::vector<int>> in_Act;
    
    Network();
    Network(int N_NP, int N_normal);
    Network(std::string /*filename*/,std::vector<std::string> /*musclenames*/, std::vector<std::vector<std::string>> /*mnnames*/);
    
    void initialize(int /*N_NP*/, int /*N_normal*/);
    
    static myvec create_scalarV(int N, double k);
    static void set_para(myvec &to,double value,int start,int end);
    static void assign_para(myvec &to,myvec from);
    myvec genInitialCond() const;
    
    void setConnE(const int from,const int to, double *w);
    void setConnI(const int from,const int to, double *w);
    void setConnC(const int from,const int to,double *w);
    void setConnE(const int from,const int to, double *w, double *slope);
    void setConnI(const int from,const int to, double *w, double *slope);
    void resetConn(connection_matrix&, const int, const int,  double *w);
    void resetConnE(const int from,const int to, double *w){resetConn(connE,from,to,w);};
    void resetConnI(const int from,const int to, double *w){resetConn(connI,from,to,w);};
    
    double* getConn(connection_matrix& , const int , const int );
    double* getConnE(const int, const int);
    double* getConnI(const int, const int);
    
    
    void setDriveE(const int to,  double *weight,  double *offset);
    void setDriveI(const int to,  double *weight,  double *offset);

    void setFeedbackIa(const int to, const int fromleg, const int frommg,  double* weight);
    void setFeedbackIb(const int to, const int fromleg, const int frommg,  double* weight);
    void setFeedbackII(const int to, const int fromleg, const int frommg,  double* weight);
    void setFeedbackCutaneous(const int to, const int fromleg, double* weight);
    
    std::string getName(const int /*neuronID*/) const;
    int getIndex(const std::string /*name*/) const;

    
    friend std::ostream& operator<<(std::ostream&,const Network&);
    
    //integration step
    void step(const myvec &, myvec &, double);
    
    void run(double);
    bool updateVariable(const std::string,const double);
    double getVariableValue(std::string);
    double* getVariablePointers(std::string);
    void setLscond(std::vector<LimbSensorCondition>& ls_){
        lscond=ls_;
    }
    void setPara(std::string parameter, std::string neuron_substr,double value);
};

bool allEqual(myvec);

#endif /* defined(__iCPG__Network__) */
