//
//  typedefs.h
//  aoi2limb
//
//  Created by Simon Danner on 03/07/2017.
//  Copyright © 2017 Simon Danner. All rights reserved.
//

#ifndef typedefs_h
#define typedefs_h

#include <boost/numeric/ublas/vector.hpp>
#include <array>

struct BodyTilt{
    //anterior/left positive
    double anterior_posterior_angle = 0.0;
    double anterior_posterior_velocity = 0.0;
    double left_right_angle = 0.0;
    double left_right_velocity = 0.0;
};

struct LimbSensorCondition{
    std::vector<double> Ia;
    std::vector<double> Ib;
    std::vector<double> II;
    double cutaneous=0.0;
    LimbSensorCondition(){}; 
    LimbSensorCondition(int N){
        Ia.resize(N);
        Ib.resize(N);
        II.resize(N);
    }
};

typedef boost::numeric::ublas::vector<double,std::vector<double>> myvec;
std::ostream& operator<<(std::ostream &strm, const myvec &re);
std::ostream& operator<<(std::ostream &strm, const std::vector<int> &re);

#endif /* typedefs_h */
