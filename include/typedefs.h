/*  typedefs.h
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
