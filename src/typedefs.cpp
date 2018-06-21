//
//  typedefs.cpp
//  aoi2limb
//
//  Created by Simon Danner on 03/07/2017.
//  Copyright © 2017 Simon Danner. All rights reserved.
//

#include <stdio.h>
#include "typedefs.h"

std::ostream& operator<<(std::ostream &strm, const myvec &re){
    auto prec= strm.precision();
    strm.precision(5);
    for(int i=0;i<re.size();++i){
        strm << re[i] << "\t";
    }
    strm.precision(prec);
    return strm;
}

std::ostream& operator<<(std::ostream &strm, const std::vector<int> &re){
    for(int i=0;i<re.size();++i){
        strm << re[i] << "\t";
    }
    return strm;
}
