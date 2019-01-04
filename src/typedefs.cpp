/*  typedefs.cpp
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
