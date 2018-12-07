//
// Created by darius on 30/10/18.
//

#include <iostream>
#include "MeasurementData.h"

std::istream& operator>>(std::istream& input, FunctionData& d){
    input >> d.x >> d.y;
    return input;
}
