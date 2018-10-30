//
// Created by darius on 30/10/18.
//

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include "MeasurementData.h"



std::ostream& operator<<(std::ostream& os, const MeasurementData& md){
    os << md.depth << " m \t" << md.grav << " mGal" << std::endl;
    return os;
}

std::istream& operator>>(std::istream& input, MeasurementData& md){
    input >> md.depth >> md.grav;
    // apply free air gradient correction by subtracting it from the data
    md.grav -= FREE_AIR_GRADIENT * md.depth;
    return input;
}
