//
// Created by darius on 11/11/18.
//

#include <iostream>
#include "Result.h"

std::ostream& operator<<(std::ostream& os, const Result& md){
    os << md.depth << " " << md.density << "\n";
    return os;
}