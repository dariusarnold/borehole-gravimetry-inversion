//
// Small utility functions or classes
//

#ifndef PROGRAMM_UTILS_H
#define PROGRAMM_UTILS_H

#include <iostream>


/*
 * print a vector of printable elements to ostream, element for element
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> vec){
    for (auto el : vec){
        os << el << std::endl;
    }
    return os;
}

/**
 * Heaviside step function: return 0 if arg < 0, 1 if arg >= 0
 */
template <typename T>
T heaviside(T arg){
    return arg < 0 ? 0. : 1.;
}

/**
 * Ramp function: R(x) = min(x, 0)
 */
 template <typename T>
 T ramp(T arg){
     return std::max(arg, static_cast<T>(0));
 }

#endif //PROGRAMM_UTILS_H