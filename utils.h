//
// Small utility functions or classes
//

#ifndef PROGRAMM_UTILS_H
#define PROGRAMM_UTILS_H

/*
 * print a vector of printable elements to ostream, element for element
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> vec){
    for (auto el : vec){
        os << el;
    }
    return os;
}

#endif //PROGRAMM_UTILS_H