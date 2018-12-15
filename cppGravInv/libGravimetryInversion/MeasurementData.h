//
// Struct that holds a single measurement data point
//

#ifndef PROGRAMM_MEASUREMENTDATA_H
#define PROGRAMM_MEASUREMENTDATA_H

/**
 * Struct that holds a x, f(x) function value pair
 */
struct FunctionData{
    double x;
    double y;
};


/**
 * Read two values from inputstream and save them in Data.
 * Order of values is expected to be x, f(x)
 * @param input
 * @param d
 * @return
 */
std::istream& operator>>(std::istream& input, FunctionData& d);

#endif //PROGRAMM_MEASUREMENTDATA_H
