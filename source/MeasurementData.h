//
// Struct that holds a single measurement data point
//

#ifndef PROGRAMM_MEASUREMENTDATA_H
#define PROGRAMM_MEASUREMENTDATA_H


const double FREE_AIR_GRADIENT = 0.308; // MGal/m


/**
 * Struct that holds measurement data with measurement errors read from file
 * and associates the data with a label which describes the data.
 * depth unit is Meters, with a positive Z axis down the borehole.
 * grav unit is mGal.
 * error unit is the same as grav.
 */
struct MeasurementData{
    double depth;
    double grav;
    double error;
};


/**
 * Output values of measurement data to stream.
 * @param os
 * @param md
 * @return
 */
std::ostream& operator<<(std::ostream& os, const MeasurementData& md);


/**
 * Read two values from a inputstream and save them in the Measurementdata.
 * Corrects for free air gradient
 * The order of the values is expected to be depth, gravity, error
 * @param input
 * @param md
 * @return
 */
std::istream& operator>>(std::istream& input, MeasurementData& md);


/**
 * Struct that holds a x, f(x) function value pair
 */
struct Data{
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
std::istream& operator>>(std::istream& input, Data& d);

#endif //PROGRAMM_MEASUREMENTDATA_H
