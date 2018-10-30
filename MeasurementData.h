#ifndef PROGRAMM_MEASUREMENTDATA_H
#define PROGRAMM_MEASUREMENTDATA_H

const double FREE_AIR_GRADIENT = 0.308; // MGal/m

/**
 * Struct that holds measurement data read from file and associates the data
 * with a label which describes the data.
 * depth unit is Meters, with a positive Z axis down the borehole.
 * grav unit is mGal.
 */
struct MeasurementData{
    double depth;
    double grav;
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
 * The order of the values is expected to be depth, gravity.
 * @param input
 * @param md
 * @return
 */
std::istream& operator>>(std::istream& input, MeasurementData& md);

#endif //PROGRAMM_MEASUREMENTDATA_H
