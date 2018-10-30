#ifndef PROGRAMM_GRAVIMETRYINVERSION_H
#define PROGRAMM_GRAVIMETRYINVERSION_H


class GravimetryInversion{
public:
    /**
     * Open .dat file and read measurements.
     * It is expected that data is given in depth, gravity measurement order.
     * @param filepath
     * @return Vector holding MeasurementData objects, which represent one row from the file
     */
    void read_measurements_file(const std::string& filepath);

    /*
    * print a vector of printable elements to ostream
    */
    void print_data();

private:
    std::vector<MeasurementData> data;
};




#endif //PROGRAMM_GRAVIMETRYINVERSION_H
