#include <iostream>

#include <sstream>
#include <fstream>
#include <vector>

#include "MeasurementData.h"


class MeasurementReader{
public:
    std::vector<MeasurementData> read_measurements_file(const std::string& filepath);
};


/**
 * Open .dat file and read measurements.
 * It is expected that data is given in depth, gravity measurement order.
 * @param filepath
 * @return Vector holding MeasurementData objects, which represent one row from the file
 */
std::vector<MeasurementData> MeasurementReader::read_measurements_file(const std::string& filepath) {
    std::ifstream matrixFile(filepath);
    std::vector<MeasurementData> output;
    if (matrixFile.is_open()){
        MeasurementData row;
        while (matrixFile >> row) {
            output.push_back(row);
        }
    }
    matrixFile.close();
    return output;
}


/*
 * print a vector of printable elements to ostream
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> vec){
    for (auto el : vec){
        os << el;
    }
    return os;
}


int main() {
    MeasurementReader mr;
    auto data = mr.read_measurements_file("../4_grav3.dat");
    std::cout << data;
    return 1;
}