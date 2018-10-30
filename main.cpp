#include <iostream>

#include <sstream>
#include <fstream>
#include <vector>


/**
 * Struct that holds measurement data read from file and associates the data
 * with a label which describes the data.
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
std::ostream& operator<<(std::ostream& os, const MeasurementData& md){
    os << md.depth << " " << md.grav << std::endl;
    return os;
}


/**
 * Read two values from a inputstream and save them in the Measurementdata.
 * The order of the values is expected to be depth, gravity.
 * @param input
 * @param md
 * @return
 */
std::istream& operator>>(std::istream& input, MeasurementData& md){
    input >> md.depth >> md.grav;
    return input;
}

class MeasurementReader{
public:
    std::vector<MeasurementData> read_measurements_file(const std::string& filepath);
};


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