#include <iostream>

#include <sstream>
#include <fstream>
#include <vector>


/**
 * Struct that holds measurement data read from file and associates the data
 * with a label.
 */
struct MeasurementData{
    double depth;
    double grav;

};

/**
 * Output value of measurement data to stream.
 * @param os
 * @param md
 * @return
 */
std::ostream& operator<<(std::ostream& os, const MeasurementData& md){
    os << md.depth << " " << md.grav << std::endl;
    return os;
}


/**
 * Read two values from a stringstream and save them in the Measurementdata.
 * The order of the values is the same as in the file.
 * @param input
 * @param md
 * @return
 */
std::istringstream& operator>>(std::istringstream& input, MeasurementData& md){
    input >> md.depth >> md.grav;
    return input;
}

class MeasurementReader{
public:
    std::vector<MeasurementData> read_measurements_file(std::string filepath);
};


std::vector<MeasurementData> MeasurementReader::read_measurements_file(std::string filepath) {
    std::ifstream matrixFile(filepath);
    std::string line;
    std::string num;
    std::vector<MeasurementData> output;
    if (matrixFile.is_open()){
        // read a line, ending on newline
        while (std::getline(matrixFile, line)) {
            MeasurementData row;
            // split and read into MeasurementData objects
            std::istringstream iss(line);
            iss >> row;
            //iss >> row.depth >> row.grav;
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