//
// Created by darius on 02/11/18.
//

#include <fstream>

#include "FileIO.h"
#include "MeasurementData.h"

void FileIO::writeData(const std::vector<Result> &result, const fs::path &filepath){
    std::ofstream file;
    file.open(filepath);
    for (auto depth_value_pair : result){
        file << depth_value_pair;
    }
    file.close();
}

std::pair<std::vector<double>, std::vector<double>> FileIO::readData(const fs::path &filepath) {
    std::vector<double> measurement_depths, measurement_data;
    std::ifstream file(filepath);
    if (file.is_open()){
        MeasurementData row;
        while (file >> row){
            measurement_depths.push_back(row.depth);
            measurement_data.push_back(row.grav);
        }
    }
    file.close();
    return std::make_pair(measurement_depths, measurement_data);
}
