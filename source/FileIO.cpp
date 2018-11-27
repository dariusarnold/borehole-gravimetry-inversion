//
// Created by darius on 02/11/18.
//

#include <fstream>
#include <FileIO.h>


#include "FileIO.h"
#include "Result.h"
#include "MeasurementData.h"

void FileIO::writeData(const std::vector<Result> &result, const fs::path &filepath){
    std::ofstream file;
    file.open(filepath);
    for (auto depth_value_pair : result){
        file << depth_value_pair;
    }
    file.close();
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> FileIO::readData(const fs::path &filepath) {
    std::vector<double> measurement_depths, measurement_data, measurement_error;
    std::ifstream file(filepath);
    if (file.is_open()){
        MeasurementData row;
        while (file >> row){
            measurement_depths.push_back(row.depth);
            measurement_data.push_back(row.grav);
            measurement_error.push_back(row.error);
        }
    }
    file.close();
    return std::make_tuple(measurement_depths, measurement_data, measurement_error);
}

std::pair<std::vector<double>, std::vector<double>> FileIO::readFunctionData(const fs::path &filepath) {
    std::vector<double> x, y;
    std::ifstream file(filepath);
    if (file.is_open()){
        Data row;
        while (file >> row){
            x.push_back(row.x);
            y.push_back(row.y);
        }
    }
    file.close();
    return std::make_pair(x, y);
}
