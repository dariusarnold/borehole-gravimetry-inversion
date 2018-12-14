//
// Created by darius on 02/11/18.
//

#include <fstream>
#include <iterator>

#include "FileIO.h"
#include "Result.h"
#include "MeasurementData.h"
#include <stdexcept>

void FileIO::writeData(const std::vector<Result> &result, const fs::path &filepath){
    std::ofstream file;
    file.open(filepath);
    for (auto depth_value_pair : result){
        file << depth_value_pair;
    }
    file.close();
}

std::tuple<std::vector<double>, std::vector<double>> FileIO::readData(const fs::path &filepath) {
    std::vector<double> measurement_depths, measurement_data;
    double depth, grav;
    std::ifstream file(filepath);
    std::string line;
    if (file.is_open()){
        // only depth and data
        while (std::getline(file, line)) {
            std::istringstream sstream(line);
            sstream >> depth;
            sstream >> grav;
            measurement_depths.push_back(depth);
            measurement_data.push_back(grav-FREE_AIR_GRADIENT*depth);
        }
    }else{
        std::string what = "Cant open file: " + filepath.string();
        throw std::runtime_error(what);
    }
    file.close();
    return std::tuple(measurement_depths, measurement_data);
}



std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> FileIO::readErrorData(const fs::path &filepath) {
    // depth, data, error
    double depth, grav, error;
    std::vector<double> measurement_depths, measurement_data, measurement_error;
    std::ifstream file(filepath);
    std::string line;
    if (file.is_open()){
        while (std::getline(file, line)) {
            std::istringstream sstream(line);
            sstream >> depth;
            sstream >> grav;
            sstream >> error;
            measurement_depths.push_back(depth);
            measurement_data.push_back(grav-FREE_AIR_GRADIENT*depth);
            measurement_error.push_back(error);
            }
    }else{
        std::string what = "Cant open file: " + filepath.string();
        throw std::runtime_error(what);
        }
    file.close();
    return std::tuple(measurement_depths, measurement_data, measurement_error);
}

std::pair<std::vector<double>, std::vector<double>> FileIO::readFunctionData(const fs::path &filepath) {
    std::vector<double> x, y;
    std::ifstream file(filepath);
    if (file.is_open()){
        FunctionData row;
        while (file >> row){
            x.push_back(row.x);
            y.push_back(row.y);
        }
    }
    file.close();
    return std::make_pair(x, y);
}
