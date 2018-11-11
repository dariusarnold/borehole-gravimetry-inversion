//
// Created by darius on 02/11/18.
//

#include <fstream>

#include "FileWriter.h"

void FileWriter::writeData(const std::vector<Result> &result, const fs::path &filepath){
    std::ofstream file;
    file.open(filepath);
    for (auto depth_value_pair : result){
        file << depth_value_pair;
    }
    file.close();
}