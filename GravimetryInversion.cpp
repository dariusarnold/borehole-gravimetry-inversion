//
// Class to perform a gravimetry inversion
//

#include <vector>
#include <fstream>
#include <iostream>

#include "MeasurementData.h"
#include "GravimetryInversion.h"
#include "utils.h"


void GravimetryInversion::read_measurements_file(const std::string& filepath) {
    std::ifstream matrixFile(filepath);
    std::vector<MeasurementData> output;
    if (matrixFile.is_open()){
        MeasurementData row;
        while (matrixFile >> row) {
            output.push_back(row);
        }
    }
    matrixFile.close();
    data = output;
}


void GravimetryInversion::print_data(){
    std::cout << data;
}