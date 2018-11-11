#include <iostream>
#include <experimental/filesystem>

#include "MeasurementData.h"
#include "GravimetryInversion.h"



namespace fs = std::experimental::filesystem;


/**
 * Call programm with positional arguments:
 * 1. input measurement data filepath
 * 2. number of discretization steps
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {
    if (argc > 1){
        // read command line parameters
        fs::path filepath{argv[1]};
        uint64_t steps = std::stoul(argv[2]);
        // input checking of steps, if less than 0 use default value
        steps > 0 ? GravimetryInversion::invert_data_from_file_L2_norm(filepath, steps) : GravimetryInversion::invert_data_from_file_L2_norm(filepath, 10000);
        return 0;
    }else {
        std::cout << "Call programm with positional arguments:\n"
                  << "1. input measurement data filepath\n"
                  << "2. number of discretization steps" << std::endl;
        return -1;
    }
}