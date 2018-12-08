#include <iostream>
#include <experimental/filesystem>

#include "GravimetryInversion.h"
#include "Norms.h"


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
    if (argc != 4) {
        std::cout << "Call programm with positional arguments:\n"
                  << "1. input measurement data filepath\n"
                  << "2. number of discretization steps\n"
                  << "3. Id of norm: L2 (0)\n"
                  << std::endl;
        return -1;
    }
    else{
        // read command line parameters
        fs::path filepath{argv[1]};
        uint64_t steps = std::stoul(argv[2]);
        switch (atoi(argv[3])){
            case 0:
                GravimetryInversion<L2ErrorNorm>::invert_data_from_file_with_errors(filepath, steps);
                break;
            case 1:
                GravimetryInversion<SemiErrorNorm>::invert_data_from_file_with_errors(filepath, steps);
                break;
            default:
                std::cout << "Enter valid norm id" << std::endl;
                return -1;
        }
        return 0;
    }
}