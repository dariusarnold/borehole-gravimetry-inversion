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
    if (argc != 5) {
        std::cout << "Call programm with positional arguments:\n"
                  << "1. input measurement data filepath\n"
                  << "2. number of discretization steps\n"
                  << "3. Id of norm: L2 (0), W12 (1), Semi (2)\n"
                  << "4. Lagrange multiplicator nu"
                  << std::endl;
        return -1;
    }
    else{
        // read command line parameters
        fs::path filepath{argv[1]};
        uint64_t steps = std::stoul(argv[2]);
        double nu = atof(argv[4]);
        switch (atoi(argv[3])){
            case 0:
                GravimetryInversion::invert_data_from_file_with_errors<L2ErrorNorm>(filepath, steps, nu);
                break;
            default:
                std::cout << "Enter valid norm id" << std::endl;
                return -1;
        }
        return 0;
    }
}