#include <iostream>

#include "MeasurementData.h"
#include "GravimetryInversion.h"



void gravimetry_inversion(const std::string& filepath){
    GravimetryInversion mr;
    mr.read_measurements_file(filepath);
    //mr.print_data();
    mr.calculate_gram_matrix();
    //mr.print_gram();
    mr.solve_alpha();
    //mr.print_alpha();
    // find file ending
    auto s = filepath.rfind(".dat");
    // create new string since filepath is const, append _density to file ending
    std::string out_fname = filepath;
    out_fname.replace(s, 4, std::string("_density.dat"));
    mr.write_density_distribution_to_file(out_fname);
}


int main(int argc, char* argv[]) {
    std::string filepath = argv[1];
    std::cout << filepath << std::endl;
    gravimetry_inversion(filepath);
}