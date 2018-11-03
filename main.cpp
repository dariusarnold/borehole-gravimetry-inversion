#include <iostream>

#include "MeasurementData.h"
#include "GravimetryInversion.h"


int main() {
    GravimetryInversion mr;
    mr.read_measurements_file("../4_grav3.dat");
    mr.print_data();
    mr.calculate_gram_matrix();
    mr.print_gram();
    mr.solve_alpha();
    mr.print_alpha();
    mr.write_density_distribution_to_file("../4_grav3_density.dat");
}