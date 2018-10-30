#include <iostream>

#include <sstream>
#include <fstream>
#include <vector>

#include "MeasurementData.h"
#include "GravimetryInversion.h"




int main() {
    GravimetryInversion mr;
    mr.read_measurements_file("../4_grav3.dat");
    mr.print_data();
}