#include "Norms.h"
#include "utils.h"



L2ErrorNorm::L2ErrorNorm(const std::vector<double>& depth, const std::vector<double>& data, const std::vector<double>& errors) :
        ErrorNorm(depth, data, errors){}

L2ErrorNorm::~L2ErrorNorm() {}

double L2ErrorNorm::representant_function(double zj, double z) {
    // same as L2Norm
    return -gamma*heaviside(zj-z);
}

double L2ErrorNorm::gram_entry_analytical(double zj, double zk) {
    return gamma*gamma*std::min(zj, zk);
}
