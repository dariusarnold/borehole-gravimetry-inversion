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



SemiErrorNorm::SemiErrorNorm(const std::vector<double> &depth, const std::vector<double> &data, const std::vector<double> &errors) :
    ErrorNorm(depth, data, errors){}


//SemiErrorNorm::~SemiErrorNorm() {}


double SemiErrorNorm::gram_entry_analytical(double zj, double zk) {
    //return Seminorm::gram_entry_analytical(zj, zk);
    // TODO this is just copied from SemiNorm, maybe there is a way to access that function here?
    double zmin = std::min(zj, zk);
    double gamma_square = gamma*gamma;
    return gamma_square*zj*zk + gamma_square*(1./3. * std::pow(zmin, 3) - 1./2. * std::pow(zmin, 2) * (zj+zk) + zmin*zj*zk);
}

double SemiErrorNorm::representant_function(double zj, double z) {
    // TODO this is just copied from SemiNorm, maybe there is a way to access that function here?
    return gamma/2. * (zj - z) * (zj - z) * heaviside(zj-z) - gamma * (zj + 1./2. * zj*zj);
}
