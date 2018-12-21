#include "Norms.h"
#include "utils.h"



L2_Norm::L2_Norm(const std::vector<double> &depth, const std::vector<double> &data) :
    Norm(depth, data){}

double L2_Norm::gram_entry_analytical(double zj, double zk) {
    return gamma*gamma*std::min(zj, zk);
}

double L2_Norm::representant_function(double zj, double z){
    return -gamma*heaviside(zj-z);
}



W12_Norm::W12_Norm(const std::vector<double> &depth, const std::vector<double> &data) :
        Norm(depth, data){}

double W12_Norm::gram_entry_analytical(double zj, double zk) {
    double zmin = std::min(zj, zk);
    double gamma_square = gamma*gamma;
    return gamma_square*zj*zk + gamma_square*(1./3. * std::pow(zmin, 3) - 1./2. * std::pow(zmin, 2) * (zj+zk) + zmin*zj*zk);
}

double W12_Norm::representant_function(double zj, double z) {
    return gamma/2. * (zj - z) * (zj - z) * heaviside(zj-z) - gamma * (zj + 1./2. * zj*zj);
}


Seminorm::Seminorm(const std::vector<double>& depth, const std::vector<double>& data) :
    Norm(depth, data){}

double Seminorm::gram_entry_analytical(double zj, double zk) {
    double zmin = std::min(zj, zk);
    double gamma_square = gamma*gamma;
    return gamma_square*zj*zk + gamma_square*(1./3. * std::pow(zmin, 3) - 1./2. * std::pow(zmin, 2) * (zj+zk) + zmin*zj*zk);
}

double Seminorm::representant_function(double zj, double z) {
    return gamma/2. * (zj - z) * (zj - z) * heaviside(zj-z) - gamma * (zj + 1./2. * zj*zj);
}

void Seminorm::gram_matrix_analytical() {
    // get top left gram matrix, containing Gamma_jk
    Norm::gram_matrix_analytical();
    // calculate the additional data given as (gj, 1)
    Eigen::VectorXd additional(gram_matrix.cols()+1);
    for (size_t i = 0; i < measurement_depths.size(); ++i){
        additional(i) = -gamma*measurement_depths[i];
    }
    additional(additional.size()-1) = 0.;
    // add one more row and column to the matrix
    gram_matrix.conservativeResize(gram_matrix.rows()+1, gram_matrix.cols()+1);
    // fill additional space, Gram matrix is symmetrical so the same vector is inserted twice
    gram_matrix.rightCols(1) = additional;
    gram_matrix.bottomRows(1) = additional.transpose();
}

void Seminorm::solve_for_alpha() {
    // extend data by the additional constant 0
    measurement_data.emplace_back(0.);
    Norm::solve_for_alpha();
}

std::pair<Eigen::VectorXd, Eigen::VectorXd> Seminorm::calculate_density_distribution(uint64_t num_steps) {
    std::pair<Eigen::VectorXd, Eigen::VectorXd> depth_dens = Norm::calculate_density_distribution(num_steps);
    auto density_constant = alpha(alpha.rows()-1);
    depth_dens.second.array() += density_constant;
    return depth_dens;
}