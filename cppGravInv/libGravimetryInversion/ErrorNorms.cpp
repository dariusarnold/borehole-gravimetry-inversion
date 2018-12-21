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

void SemiErrorNorm::gram_matrix_analytical(double nu) {
    Eigen::MatrixXd sigma_squared = sigma_matrix*sigma_matrix;
    // get top left gram matrix, containing Gamma_jk
    Norm::gram_matrix_analytical();
    // modify gram matrix by adding nu^-1*sigma²
    gram_matrix += 1/nu * sigma_squared;
    // calculate the additional data given as (gj, 1)
    Eigen::VectorXd additional = Eigen::VectorXd::Zero(gram_matrix.cols()+1);
    for (size_t i = 0; i < measurement_depths.size(); ++i){
        additional(i) = -gamma*measurement_depths[i];
    }
    //additional(additional.size()-1) = 0.;
    // add one more row and column to the matrix
    gram_matrix.conservativeResize(gram_matrix.rows()+1, gram_matrix.cols()+1);
    // fill additional space, Gram matrix is symmetrical so the same vector is inserted twice
    gram_matrix.rightCols(1) = additional;
    gram_matrix.bottomRows(1) = additional.transpose();
}

void SemiErrorNorm::solve_for_alpha(double nu) {
    // calc gram matrix and extend it by H = (g_j, h_k)
    gram_matrix_analytical(nu);
    // now set up the equation to be solved to calculate alpha
    //Eigen::MatrixXd term = 1/nu * sigma_squared + gram_matrix;
    auto data = measurement_data;
    data.push_back(0);
    Eigen::VectorXd data_vec = std_to_eigen(data);
    alpha = gram_matrix.colPivHouseholderQr().solve(data_vec);
}

double SemiErrorNorm::calculate_misfit(double nu) {
    double betrag = (1/nu * sigma_matrix * alpha.head(alpha.rows()-1)).norm();
    return betrag*betrag;
}

std::pair<Eigen::VectorXd, Eigen::VectorXd> SemiErrorNorm::calculate_density_distribution(uint64_t num_steps) {
    // TODO this is the same as in SemiNorm, maybe we can derive from that as well an inherit it?
    std::pair<Eigen::VectorXd, Eigen::VectorXd> depth_dens = ErrorNorm::calculate_density_distribution(num_steps);
    auto density_constant = alpha(alpha.rows()-1);
    depth_dens.second.array() += density_constant;
    return depth_dens;
}
