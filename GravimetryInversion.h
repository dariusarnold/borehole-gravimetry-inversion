#ifndef PROGRAMM_GRAVIMETRYINVERSION_H
#define PROGRAMM_GRAVIMETRYINVERSION_H


#include <functional>
#include <vector>
#include <experimental/filesystem>
#include <Eigen/Dense>


// set type for iterating over vectors from GravimetryInversion
typedef std::vector<double>::size_type vec_size_t;

namespace fs = std::experimental::filesystem;

/**
 * Class to perform a gravimetry inversion on borehole data.
 */
class GravimetryInversion{
public:
    GravimetryInversion(uint64_t discretization_steps=10000);

    /**
     * Do an inversion on data read from file and save the result, a discretized density distribution
     * in a file in the same path were data was read from.
     * @param filepath Path to file containing data in the following format:
     * no header, on column depth in meter, one col gravity measured in mGal, tab separated.
     * One depth/gravity pair per line, line ending \n
     */
    static void invert_data_from_file_L2_norm(fs::path &filepath, uint64_t steps);


    /**
     * Open .dat file and read measurements.
     * It is expected that data is given in depth, gravity measurement order.
     * @param filepath
     * @return Vector holding MeasurementData objects, which represent one row from the file
     */
    void read_measurements_file(const fs::path& filepath);

    /*
    * print a vector of printable elements to ostream
    */
    void print_data();

    /**
     * print gram matrix to cout
     */
     void print_gram();

     /**
      * Print coefficients of alpha to cout
      */
     void print_alpha();

    /**
     * Calculate the gram matrix Gamma_jk = integral_0^L g_j(z) * g_k(z) dz
     * with g_j, g_z Representants
     */
    void calculate_gram_matrix_L2_norm();

    /**
     * Solve the equation system given by d_j = Gamma_jk alpha for vector alpha
     */
     void solve_alpha();

     void calculate_density_distribution();

     /**
      * Discretize density distribution and save it to .txt file
      */
     void write_density_distribution_to_file(const fs::path& filepath);

private:
    const double LOWER_LIMIT = 0;   // m, lower limit of integral
    uint64_t discretization_steps;  // discretization steps during integration
    const double gamma = 0.08382;
    std::vector<MeasurementData> data;
    Eigen::MatrixXd gram_matrix;
    std::vector<double> alpha;
    std::vector<double> density;
    std::vector<double> depth_meters;
};


#endif //PROGRAMM_GRAVIMETRYINVERSION_H
