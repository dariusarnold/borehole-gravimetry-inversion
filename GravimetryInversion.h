#ifndef PROGRAMM_GRAVIMETRYINVERSION_H
#define PROGRAMM_GRAVIMETRYINVERSION_H


#include <experimental/filesystem>
#include <Eigen/Dense>

#include "Norms.h"


// set type for iterating over vectors from GravimetryInversion
typedef std::vector<double>::size_type vec_size_t;

namespace fs = std::experimental::filesystem;

/**
 * Class to perform a gravimetry inversion on borehole data.
 */
class GravimetryInversion{
public:
    /**
     * Constructor that takes a pointer to the norm which should be used.
     * Depending on the norm, the analytical solution for the Gram matrix and the way
     * to calculate the density distribution changes.
     * @param _norm Pointer to norm to use. Norm knows how to calculate Gram matrix and density
     * @param _discretization_steps number of discretization steps to use for density
     */
    explicit GravimetryInversion(std::unique_ptr<Norm> _norm, uint64_t _discretization_steps=10000);


    /**
     * Do an inversion using the norm type given as a template on data read from file and save the result,
     * a discretized density distribution in a file in the same path were data was read from.
     * @param filepath Path to file containing data in the following format:
     * no header, on column depth in meter, one col gravity measured in mGal, tab separated.
     * One depth/gravity pair per line, line ending \n
     * @param steps Number of discretization steps to use for density distribution
     */
    template <typename Norm_Type>
    static void invert_data_from_file(fs::path &filepath, uint64_t steps) {
        GravimetryInversion mr(std::unique_ptr<Norm_Type>(new Norm_Type), steps);
        mr.read_measurements_file(filepath);
        mr.calculate_gram_matrix();
        mr.solve_alpha();
        mr.calculate_density_distribution();
        filepath.replace_extension({".dens"});
        mr.write_density_distribution_to_file(filepath);
    }


    /**
     * Open .dat file and read measurements.
     * It is expected that data is given in depth, gravity measurement order.
     * @param filepath
     * @return Vector holding MeasurementData objects, which represent one row from the file
     */
    void read_measurements_file(const fs::path& filepath);

    /**
     * Calculate the gram matrix
     */
    void calculate_gram_matrix();

    /**
     * Solve the equation system given by d_j = Gamma_jk alpha for vector alpha
     */
     void solve_alpha();

     /**
      * Calculate density distribution from the representants and the alpha coefficients.
      */
     void calculate_density_distribution();

     /**
      * Discretize density distribution and save it to .txt file
      */
     void write_density_distribution_to_file(const fs::path& filepath);

private:
    std::unique_ptr<Norm> norm;
    uint64_t discretization_steps;          // discretization steps during integration
    std::vector<double> measurement_depths; // holds measurement depths read from file
    std::vector<double> measurement_data;   // holds measurement data (gravity acceleration) read from file
    std::vector<Result> result;             // holds depth/density distribution resulting from inversion
};


#endif //PROGRAMM_GRAVIMETRYINVERSION_H
