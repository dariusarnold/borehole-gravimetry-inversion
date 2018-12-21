#ifndef PROGRAMM_GRAVIMETRYINVERSION_H
#define PROGRAMM_GRAVIMETRYINVERSION_H


#include <experimental/filesystem>
#include <Eigen/Dense>
#include "FileIO.h"
#include "ModelParameters.h"


namespace fs = std::experimental::filesystem;

/**
 * Class to perform a gravimetry inversion on borehole data.
 */
 template <typename Norm_Type>
class GravimetryInversion{
public:
    /**
     * Constructor that takes a pointer to the norm which should be used.
     * Depending on the norm, the analytical solution for the Gram matrix and the way
     * to calculate the density distribution changes.
     * @param _norm Pointer to norm to use. Norm knows how to calculate Gram matrix and density
     * @param _discretization_steps number of discretization steps to use for density
     */
    explicit GravimetryInversion(std::unique_ptr<Norm_Type> _norm) :
        norm(std::move(_norm))
        {}

    /**
     * Do an inversion using the norm type given as a template on data read from file and save the result,
     * a discretized density distribution in a file in the same path were data was read from.
     * @param filepath Path to file containing data in the following format:
     * no header, on column depth in meter, one col gravity measured in mGal, tab separated.
     * One depth/gravity pair per line, line ending \n
     * @param steps Number of discretization steps to use for density distribution
     */
    static std::tuple<std::pair<Eigen::VectorXd, Eigen::VectorXd>, ModelParameters> invert_data_from_file(fs::path& filepath, uint64_t steps) {
        // read data from file
        auto [measurement_depths, measurement_data] = read_measurement_data_no_errors(filepath);
        auto _norm = std::make_unique<Norm_Type>(measurement_depths, measurement_data);
        GravimetryInversion gi(std::move(_norm));
        return gi.norm->do_work(steps);
    }

    static std::tuple<std::pair<Eigen::VectorXd, Eigen::VectorXd>, ModelParameters> invert_data_from_file_with_errors(fs::path& filepath, uint64_t steps, double nu=-1){
        // read data from file
        auto [measurement_depths, measurement_data, measurement_errors] = read_measurement_data_with_errors(filepath);
        // create a norm instance using this data
        auto _norm = std::make_unique<Norm_Type>(measurement_depths, measurement_data, measurement_errors);
        GravimetryInversion gi(std::move(_norm));
        if (nu < 0) {
            return gi.norm->do_work(steps);
        }else{
            return gi.norm->do_work(steps, nu);
        }
    }

    /**
     * Do an interpolation using an interpolating norm. Points are read from the file as x, y values
     * The interpolation function is discretized in the interval a, b and saved to a file.
     * @tparam Norm_Type Which norm to use
     * @param filepath Filepath input data
     * @param steps number of steps for discretization
     * @param a lower boundary of evaluation interval
     * @param b upper boundary of evaluation interval
     */
    static std::tuple<std::pair<Eigen::VectorXd, Eigen::VectorXd>, ModelParameters> interpolate_data_from_file(fs::path& filepath, uint64_t steps, double a, double b){
        auto [x_values, y_values] = read_function_data(filepath);
        GravimetryInversion gi(std::unique_ptr<Norm_Type>(new Norm_Type(a, b, x_values, y_values)));
        return gi.norm->do_work(steps);
    }


private:
    /**
     * Read file containing measurement data without errors
     */
    static std::tuple<std::vector<double>, std::vector<double>> read_measurement_data_no_errors(const fs::path& filepath){
        FileIO fr;
        return fr.readData(filepath);
    }

     /**
     * Read file containing measurement data with errors
     */
     static std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> read_measurement_data_with_errors(const fs::path& filepath){
         FileIO fr;
         return fr.readErrorData(filepath);
     }

     /**
      * Read file containing data for interpolation, one x, f(x) pair per line
      */
     static std::tuple<std::vector<double>, std::vector<double>> read_function_data(const fs::path& filepath){
         FileIO fr;
         return fr.readFunctionData(filepath);
     }

    std::unique_ptr<Norm_Type> norm;
};


#endif //PROGRAMM_GRAVIMETRYINVERSION_H
