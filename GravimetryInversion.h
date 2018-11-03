#ifndef PROGRAMM_GRAVIMETRYINVERSION_H
#define PROGRAMM_GRAVIMETRYINVERSION_H


#include <functional>
#include <vector>

#include <Eigen/Dense>

/**
 * Representant function with the form of a step function
 */
struct Representant{
    explicit Representant(double zj);
    double zj;      // depth of step in meters
    const double gamma = 0.08382; //mGal * cm³ / (m g)
    double operator()(double z) const;  // calculate value of Representant at depth z
};

/**
 * Multiply two Representants a and b, returns a function that evaluates the multiplication.
 * Eg. a * b results in a new function, that can be used to evaluate a * b at any point
 * @param a Representant
 * @param b Representant
 * @return function objects used to evaluate the result of a*b
 */
 template <typename Callable>
std::function<double(double)> operator*(const Callable& a, const Callable& b){
    return [&a, &b](double arg) { return a(arg) * b(arg);};
}

/**
 * Calculate integral of f from the lower to the upper limit, using the given amount of steps
 * for discretization.
 * @tparam Func Callable type that overrides operator(). The given argument is the
 * position at which the function should be evaluated.
 * @param f Function to integrate.
 * @param lower_limit lower integral limit
 * @param upper_limit upper integral limit
 * @param steps amount of steps for discretisation
 * @return
 */
struct Integrator{
    double operator()(const std::function<double(double)>& f, double lower_limit, double upper_limit, uint32_t steps){
        double result = 0.;
        double stepsize = (upper_limit - lower_limit) / steps;
        if (stepsize <= 0) return 0.;    // upper limit below lower limit
        for (uint32_t i = 0; i < steps; ++i){
            result += f(lower_limit + i * stepsize) * stepsize;
        }
        return result;
    }
};


/**
 * Class to perform a gravimetry inversion on borehole data.
 */
class GravimetryInversion{
public:
    /**
     * Open .dat file and read measurements.
     * It is expected that data is given in depth, gravity measurement order.
     * @param filepath
     * @return Vector holding MeasurementData objects, which represent one row from the file
     */
    void read_measurements_file(const std::string& filepath);

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
    void calculate_gram_matrix();

    /**
     * Solve the equation system given by d_j = Gamma_jk alpha for vector alpha
     */
     void solve_alpha();

     /**
      * Discretize density distribution and save it to .txt file
      */
     void write_density_distribution_to_file(const std::string& filepath);

private:
    const double LOWER_LIMIT = 0; // m, lower limit of integral
    const uint32_t INTEGRAL_STEPS = 10000;  // discretization steps during integration
    std::vector<MeasurementData> data;
    std::vector<Representant> representant_functions;
    Eigen::MatrixXd gram_matrix;
    std::vector<double> alpha;

};


#endif //PROGRAMM_GRAVIMETRYINVERSION_H
