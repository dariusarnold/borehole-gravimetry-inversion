//
// Created by darius on 11/11/18.
//

#ifndef GRAVITYINVERSION_REPRESENTANT_H
#define GRAVITYINVERSION_REPRESENTANT_H

#include <functional>


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
 * @param func1 Representant
 * @param func2 Representant
 * @return function objects used to evaluate the result of a*b
 */
template <typename Callable>
std::function<double(double)> operator*(const Callable& func1, const Callable& func2){
    return [&func1, &func2](double arg) { return func1(arg) * func2(arg);};
}


/**
 * Add two Representants, resulting in a new function used to evaluate the result
 * @param func1
 * @param func2
 * @return
 */
std::function<double(double)> operator+(const Representant& func1, const Representant& func2);


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


#endif //GRAVITYINVERSION_REPRESENTANT_H
