#ifndef GRAVITYINVERSION_MODELPARAMETERS_H
#define GRAVITYINVERSION_MODELPARAMETERS_H


/*
 * Describes the parameters of a model resulting from the GravimetryInversion
 */
struct ModelParameters{
    // The lagrange multiplicator used to calculate the model
    double nu;
    // The misfit of the model given as Chi²/N
    double misfit;
    // The norm of the model, which characterizes the quality/fit of the model
    double norm;
};


#endif //GRAVITYINVERSION_MODELPARAMETERS_H
