//
// Created by darius on 11/11/18.
//

#include "Representant.h"
#include "utils.h"


Representant_L2_Norm::Representant_L2_Norm(double _zj) : zj(_zj) {}


double Representant_L2_Norm::operator()(double z) const{
    return -gamma * heaviside(zj -z);
}
