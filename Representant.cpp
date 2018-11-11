//
// Created by darius on 11/11/18.
//

#include "Representant.h"
#include "utils.h"


Representant::Representant(double _zj) : zj(_zj) {}


double Representant::operator()(double z) const{
    return -gamma * heaviside(zj -z);
}
