#include "ModelParameters.h"


std::ostream& operator<<(std::ostream& os, const ModelParameters& mp){
    os
    << "Nu: " << mp.nu << "\n"
    << "Misfit squared/N: " << mp.misfit << "\n"
    << "Norm: " << mp.norm;
    return os;
    }