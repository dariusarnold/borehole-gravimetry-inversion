//
// Created by darius on 11/11/18.
//

#ifndef GRAVITYINVERSION_RESULT_H
#define GRAVITYINVERSION_RESULT_H


struct Result {
    double depth;
    double density;
};


/**
 * Output depth/density pair to standard stream
 * @param os
 * @param res
 * @return
 */
std::ostream& operator <<(std::ostream& os, const Result& res);


#endif //GRAVITYINVERSION_RESULT_H
