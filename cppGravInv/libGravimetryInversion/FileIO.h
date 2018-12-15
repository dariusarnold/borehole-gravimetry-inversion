//
// Created by darius on 02/11/18.
//

#ifndef GRAVITYINVERSION_FILEWRITER_H
#define GRAVITYINVERSION_FILEWRITER_H


#include <experimental/filesystem>
#include <fstream>
#include <iterator>

// forward declaration
struct Result;

namespace fs = std::experimental::filesystem;


class FileIO {
public:
     /**
      * Save depth/value pairs to file separated by a single whitespace.
      * Every pair is then separated by a newline.
      * No header or other data description is saved.
      * @param result Vector of Results structs, containing depth/density pairs
      * @param filepath full path with filename in which file is created
      */
    void writeData(const std::vector<Result>& result, const fs::path& filepath);

    /**
     * Write the discretized version of the interpolated function into a file
     * @param x all values of x
     * @param f_of_x the corresponding values f(x)
     * @param filepath full path with filename in which file is created
     */
    void writeData(const std::vector<double>& x, const std::vector<double>& f_of_x, const fs::path& filepath);

    /**
     * Read depth, value pairs from file and return them in two vectors
     * @param filepath full path with filename
     * @return first vector contains all depth values, second vector all measurement values
     */
    std::tuple<std::vector<double>, std::vector<double>> readData(const fs::path& filepath);

    /**
     * Read depth, value, error pairs from file and return them in three vectors
     * @param filepath full path with filename
     * @return first vector contains all depth values, second vector all measurement values, third all error values
     */
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> readErrorData(const fs::path& filepath);

    /**
     * Read x, f(x) pairs from file and return them in two vectors
     * @param filepath full path with filename
     * @return first vector contains all x values, second vector contains all f(x) values
     */
    std::pair<std::vector<double>, std::vector<double>> readFunctionData(const fs::path& filepath);

private:
    const double FREE_AIR_GRADIENT = 0.308; // MGal/m
};

#endif //GRAVITYINVERSION_FILEWRITER_H
