//
// Created by darius on 02/11/18.
//

#ifndef GRAVITYINVERSION_FILEWRITER_H
#define GRAVITYINVERSION_FILEWRITER_H


#include <experimental/filesystem>
#include "Result.h"

namespace fs = std::experimental::filesystem;


class FileWriter {
public:
     /**
      * Save depth/value pairs to file separated by a single whitespace.
      * Every pair is then separated by a newline.
      * No header or other data description is saved.
      * @param result Vector of Results structs, containing depth/density pairs
      * @param filepath
      */
    void writeData(const std::vector<Result>& result, const fs::path& filepath);
};

#endif //GRAVITYINVERSION_FILEWRITER_H
