//
// Created by darius on 02/11/18.
//

#ifndef GRAVITYINVERSION_FILEWRITER_H
#define GRAVITYINVERSION_FILEWRITER_H


#include <fstream>

class FileWriter {
public:
     /**
      * Save depth/value pairs to file separated by a single whitespace.
      * Every pair is then separated by a newline.
      * No header or other data description is saved.
      * @tparam T: T is an iterable  sequence of values
      * @param depth
      * @param data
      * @param filepath
      */
    template <typename T>
    void writeData(const T& depth, const T& data, const std::string& filepath){
        std::ofstream file;
        file.open(filepath);
        for (auto de = depth.begin(), da = data.begin(), e = depth.end(); de != e; ++de, ++da){
            file << *de << " " << *da << "\n";
        }
        file.close();
    }
};

#endif //GRAVITYINVERSION_FILEWRITER_H
