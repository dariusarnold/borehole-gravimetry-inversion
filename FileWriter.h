//
// Created by darius on 02/11/18.
//

#ifndef GRAVITYINVERSION_FILEWRITER_H
#define GRAVITYINVERSION_FILEWRITER_H


#include <fstream>

class FileWriter {
public:
    /**
     * It is expected that T is a list of values with one value per element
     * @tparam T
     * @param data
     * @return
     */
    template <typename T>
    void writeData(const T& depth, const T& data, const std::string& filepath){
        std::ofstream file;
        file.open(filepath);
        for (auto de = depth.begin(), da = data.begin(), e = depth.end(); de != e; ++de, ++da){
            file << *de << "\t" << *da << "\n";
        }
        file.close();
    }
};

#endif //GRAVITYINVERSION_FILEWRITER_H
