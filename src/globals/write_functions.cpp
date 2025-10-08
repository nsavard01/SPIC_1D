#include "globals/write_functions.hpp"
#include <cstring>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <dirent.h>
#include <cmath>

template <typename T>
void write_vector_to_binary_file(const std::vector<T>&  data, size_t count, const std::string& filename, size_t paddingBytes, char paddingValue, bool append) {
    size_t size_of_T = sizeof(T);
    size_t total_size = 2 * paddingBytes + count * size_of_T;

    std::vector<char> byte_stream(total_size);

    if (paddingBytes > 0) {
        std::fill(byte_stream.begin(), byte_stream.begin() + paddingBytes, paddingValue);
    }

    std::memcpy(byte_stream.data() + paddingBytes, data.data(), count * size_of_T);

    if (paddingBytes > 0) {
        std::fill(byte_stream.end() - paddingBytes, byte_stream.end(), paddingValue);
    }

    std::ios_base::openmode mode = std::ios::binary | (append ? std::ios::app : std::ios::trunc);
    std::ofstream outFile(filename, mode);
    if (!outFile) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    outFile.write(byte_stream.data(), byte_stream.size());
    outFile.close();
}




template void write_vector_to_binary_file<double>(const std::vector<double>&, size_t, const std::string&, size_t, char, bool);
template void write_vector_to_binary_file<int>(const std::vector<int>&, size_t, const std::string&, size_t, char, bool);
