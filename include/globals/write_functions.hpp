#pragma once

#include <string>
#include <vector>

template <typename T>
void write_vector_to_binary_file(const std::vector<T>& vec, size_t count, const std::string& filename, size_t paddingBytes = 0, char paddingValue = 0x00, bool append = false);