#ifndef FNVHASH_H
#define FNVHASH_H

#include <iostream>
#include <fstream>
#include <string>


class FnvHash {
public:

  void add(const unsigned char* str, size_t size) {
    for (size_t c = 0; c < size; c++) {
      value ^= str[c];
      value *= prime;
    }
  }

  void write(const std::string& filename) {
    std::ofstream file(filename);
    file << value << std::endl;
    file.close();
  }

  void write(const std::string& filename, const std::string& suffix) {
    std::ofstream file(filename);
    file << value << " | " << suffix << std::endl;
    file.close();
  }

  unsigned long value = 2166136261u;

  friend std::ostream& operator << (std::ostream& os, const FnvHash& hash) {
    os << hash.value;
    return os;
  }

private:
  static const size_t prime = 16777619u;
};

#endif